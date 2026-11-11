#!/usr/bin/env python3
"""Host regressions: Fortran CSR helpers, actual CUDA ops via serial emulation,
and the MPI solver driver using deterministic CPU substitutes for cuSPARSE.
No GPU performance or CUDA runtime correctness is established by these tests.
"""
from pathlib import Path
import argparse
import os
import json
import re
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
APP = HERE.parent
NTPOLY = APP.parent.parent


def run(args, **kwargs):
    result = subprocess.run([str(a) for a in args], text=True, capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f"{args}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ntpoly-build', type=Path, default=NTPOLY / 'BuildCheck',
                        help='CPU NTPoly build containing include/ and lib/libNTPoly.a')
    options = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="ntg-host-checks-") as directory:
        work = Path(directory)
        # Exercise the real CMake custom command with a fake nvcc. No CUDA
        # compilation happens here: this checks only dependency/rebuild behavior.
        fixture = work / "cmake-fixture"
        fixture.mkdir()
        fixture_app = fixture / "app"
        fixture_app.mkdir()
        (fixture_app / 'tests').mkdir()
        shutil.copyfile(HERE / 'device_partials_test.cpp', fixture_app / 'tests/device_partials_test.cpp')
        shutil.copyfile(HERE / 'resident_panels_test.cpp', fixture_app / 'tests/resident_panels_test.cpp')
        shutil.copyfile(HERE / 'nccl_transport_test.cpp', fixture_app / 'tests/nccl_transport_test.cpp')
        shutil.copyfile(HERE / 'test_identity.f90', fixture_app / 'tests/test_identity.f90')
        for name in ("CMakeLists.txt", "yukawa_gemm.f90", "trs2_reference.F90", "solver_csr.f90", "solver_nccl.f90", "generate_yukawa.f90", "cusparse_gemm.cu", "nccl_panels.cu", "resident_panels.cu"):
            shutil.copyfile(APP / name, fixture_app / name)
        (fixture / "CMakeLists.txt").write_text("cmake_minimum_required(VERSION 3.20)\n"
            "project(DependencyCheck LANGUAGES CXX Fortran)\nadd_library(NTPoly INTERFACE)\nadd_subdirectory(app)\n")
        fake_cuda = fixture / "cuda"
        (fake_cuda / "bin").mkdir(parents=True)
        (fake_cuda / "lib64").mkdir()
        fake_nvcc = fake_cuda / "bin/nvcc"
        fake_nvcc.write_text("#!/usr/bin/env python3\nimport pathlib,sys\n"
            "out=pathlib.Path(sys.argv[sys.argv.index('-o')+1])\n"
            "out.write_bytes(b'fake object')\n"
            "with out.with_suffix('.calls').open('a') as f: f.write('compile\\n')\n")
        fake_nvcc.chmod(0o755)
        fake_lib = fake_cuda / "lib64/libcusparse.so"
        fake_lib.touch()
        cmake_build = fixture / "build"
        run(["cmake", "-S", fixture, "-B", cmake_build, "-G", "Ninja",
             f"-DNVCC_BIN={fake_nvcc}", f"-DCUSPARSE_LIB={fake_lib}"])
        assert not (cmake_build / "app/cusparse_gemm.o").exists()
        run(["ninja", "-C", cmake_build, "app/cusparse_gemm.o"])
        with (fixture_app / "cusparse_gemm.cu").open("a") as f:
            f.write("\n// exercise source dependency\n")
        run(["ninja", "-C", cmake_build, "app/cusparse_gemm.o"])
        assert (cmake_build / "app/cusparse_gemm.calls").read_text().count("compile") == 2
        print("CUDA object rebuilds when its source changes (mock nvcc): PASS")
        run(["gfortran", "-g", "-fcheck=all", "-J", work, "-I", work,
             APP / "solver_csr.f90", HERE / "test_solver_csr.f90", "-o", work / "csr"])
        print(run([work / "csr"]).strip())
        run(["gfortran", "-g", "-fcheck=all", "-ftrapv", "-J", work,
             NTPOLY / "Source/Fortran/GPUPlanningModule.F90", HERE / "test_gpu_planning.f90",
             "-o", work / "gpu-planning"])
        print(run([work / "gpu-planning"]).strip())
        source = (APP / "cusparse_gemm.cu").read_text()
        ops = source[source.index("#define OPS_SLOTS"):source.index('} // extern "C"')]
        # Only execution is substituted: compile the actual kernel and operation bodies.
        ops = re.sub(r"(\w+)<<<([^,]+),\s*([^>]+)>>>(\(.*?\));",
                     lambda m: f"launch({m[2]}, {m[3]}, [&]() {{ {m[1]}{m[4]}; }});",
                     ops, flags=re.S)
        emulated = '#include "host_cuda_shim.hpp"\nextern "C" {\n' + ops + '\n}\n'
        (work / "ops.cpp").write_text(emulated)
        (work / "ops_test.cpp").write_text(emulated + (HERE / "test_gpu_ops.cpp").read_text())
        run(["clang++", "-std=c++17", "-fsanitize=address,undefined", "-g", "-I", HERE,
             work / "ops_test.cpp", "-o", work / "ops-test"])
        print(run([work / "ops-test"]).strip())
        cache = source[source.index("// Persistent A_local"):source.index("// Persistent threshold-pipeline")]
        kernel = source[source.index("__global__ void subtract_one_kernel"):source.index("__global__ void threshold_mark_kernel")]
        helpers = source[source.index("static int ensure_aloc_uploaded"):source.index("// Ensure persistent C row_ptr")]
        helpers = re.sub(r"(\w+)<<<([^,]+),\s*([^>]+)>>>(\(.*?\));",
                         lambda m: f"launch({m[2]}, {m[3]}, [&]() {{ {m[1]}{m[4]}; }});",
                         helpers, flags=re.S)
        (work / "cache.cpp").write_text('#include "host_cuda_shim.hpp"\n#include "test_upload_cache.hpp"\n' +
                                       cache + kernel + helpers + (HERE / "test_upload_cache.cpp").read_text())
        run(["clang++", "-std=c++17", "-fsanitize=address,undefined", "-g", "-I", HERE,
             work / "cache.cpp", "-o", work / "cache-test"])
        print(run([work / "cache-test"]).strip())
        # Link the real Fortran driver to the existing local NTPoly CPU build.
        build = options.ntpoly_build.resolve()
        if not (build / "lib/libNTPoly.a").exists():
            raise RuntimeError("Build NTPoly's CPU library in BuildCheck before MPI driver checks")
        run(["clang++", "-std=c++17", "-g", "-I", HERE, "-c", work / "ops.cpp", "-o", work / "ops.o"])
        run(["clang++", "-std=c++17", "-g", "-c", HERE / "host_spgemm.cpp", "-o", work / "spgemm.o"])
        run(["clang++", "-std=c++17", "-x", "c++", "-c", APP / "nccl_panels.cu", "-o", work / "nccl-stub.o"])
        # This module makes the driver's device selection a no-op in host tests.
        (work / "openacc.f90").write_text("""MODULE openacc
  INTEGER, PARAMETER :: ACC_DEVICE_NVIDIA=1
CONTAINS
  SUBROUTINE ACC_SET_DEVICE_NUM(device,kind)
    INTEGER :: device,kind
  END SUBROUTINE
END MODULE
""")
        run(["gfortran", "-c", "-J", work, work / "openacc.f90", "-o", work / "acc.o"])
        run(["mpifort", "-g", "-fcheck=all", "-ffree-line-length-none", "-I", work,
             "-I", build / "include", "-J", work, "-c", APP / "solver_csr.f90", "-o", work / "csr.o"])
        run(["mpifort", "-g", "-fcheck=all", "-ffree-line-length-none", "-J", work,
             "-c", APP / "solver_nccl.f90", "-o", work / "nccl.o"])
        run(["mpifort", "-g", "-fcheck=all", "-ffree-line-length-none", "-I", work,
             "-I", build / "include", "-c", APP / "yukawa_gemm.f90", "-o", work / "driver.o"])
        run(["mpifort", work / "driver.o", work / "csr.o", work / "nccl.o", work / "nccl-stub.o", work / "acc.o", work / "ops.o",
             work / "spgemm.o", build / "lib/libNTPoly.a", "-lc++", "-fopenmp",
             "-framework", "Accelerate", "-o", work / "solver"])
        # Odd dimension exercises padded NTPoly blocks; off-diagonal blocks are nonsymmetric.
        h = [[-1.1, .03, .02, .04, .01], [.03, -.8, .05, .01, .02],
             [.02, .05, .7, .03, .01], [.04, .01, .03, 1.0, .02],
             [.01, .02, .01, .02, 1.2]]
        matrix = work / "h.mtx"
        matrix.write_text("%%MatrixMarket matrix coordinate real general\n5 5 25\n" +
                          "".join(f"{r+1} {c+1} {h[r][c]:.17g}\n" for r in range(5) for c in range(5)))
        env = dict(os.environ, OMP_NUM_THREADS="1")
        run(["mpifort", "-g", "-fcheck=all", "-I", build / "include", "-J", work,
             HERE / "test_identity.f90", build / "lib/libNTPoly.a", "-fopenmp",
             "-framework", "Accelerate", "-o", work / "identity"])
        for ranks in (1, 2, 4):
            print(run(["mpirun", "-n", ranks, work / "identity"], env=env, timeout=30).strip())
        for ranks, side, balance, diagonal in [(1,1,False,False), (2,1,False,False), (4,2,False,False), (4,2,True,False), (4,2,False,True)]:
            case = work / f"np{ranks}-lb{balance}-diag{diagonal}"
            case.mkdir()
            case_matrix = matrix
            if diagonal:
                case_matrix = case / "diagonal.mtx"
                case_matrix.write_text("%%MatrixMarket matrix coordinate real general\n5 5 5\n" +
                                       "".join(f"{i+1} {i+1} {h[i][i]}\n" for i in range(5)))
            command = ["mpirun", "-n", ranks, work / "solver", "--input", case_matrix,
                       "--mode", "fixed", "--solver", "trs2", "--nelec", "2", "--max-iter", "60",
                       "--conv-diff", "1e-10", "--threshold", "1e-8", "--prows", side,
                       "--pcols", ranks//side, "--pslices", "1", "--chunk-work", "30"]
            if balance:
                command.append("--load-balance")
            for resident in (False, True, "dense"):
                output = run(command + (["--gpu-dense"] if resident=="dense" else ["--device-partials"] if resident else []), env=env, cwd=case, timeout=60)
                if "gpu energy" not in output:
                    raise RuntimeError(output)
                # NTPoly logs its iteration count; enforce parity in these deterministic cases.
                cpu_iterations = re.search(r"Total Iterations:\s*(\d+)", output)
                gpu_iterations = re.search(r"iterations\s*:\s*(\d+)", output)
                if cpu_iterations is None or gpu_iterations is None:
                    raise RuntimeError("Missing iteration counts:\n"+output)
                assert cpu_iterations[1] == gpu_iterations[1], output
                print(f"MPI TRS2 np={ranks}, load_balance={balance}, diagonal={diagonal}, device_partials={resident}: PASS")
        repeated = run(command + ['--solver-reps', '3'], env=env, cwd=case, timeout=60)
        expected = re.search(r'Total Iterations:\s*(\d+)', repeated)[1]
        assert re.findall(r'iterations\s*:\s*(\d+)', repeated) == [expected]*3, repeated
        assert 'gpu_chain_with_setup_seconds:' in repeated
        print('Three repeated solves restore X0 and retain matrix/iteration parity: PASS')
        failure = subprocess.run([str(a) for a in command], text=True, capture_output=True,
                                 env=dict(env, NTG_TEST_FAIL_SPGEMM="1"), cwd=case, timeout=30)
        assert failure.returncode != 0 and "partial SpGEMM failed" in failure.stdout + failure.stderr
        print("Injected partial-product failure aborts the solver: PASS")
        for op in ("scale", "increment", "copy", "trace", "dot", "norm", "sigma", "grandsum"):
            summary = work / f"{op}.json"
            run(["mpirun", "-n", "4", work / "solver", "--input", matrix, "--mode", "fixed",
                 "--op", op, "--reps", "3", "--prows", "2", "--pcols", "2", "--pslices", "1",
                 "--out", summary], env=env, cwd=work, timeout=30)
            result = json.loads(summary.read_text())
            assert result["rel_err"] < 1e-12, result
        print("P1/P2 operation benchmarks on four MPI ranks: PASS")


if __name__ == "__main__":
    main()
