#pragma once
#include <cuda_runtime.h>
#include <thrust/system/cuda/execution_policy.h>
#include <algorithm>
#include <cstdint>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>

// Weak lookup preserves standalone CUDA primitives without an NCCL context.
extern "C" {
int solver_nccl_workspace_get_kind_(int,void**) __attribute__((weak));
int solver_nccl_workspace_attach_kind_(int,void*,void(*)(void*),uint64_t(*)(void*)) __attribute__((weak));
int solver_nccl_workspace_clear_kind_(int) __attribute__((weak));
}
namespace ntpoly_gpu {
constexpr int real_workspace=1,panel_workspace=2;
inline uint64_t allocation_count=0;
struct Buffer {
    void* data=nullptr;
    size_t capacity=0;
    Buffer()=default;
    Buffer(const Buffer&)=delete;
    Buffer& operator=(const Buffer&)=delete;
    ~Buffer(){if(data)cudaFree(data);}
    int reserve(size_t bytes) {
        bytes=std::max<size_t>(bytes,8);
        if(bytes<=capacity)return 0;
        if(data){if(cudaFree(data)!=cudaSuccess)return -1;data=nullptr;capacity=0;}
        if(cudaMalloc(&data,bytes)!=cudaSuccess)return -1;
        capacity=bytes;++allocation_count;return 0;
    }
    template<class T>T* as(){return static_cast<T*>(data);}
};
struct ScratchPool {
    struct Block {std::unique_ptr<Buffer> buffer;bool used=false;};
    std::vector<Block> blocks;
    char* acquire(size_t bytes) {
        Block* selected=nullptr;
        for(auto& b:blocks)if(!b.used && b.buffer->capacity>=bytes)
            if(!selected || b.buffer->capacity<selected->buffer->capacity)selected=&b;
        if(!selected){
            for(auto& b:blocks)if(!b.used){selected=&b;break;}
            if(!selected){blocks.push_back({std::make_unique<Buffer>(),false});selected=&blocks.back();}
        }
        if(selected->buffer->reserve(bytes)!=0)throw std::bad_alloc();
        selected->used=true;return selected->buffer->as<char>();
    }
    void release(char* pointer) {
        for(auto& b:blocks)if(b.buffer->data==pointer && b.used){b.used=false;return;}
        throw std::runtime_error("invalid GPU temporary release");
    }
    uint64_t bytes()const {
        uint64_t result=0;for(const auto& b:blocks)result+=b.buffer->capacity;return result;
    }
};
struct Allocator {
    using value_type=char;
    ScratchPool* pool;
    char* allocate(std::ptrdiff_t bytes){if(bytes<0)throw std::bad_alloc();return pool->acquire(size_t(bytes));}
    void deallocate(char* pointer,size_t){pool->release(pointer);}
};
template<class T>T* workspace(int kind) {
    if(!solver_nccl_workspace_get_kind_ || !solver_nccl_workspace_attach_kind_)return nullptr;
    void* pointer=nullptr;
    if(solver_nccl_workspace_get_kind_(kind,&pointer)!=0)return nullptr;
    if(!pointer){
        auto owner=std::unique_ptr<T>(new(std::nothrow) T);
        if(!owner)return nullptr;
        auto destroy=[](void* p){delete static_cast<T*>(p);};
        auto bytes=[](void* p)->uint64_t{return static_cast<T*>(p)->bytes();};
        if(solver_nccl_workspace_attach_kind_(kind,owner.get(),destroy,bytes)!=0)return nullptr;
        pointer=owner.release();
    }
    return static_cast<T*>(pointer);
}
}
