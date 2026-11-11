// Included after the actual operation implementations, with serial launches.
static void fill(int slot,const std::vector<int>& cols,const std::vector<double>& vals) {
    int m=1,n=cols.size(),err,*rp,*ci; double *v;
    gpu_ops_alloc_(&slot,&m,&n,&rp,&ci,&v,&err);
    assert(err==0); rp[0]=0; rp[1]=n;
    for(int j=0;j<n;++j) { ci[j]=cols[j]; v[j]=vals[j]; }
}
int main() {
    int init_error, device;
    gpu_ops_init_(&init_error, &device);
    assert(init_error==0);
    int a=1,b=2,out;
    double alpha=1,threshold=0;
    fill(a,{1},{1}); fill(b,{2,3,4},{2,3,4});
    gpu_ops_increment_threshold_(&a,&b,&alpha,&threshold,&out);
    assert(out==4);
    fill(b,{5,6},{5,6});
    gpu_ops_increment_threshold_(&a,&b,&alpha,&threshold,&out);
    assert(out==6); // Former scratch-capacity overrun; checked by ASan.
    int *rp,*ci,m,n; double *v;
    gpu_ops_get_(&a,&rp,&ci,&v,&m,&n);
    for(int j=0;j<6;++j) { assert(ci[j]==j+1); assert(v[j]==j+1); }
    int columns=8,diagonal=2;double lower[8],upper[8];
    assert(gpu_ops_bounds_(&a,&columns,&diagonal,lower,upper)==0);
    for(int j=0;j<8;++j) {
        double value=j<6?double(j+1):0;
        assert(upper[j]==value && lower[j]==(j==2?value:-value));
    }
    fill(a,{1,2,3},{1,2,3}); fill(b,{1,2,4},{-1,-1.99,.1});
    threshold=.1;
    gpu_ops_increment_threshold_(&a,&b,&alpha,&threshold,&out);
    assert(out==1);
    gpu_ops_get_(&a,&rp,&ci,&v,&m,&n);
    assert(ci[0]==3 && v[0]==3); // exact cancellation, below/equal threshold
    int zero=0,err;
    gpu_ops_alloc_(&a,&m,&zero,&rp,&ci,&v,&err);
    assert(err==0 && rp[0]==0 && rp[1]==0);
    gpu_ops_copy_(&b,&a);
    gpu_ops_get_(&b,&rp,&ci,&v,&m,&n);
    assert(n==0 && rp[1]==0);
    // Repeated changing dimensions exercise row-pointer ownership too.
    for(int rows: {3,1,5,2,8,1}) {
        for(int slot: {a,b}) {
            gpu_ops_alloc_(&slot,&rows,&zero,&rp,&ci,&v,&err);
            assert(err==0);
        }
        gpu_ops_increment_threshold_(&a,&b,&alpha,&threshold,&out);
        assert(out==0);
    }
    // Library owners must coexist beyond four matrices and release independently.
    std::vector<int> owners(32);
    for(int &owner:owners) gpu_ops_create_(&owner);
    for(int j=0;j<32;++j) { assert(owners[j]>4); fill(owners[j],{1},{double(j+1)}); }
    for(int j=0;j<32;++j) {
        gpu_ops_get_(&owners[j],&rp,&ci,&v,&m,&n);
        assert(n==1 && v[0]==double(j+1));
    }
    int reused=owners[7]; gpu_ops_destroy_(&owners[7]); assert(owners[7]==0);
    gpu_ops_create_(&owners[7]); assert(owners[7]==reused);
    fill(owners[7],{1},{77});
    for(int &owner:owners) { gpu_ops_destroy_(&owner); assert(owner==0); gpu_ops_destroy_(&owner); }
    puts("Operation capacity, threshold, and empty-slot checks: PASS");
}
