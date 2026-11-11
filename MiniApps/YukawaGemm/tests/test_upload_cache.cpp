int main() {
    int rp[]={0,1,2}, ci[]={1,3}; double v[]={2,3};
    assert(ensure_aloc_uploaded(2,3,2,rp,ci,v)==0);
    assert(d_Aloc_col_ind[1]==2 && d_Aloc_val[1]==3);
    ci[1]=2; v[1]=7;
    g_spgemm_ab_warm=true;
    assert(ensure_aloc_uploaded(2,3,2,rp,ci,v)==0);
    assert(d_Aloc_col_ind[1]==1 && d_Aloc_val[1]==7 && !g_spgemm_ab_warm);
    assert(ensure_aloc_uploaded(2,4,2,rp,ci,v)==0);
    assert(g_desc_Aloc->columns==4);
    int br[]={0,1,1,2}, bi[]={1,2}; double bv[]={4,5};
    assert(ensure_b_uploaded(3,2,2,br,bi,bv)==0);
    br[1]=0; br[2]=1; bv[0]=8;
    assert(ensure_b_uploaded(3,2,2,br,bi,bv)==0);
    assert(d_B_row_ptr[1]==0 && d_B_row_ptr[2]==1 && d_B_val[0]==8);
    puts("Actual upload helpers: changed values, pattern, and shared dimension: PASS");
}
