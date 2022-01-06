// ? LAPACK i*amax()
// ? http://www.netlib.org/lapack/explore-html/d0/d73/group__aux__blas_ga8b093706b1e37c2ddfdd0a830fe21159.html#ga8b093706b1e37c2ddfdd0a830fe21159

// ? original comment
// ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
// jack dongarra, linpack, 3/11/78.
// modified 3/93 to return if incx .le. 0.
// modified 12/3/93, array(1) declarations changed to array(*)

isamax :: i_amax;
idamax :: i_amax;
icamax :: i_amax;
izamax :: i_amax;

i_amax :: (
    $$n: int,
    $$c: $C/ArrayType1D,
    $$incx: int,
) -> int {
    T :: #run get_number_inner_type(C.data_type);
    smax : T = ---;

    if n < 1 || incx < 0 then return -1;
    if n == 1 then return 0;
    
    res := 0;
    if incx == 1 {
        // code for increment equal to 1
        smax = _abs1(get(cx,0));
        for i:1..n-1 {
            if _abs1(get(cx,i)) > smax {
                res = i;
                smax = _abs1(get(cx,i));
            }
        }
    } else {
        // code for increment no equal to 1
        ix = 0;
        smax = _abs1(get(cx,0));
        ix += incx;
        for i:1..n-1 {
            if _abs1(get(cx,ix)) > smax {
                res = i;
                smax = _abs1(get(cx,ix));
            }
            ix += incx;
        }
    }

    return res;
}