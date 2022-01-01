// ? LAPACK *geru and *ger
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gab21c402efadfa2023cfbc06911506e42.html#gab21c402efadfa2023cfbc06911506e42

// ? original comment cgeru
// CGERU  performs the rank 1 operation
//
//    A := alpha*x*y**T + A,
//
// where alpha is a scalar, x is an m element vector, y is an n element
// vector and A is an m by n matrix.

// Level 2 Blas routine.
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

cgerc :: #bake_constants _geru(conjugate=true);
zgerc :: #bake_constants _geru(conjugate=true);
cgeru :: #bake_constants _geru(conjugate=false);
zgeru :: #bake_constants _geru(conjugate=false);
sger  :: #bake_constants _geru(conjugate=false);
dger  :: #bake_constants _geru(conjugate=false);

_geru :: (
    $$m: int,
    $$n: int,
    $$alpha: $T,
    x: $X/ArrayType1D,
    $$incx: int,
    y: X,
    $$incy: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $conjugate: bool
) {
    zero :: #run zero(T);
    one  :: #run one(T);

    #if CHECKS {
        #run {
            assert(T == get_number_inner_type(A.data_type));
            assert(T == get_number_inner_type(X.data_type));
        }
        assert(m >= 0);
        assert(n >= 0);
        assert(lda >= max(1,m)); // 1?
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible

    if m == 0 || n == 0 || alpha == zero then return;

    // start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.
    temp : T = ---;
    jy := ifx incy > 0 then 0 else (n-1)*(-incy);

    if incx == 1 {
        for j:0..n-1 {
            if get(y,iy) != zero {
                #if conjugate {
                    temp = alpha * conjugate(get(y,iy));
                } else {
                    temp = alpha * get(y,iy);
                }
                for i:0..m-1 {
                    set(a,i,j, get(a,i,j) + get(x,i)*temp);
                }
            }
            jy += incy;
        }
    } else {
        ix : int = ---;
        kx := ifx incx > 0 then 0 else (m-1)*(-incx);

        for j:0..n-1 {
            if get(y,iy) != zero {
                #if conjugate {
                    temp = alpha * conjugate(get(y,iy));
                } else {
                    temp = alpha * get(y,iy);
                }
                ix = kx;
                for i:0..m-1 {
                    set(a,i,j, get(a,i,j) + get(x,ix)*temp);
                    ix += incx;
                }
            }
            jy += incy;
        }
    }
}