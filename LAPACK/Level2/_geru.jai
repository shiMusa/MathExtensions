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

cgerc :: #bake_arguments _geru(is_conjugate=true);
zgerc :: #bake_arguments _geru(is_conjugate=true);
cgeru :: #bake_arguments _geru(is_conjugate=false);
zgeru :: #bake_arguments _geru(is_conjugate=false);
sger  :: #bake_arguments _geru(is_conjugate=false);
dger  :: #bake_arguments _geru(is_conjugate=false);

_geru :: (
    $$m: int,
    $$n: int,
    alpha: $T,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $is_conjugate: bool
) {
    zero_t :: #run zero(T);
    one_t  :: #run one(T);

    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == X.data_type);
            assert(T == Y.data_type);
        }
        assert(m >= 0);
        assert(n >= 0);
        assert(lda >= max(1,m)); // 1?
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible

    if m == 0 || n == 0 || alpha == zero_t then return;

    // start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.
    temp : T = ---;
    jy := ifx incy > 0 then 1 else 1 - (n-1)*incy;

    if incx == 1 {
        for j:1..n {
            if fget(y, jy) != zero_t {
                #if is_conjugate {
                    temp = alpha * conjugate(fget(y, jy));
                } else {
                    temp = alpha * fget(y, jy);
                }
                for i:1..m {
                    fset(a, i, j, fget(a.*, i, j) + fget(x, i)*temp);
                }
            }
            jy += incy;
        }
    } else {
        ix : int = ---;
        kx := ifx incx > 0 then 1 else 1 - (m-1)*incx;

        for j:1..n {
            if fget(y, jy) != zero_t {
                #if is_conjugate {
                    temp = alpha * conjugate(fget(y, jy));
                } else {
                    temp = alpha * fget(y, jy);
                }
                ix = kx;
                for i:1..m {
                    fset(a, i, j, fget(a.*, i, j) + fget(x, ix)*temp);
                    ix += incx;
                }
            }
            jy += incy;
        }
    }
}
