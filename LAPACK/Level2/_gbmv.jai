// ? LAPACK *gbmv
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga50fe7a70e7ae8f9ede50acd6747510d7.html#ga50fe7a70e7ae8f9ede50acd6747510d7

// "gbmv" = general band matrix * vector

// ? original comment cgbmv()
// CGBMV  performs one of the matrix-vector operations
// 
// y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
// 
// y := alpha*A**H*x + beta*y,
// 
// where alpha and beta are scalars, x and y are vectors and A is an
// m by n band matrix, with kl sub-diagonals and ku super-diagonals.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//  Jack Dongarra, Argonne National Lab.
//  Jeremy Du Croz, Nag Central Office.
//  Sven Hammarling, Nag Central Office.
//  Richard Hanson, Sandia National Labs.


sgbmv :: _gbmv;
dgbmv :: _gbmv;
cgbmv :: _gbmv;
zgbmv :: _gbmv;

_gbmv :: (
    $$trans: u8, 
    $$m: int,             // rows of a
    $$n: int,             // cols of a
    $$kl: int,            // sub-diagonals of a, 0 <= kl
    $$ku: int,            // super-diagonals of a, 0 <= ku
    alpha: $T,          // scalar
    a: $A/ArrayType2D,  // band matrix of dimension lda x n
    $$lda: int,           // first dimensions of a, kl+ku+1 <= lda
    x: $X/ArrayType1D,  // vector
    $$incx: int,          // increment of x, incx != 0
    beta: T,            // scalar
    y: *$Y/ArrayType1D, // vector, will be overwritten 
    $$incy: int           // increment of y, incy != 0
) {
    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == X.data_type);
            assert(T == Y.data_type);
        }
        assert(m >= 0);
        assert(n >= 0);
        assert(kl >= 0);
        assert(ku >= 0);
        assert(lda >= kl+ku);
        if char_eq(trans, #char "n") {
            assert((n-1)*abs(incx) >= dim(x));
            assert((m-1)*abs(incy) >= dim(y));
        } else {
            assert((m-1)*abs(incx) >= dim(x));
            assert((n-1)*abs(incy) >= dim(y));
        }
        assert(incx != 0);
        assert(incy != 0);
    }

    temp: T = ---;
    ix, iy, jx, jy, k, kup1, kx, ky, lenx, leny: int = ---;
    
    zero :: #run zero(T);
    one  :: #run one(T);

    // quick return if possible
    if m == 0 || n == 0 || (alpha == zero && beta == one) then return;

    #if is_consant(trans) {
        noconj :: #run char_eq(trans, #char "t");
    } else {
        noconj := char_eq(trans, #char "t");
    }

    // set lenx and leny, the lengts of the vectors x and y and set up the start points in x and y
    if char_eq(trans, #char "n") {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    kx = ifx incx > 0 then 0 else (lenx-1)*(-incx);
    ky = ifx incy > 0 then 0 else (leny-1)*(-incy);

    // start the operations. in this version the elements of a are 
    // accessed sequentially with one pass through the band part of a.

    // first form y := beta*y.

    if beta != one {
        if incy == 1 {
            if beta == zero {
                for i:0..leny-1 {
                    set(y, i, zero);
                }
            } else {
                for i:0..leny-1 {
                    set(y,i,get(y,i)*beta);
                }
            }
        } else {
            iy = ky;
            if beta == zero {
                for i:0..leny-1 {
                    set(y,iy,zero);
                    iy += incy;
                }
            } else {
                for i:0..leny-1 {
                    set(y,iy,beta*get(y,iy));
                    iy += incy;
                }
            }
        }
    }

    if alpha == zero then return;
    kyp1 = ku + 1;
    if char_eq(trans, #char "n") {
        // form y := alpha*A*x + y

        jx = kx;
        if incy == 1 {
            for j:0..n-1 {
                tmp = alpha*get(x,jx);
                k = kup1 - j;
                for i:max(0,j-ku)..min(m,j+kl)-1 {
                    set(y,i,get(y,i) + temp * get(a,k+i,j));
                }
                jx += incx;
            }
        } else {
            for j:0..n-1 {
                temp = alpha*get(x,jx);
                iy = ky;
                k = kup1 - j;
                for i:max(1,j-ku)..min(m,j+kl)-1 {
                    set(y,iy, get(y,iy) + temp*get(a,k+i,j));
                    iy += incy;
                }
                jx += incx;
                if j > ku then ky += incy;
            }
        }
    } else {
        // form y := alpha*A**T*x + y  or  y := alpha*A**H*x + y

        jy = ky;
        if incx == 1 {
            for j:0..n-1 {
                temp = zero;
                k = kup1 - j;
                if noconj {
                    for i:max(0,j-ku)..min(m,j+kl)-1 {
                        temp += get(a,k+i,j)*get(x,i);
                    }
                } else {
                    for i:max(0,j-ku)..min(m,j+kl)-1 {
                        temp += conjugate(get(a,k+i,j))*get(x,i);
                    }
                }
                set(y,jy, get(y,jy) + alpha*temp);
                jy += incy;
            }
        } else {
            for j:0..n-1 {
                temp = zero;
                ix = kx;
                k = kup1 - j;
                if noconj {
                    for i:max(0,j-ku)..min(m,j+kl)-1 {
                        temp += get(a,k+i,j)*x(ix);
                        ix += incx;
                    }
                } else {
                    for i:max(0,j-ku)..min(m,j+kl)-1 {
                        temp += conjugate(get(a,k+i,j))*get(x,ix);
                        ix += incx;
                    }
                }
                set(y,jy, get(y,jy) + alpha*temp);
                jy += incy;
                if j > ku then kx += incx;
            }
        }
    }
}