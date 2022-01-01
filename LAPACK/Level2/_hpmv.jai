// ? LAPACK *hpmv() and *spmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gad1dc14308bf0c1c8d9ba6ee068ac4d60.html#gad1dc14308bf0c1c8d9ba6ee068ac4d60

// ? hpmv = hermitian packed matrix * vector

// ? original comment
// CHPMV  performs the matrix-vector operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n hermitian matrix, supplied in packed form.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

sspmv :: _hpmv; // conjugate has no effect on real numbers
dspmv :: _hpmv; // conjugate has no effect on real numbers
chpmv :: _hpmv;
zhpmv :: _hpmv;

_hpmv :: (
    $$uplo: u8,
    $$n: int,
    $$alpha: $T,
    $$ap: $A/ArrayType1D,
    $$x: $X/ArrayType1D,
    $$incx: int,
    $$beta: T,
    $$y: *$Y/ArrayType1D,
    $$incy: int
) {
    zero :: #run zero(T);
    one  :: #run one(T);

    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == X.data_type);
            assert(T == Y.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(n >= 0);
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || (alpha == zero && beta == one) then return;

    // Set up the start points in x and y.
    kx := ifx incx > 0 then 0 else (n-1)*(-incx);
    ky := ifx incy > 0 then 0 else (n-1)*(-incy);

    temp1, temp2 : T = ---;
    ix, iy, jx, jy, k, kk, kx, ky : int = ---;

    // Start the operations. In this version the elements of the array AP
    // are accessed sequentially with one pass through AP.

    // First form  y := beta*y

    if beta != one {
        if incy == 1 {
            if beta == zero {
                for i:0..n-1 {
                    set(y,i, zero);
                }
            } else {
                for i:0..n-1 {
                    set(y,i, beta*get(y,i));
                }
            }
        } else {
            iy = ky;
            if beta == zero {
                for i:0..n-1 {
                    set(y,iy, zero);
                    iy += incy;
                }
            } else {
                for i:0..n-1 {
                    set(y,iy, beta*get(y,iy));
                    iy += incy;
                }
            }
        }
    }
    if alpha == zero then return;
    
    kk = 0;
    if char_eq(uplo, #char "u") {
        // Form  y  when AP contains the upper triangle.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                temp1 = alpha*get(x,j);
                temp2 = zero;
                k = kk;
                for i:0..j-1 {
                    set(y,i, get(y,i) + temp1*get(ap,k));
                    temp2 += conjugate(get(ap,k))*get(x,i);
                    k += 1;
                }
                set(y,j, get(y,j) + temp1*real_nc(get(ap,kk+j-1)) + alpha*temp2);
                kk += j;
            }
        } else {
            jx = kx;
            jy = ky;
            for j:0..n-1 {
                temp1 = alpha*get(x,jx);
                temp2 = zero;
                ix = kx;
                iy = ky;
                for k:kk..kk+j-2 {
                    set(y,iy, get(y,iy) + temp1*get(ap,k));
                    temp2 += conjugate(get(ap,k))*get(x,ix);
                    ix += incx;
                    iy += incy;
                }
                set(y,jy, get(y,jy) + temp1*real_nc(get(ap,kk+j-1)) + alpha*temp2);
                jx += incx;
                jy += incy;
                kk += j;
            }
        }
    } else {
        // Form  y  when AP containse the lower triangle.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                temp1 = alpha*get(x,j);
                temp2 = zero;
                set(y,j, get(y,j) + temp1*real_nc(get(ap,kk)));
                k = kk + 1;
                for i:j+1..n-1 {
                    set(y,i, get(y,i) + temp1*get(ap,k));
                    temp2 += conjugate(get(ap,k))*get(x,i);
                    k += 1;
                }
                set(y,j, get(y,j) + alpha * temp2);
                kk += (n-j+1); // ! +1 ?
            }
        } else {
            jx = kx;
            jy = ky;
            for j:0..n-1 {
                temp1 = alpha*get(x,jx);
                temp2 = zero;
                set(y,jy, get(y,jy) + temp1*real_nc(get(ap,kk)));
                ix = jx;
                iy = jy;
                for k:kk+1..kk+n-j {
                    ix += incx;
                    iy += incy;
                    set(y,iy, get(y,iy) + temp1*get(ap,k));
                    temp2 += conjugate(get(ap,k))*get(x,ix);
                }
                set(y,jy, get(y,jy) + alpha*temp2);
                jx += incx;
                jy += incy;
                kk += n-j+1;
            }
        }
    }
}