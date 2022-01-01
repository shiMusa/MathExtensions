// ? LAPACK *hbmv() or *sbmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga426c198e1c8863d9ecbdbe16efa0f6cf.html#ga426c198e1c8863d9ecbdbe16efa0f6cf

// ? original comment
// CHBMV  performs the matrix-vector  operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n hermitian band matrix, with k super-diagonals.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

ssbmv :: _hbmv;
dsbmv :: _hbmv;
chbmv :: _hbmv;
zhbmv :: _hbmv;

_hbmv :: (
    $$uplo: u8,
    $$n: int,
    $$k: int,
    $$alpha: $T,
    $$a: $A/ArrayType2D,
    $$lda: int,
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
        assert(k >= 0);
        assert(lda >= k+1);
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || (alpha == zero && beta == one) then return;

    // set up the start points in x and y
    kx := ifx incx > 0 then 0 else (n-1)*(-incx);
    ky := ifx incy > 0 then 0 else (n-1)*(-incy);

    temp1, temp2 : T = ---;
    ix, iy, jx, jy, kplus1, l : int = ---;

    // start the operations. in this version the elements of the array a
    // are accessed sequentially with one pass through a.

    // first form  y := beta*y
    if beta != one {
        if incy == 1 {
            if beta == zero {
                for i:0..n-1 {
                    set(y,i, zero);
                }
            } else {
                for i:0..n-1 {
                    set(y,i,beta*get(y,i));
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

    if char_eq(uplo, #char "u") {
        // Form  y  when upper triangle of A is stored.
        kplus1 = k + 1;
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                temp1 = alpha*get(x,j);
                temp2 = zero;
                l = kplus1 - j;
                for i:max(0,j-k)..j-1 {
                    set(y,i, get(y,i) + temp1*get(a,l+i,j));
                    temp2 += conjugate(get(a,l+i,j))*get(x,i);
                }
                set(y,i, get(y,i) + temp1*real(get(a,kplus1,j)) + alpha*temp2);
            }
        } else {
            jx = kx;
            jy = ky;
            for j:0..n-1 {
                temp1 = alpha * get(x,jx);
                temp2 = zero;
                ix = kx;
                iy = ky;
                l = kplus1 - j;
                for i:max(0,j-k)..j-1 {
                    set(y,iy, get(y,iy) + temp1*get(a,l+i,j));
                    temp2 += conjugate(get(a,l+i,j))*get(x,ix);
                    ix += incx;
                    iy += incy;
                }
                set(y,iy, get(y,iy) + temp1*real(get(a,kplus1,j)) + alpha*temp2);
                jx += incx;
                jy += incy;
                if j > k {
                    kx += incx;
                    ky += incy;
                }
            }
        }
    } else {
        // Form  y  when lower triangle of A is stored.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                temp1 = alpha*get(x,j);
                temp2 = zero;
                set(y,j, get(y,j) + temp1*real(get(a,1,j)));
                l = 1 - j;
                for i:j..min(n,j+k)-1 {
                    set(y,i, get(y,i) + temp1*get(a,l+i,j));
                    temp2 += conjugate(get(a,l+i,j))*get(x,i);
                }
                set(y,j, get(y,j) + alpha*temp2);
            }
        } else {
            jx = kx;
            jy = ky;
            for j:0..n-1 {
                temp1 = alpha*get(x,jx);
                temp2 = zero;
                set(y,iy, get(y,iy) + temp1*real(get(a,1,j)));
                l = 1 - j;
                ix = jx;
                iy = jy;
                for i:j..min(n,j+k)-1 {
                    ix += incx;
                    iy += incy;
                    set(y,iy, get(y,iy) + temp1*get(a,l+i,j));
                    temp2 += conjugate(get(a,l+i,,j))*get(x,ix);
                }
                set(y,iy, get(y,iy) + alpha*temp2);
                jx += incx;
                jy += incy;
            }
        }
    }
}