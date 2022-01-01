// ? LAPACK *hemv() and *symv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gac23d51cad7f51484371119bb6e3fb1f3.html#gac23d51cad7f51484371119bb6e3fb1f3

// ? original comment
// CHEMV  performs the matrix-vector  operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n hermitian matrix.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.


// symv = symmetric matrix * vector
// hemv = hermitian matrix * vector
ssymv :: _hemv; // conjugate() has no effect on real numbers
dsymv :: _hemv; // conjugate() has no effect on real numbers
chemv :: _hemv;
zhemv :: _hemv;


_hemv :: (
    $$uplo: u8,
    $$n: int,
    $$alpha: $T,
    $$a: $A/ArrayType2D,
    $$lda: int,
    $$x: $X/ArrayType1D,
    $$incx: int,
    $$beta: T,
    $$y: *X,
    $$incy: int
) {
    zero :: #run zero(T);
    one  :: #run one(T);

    #if CHECKS {
        #run {
            assert(T == get_number_inner_type(A.data_type));
            assert(T == get_number_inner_type(X.data_type));
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(n >= 0);
        assert(lda >= max(1,n)); // 1 offset?
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || (alpha == zero && beta == one) then return;

    // set up the start points in x and y.
    kx := ifx incx > 0 then 0 else (n-1)*(-incx);
    ky := ifx incy > 0 then 0 else (n-1)*(-incy);

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through the triangular part

    // first form  y := beta*y
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

    if char_eq(uplo, #char "u") {
        // form  y  when A is stored in upper triangle.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                temp1 = alpha*get(x,i);
                temp2 = zero;
                for i:0..j-1 {
                    set(y,i, get(y,i) + temp1*get(a,i,j));
                    temp2 += conjugate(get(a,i,j))*get(x,i);
                }
                set(y,j, get(y,j) + temp1*real(get(a,j,j)) + alpha*temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            for j:0..n-1 {
                temp1 = alpha*get(x,jx);
                temp2 = zero;
                ix = kx;
                iy = ky;
                for i:0..j-1 {
                    set(y,iy, get(y,iy) + temp1*get(a,i,j));
                    temp2 += conjugate(get(a,i,j))*get(x,ix);
                    ix += incx;
                    iy += incy;
                }
                set(y,jy, get(y,jy) + temp1*real(get(a,j,j)) + alpha*temp2);
                jx += incx;
                jy += incy;
            }
        }
    } else {
        // form  y  when A is stored in lower triangle.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                temp1 = alpha*get(x,j);
                temp2 = zero;
                set(y,j, get(y,j) + temp1 * real(get(a,j,j)));
                for i:j..n-1 {
                    set(y,i, get(y,i) + temp1*get(a,i,j));
                    temp2 += conjugate(get(a,i,j))*get(x,i);
                }
                set(y,j, get(y,j) + alpha*temp2);
            }
        } else {
            jx = kx;
            jy = ky;
            for j:0..n-1 {
                temp1 = alpha*get(x,jx);
                temp2 = zero;
                set(y,jy, get(y,jy) + temp1*real(get(a,j,j)));
                ix = jx;
                iy = jy;
                for j:0..n-1 {
                    ix += incx;
                    iy += incy;
                    set(y,iy, get(y,iy) + temp1*get(a,i,j));
                    temp2 += conjugate(get(a,i,j))*get(x,ix);
                }
                set(y,iy, get(y,iy) + alpha*temp2);
                jx += incx;
                jy += incy;
            }
        }
    }
}