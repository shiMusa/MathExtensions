
// ? LAPACK dtrsv
// ? http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_gad2a01dd62718b28e35b752dbad8474ab.html#gad2a01dd62718b28e35b752dbad8474ab

// ? LAPACK ctrsv
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga38ae7f5a40e16fd315011faddf4412cc.html#ga38ae7f5a40e16fd315011faddf4412cc


// ? Origianl comment ctrsv
// CTRSV  solves _one of the systems of equations
//
// A*x = b,   or   A**T*x = b, or A**H*x = b,
//
// where b and x are n element vectors and A is an n by n unit, or
// non-unit, upper or lower triangular matrix.
//
// No test for singularity or near-singularity is included in this
// routine. Such tests must be performed before calling this routine.


// UPLO is CHARACTER*1
// On entry, UPLO specifies whether the matrix is an upper or
// lower triangular matrix as follows:
//
//     UPLO = 'U' or 'u'   A is an upper triangular matrix.
//
//     UPLO = 'L' or 'l'   A is a lower triangular matrix.

// TRANS is CHARACTER*1
// On entry, TRANS specifies the equations to be solved as
// follows:
//
//     TRANS = 'N' or 'n'   A*x = b.
//
//     TRANS = 'T' or 't'   A**T*x = b.
//
//     TRANS = 'C' or 'c'   A**H*x = b.

	
// DIAG is CHARACTER*1
// On entry, DIAG specifies whether or not A is unit
// triangular as follows:
//
//     DIAG = 'U' or 'u'   A is assumed to be unit triangular.
//
//     DIAG = 'N' or 'n'   A is not assumed to be unit
//                         triangular.

	
// N is INTEGER
// On entry, N specifies the order of the matrix A.
// N must be at least _zero.

// A is COMPLEX array, dimension ( LDA, N )
//     Before entry with  UPLO = 'U' or 'u', the leading n by n
//     upper triangular part of the array A must contain the upper
//     triangular matrix and the strictly lower triangular part of
//     A is not referenced.
//     Before entry with UPLO = 'L' or 'l', the leading n by n
//     lower triangular part of the array A must contain the lower
//     triangular matrix and the strictly upper triangular part of
//     A is not referenced.
//     Note that when  DIAG = 'U' or 'u', the diagonal elements of
//     A are not referenced either, but are assumed to be unity.

// LDA is INTEGER
//         On entry, LDA specifies the first dimension of A as declared
//         in the calling (sub) program. LDA must be at least
//         max( 1, n ).

// X is COMPLEX array, dimension at least
//         ( 1 + ( n - 1 )*abs( INCX ) ).
//         Before entry, the incremented array X must contain the n
//         element right-hand side vector b. On exit, X is overwritten
//         with the solution vector x.

// INCX is INTEGER
// On entry, INCX specifies the increment for the elements of
// X. INCX must not be _zero.
//

// -- Written on 22-October-1986.
// Jack Dongarra, Argonne National Lab.
// Jeremy Du Croz, Nag Central Office.
// Sven Hammarling, Nag Central Office.
// Richard Hanson, Sandia National Labs.

// * character is u8 in Jai, use it via e.g. #char "c"  - F. Fehse, 2021-11-17

strsv :: _trsv;
dtrsv :: _trsv;
ctrsv :: _trsv;
ztrsv :: _trsv;

_trsv :: (
    $$uplo: u8, 
    $$trans: u8, 
    $$diag: u8, 
    $$n: int, 
    a: $A/ArrayType2D, 
    $$lda: int, 
    x: *$X/ArrayType1D, 
    $$incx: int
) {
    T :: A.data_type;
    _zero :: #run zero(T);
    _one  :: #run one(T);
    
    #if CHECKS {
        #run {
            assert(T == X.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(char_eq(trans, #char "n", #char "t", #char "c"));
        assert(char_eq(diag, #char "u", #char "n"));
        assert(n >= 0);
        assert(lda >= max(1,n));
        assert(incx != 0);
    }

    // Quick return if possible
    if n == 0 then return;

    noconj := char_eq(trans, #char"t");
    nounit := char_eq(diag , #char"n");

    // Set up the start point in X if the increment is not unity. This 
    // will be (n-1)*incx too small for descending loops.
    kx := ifx incx > 0 then 1 else 1 - (n-1)*incx;

    tmp : T = ---;
    ix, jx : int = ---;

    // Start the operations. In this version the elements of A are 
    // accessed sequentially with _one pass through A.
    if char_eq(trans, #char"n") {
        // Form x := inv(A)*x.
        if char_eq(uplo, #char"u") {
            if incx == 1 {
                for <j: 1..n {
                    if fget(x, j)!= _zero {
                        if nounit then fset(x, j, fget(x.*, j) / fget(a, j, j));
                        tmp = fget(x, j);
                        if j > 1 for <i: 1..j-1 {
                            fset(x, i, fget(x.*, i) - tmp * fget(a, i, j));
                        }
                    }
                }
            } else {
                jx = kx + (n-1)*incx;
                for <j: 1..n {
                    if fget(x, jx)!= _zero {
                        if nounit then fset(x, jx, fget(x.*, jx) / fget(a, j, j));
                        tmp = fget(x, jx);
                        ix = jx;
                        if j > 1 for <i: 1..j-1 {
                            ix -= incx;
                            fset(x, ix, fget(x.*, ix) - tmp * fget(a, i, j));
                        }
                    }
                    jx -= incx;
                }
            }
        } else {
            if incx == 1 {
                for j: 1..n {
                    if fget(x, j)!= _zero {
                        if nounit then fset(x, j, fget(x.*, j) / fget(a, j, j));
                        tmp = fget(x, j);
                        for i: j+1..n {
                            fset(x, i, fget(x.*, i) - tmp * fget(a, i, j));
                        }
                    }
                }
            } else {
                jx = kx;
                for j: 1..n {
                    if fget(x, jx)!= _zero {
                        if nounit then fset(x, jx, fget(x.*, jx) / fget(a, j, j));
                        tmp = fget(x, jx);
                        ix = jx;
                        for i: j+1..n {
                            ix += incx;
                            fset(x, ix, fget(x.*, ix) - tmp * fget(a, i, j));
                        }
                    }
                    jx += incx;
                }
            }
        }
    } else {
        // Form x := inv( A**T )*x  or  x := inv( A**H )*x.
        if char_eq(uplo, #char"u") {
            if incx == 1 {
                for j: 1..n {
                    tmp = fget(x, j);
                    if noconj {
                        for i: 1..j-1 {
                            tmp -= fget(a, i, j)*fget(x, i);
                        }
                        if nounit then tmp /= fget(a, j, j);
                    } else {
                        for i: 1..j-1 {
                            tmp -= conjugate(fget(a, i, j))*fget(x, i);
                        }
                        if nounit then tmp /= conjugate(fget(a, j, j));
                    }
                    fset(x, j, tmp);
                }
            } else {
                jx = kx;
                for j: 1..n {
                    ix = kx;
                    tmp = fget(x, jx);
                    if noconj {
                        for i: 1..j-1 {
                            tmp -= fget(a, i, j)*fget(x, ix);
                            ix += incx;
                        }
                        if nounit then tmp /= fget(a, j, j);
                    } else {
                        for i: 1..j-1 {
                            tmp -= conjugate(fget(a, i, j))*fget(x, ix);
                            ix += incx;
                        }
                        if nounit then tmp /= conjugate(fget(a, j, j));
                    }
                    fset(x, jx, tmp);
                    jx += incx;
                }
            }
        } else {
            if incx == 1 {
                for <j: 1..n {
                    tmp = fget(x, j);
                    if noconj {
                        for <i: j+1..n {
                            tmp -= fget(a, i, j)*fget(x, i);
                        }
                        if nounit then tmp /= fget(a, j, j);
                    } else {
                        for <i: j+1..n {
                            tmp -= conjugate(fget(a, i, j))*fget(x, i);
                        }
                        if nounit then tmp /= conjugate(fget(a, j, j));
                    }
                    fset(x, j, tmp);
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j: 1..n {
                    ix = kx;
                    tmp = fget(x, jx);
                    if noconj {
                        for <i: j+1..n {
                            tmp -= fget(a, i, j)*fget(x, ix);
                            ix -= incx;
                        }
                        if nounit then tmp /= fget(a, j, j);
                    } else {
                        for <i: j+1..n {
                            tmp -= conjugate(fget(a, i, j))*fget(x, ix);
                            ix -= incx;
                        }
                        if nounit then tmp /= conjugate(fget(a, j, j));
                    }
                    fset(x, jx, tmp);
                    jx -= incx;
                }
            }
        }
    }
}
