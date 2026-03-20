// ? LAPACK *gtts2()
// ? http://www.netlib.org/lapack/explore-html/
//
// ? original comment sgtts2/cgtts2/dgtts2/zgtts2
// *GTTS2 solves a system of linear equations with a tridiagonal matrix
// using the LU factorization computed by *GTTRF.
//
// In LAPACK, *GTTRS calls *GTTS2 with ITRANS decoded from TRANS.
//
sgtts2 :: #bake_arguments _gtts2(RT=float32, MT=float32);
dgtts2 :: #bake_arguments _gtts2(RT=float64, MT=float64);
cgtts2 :: #bake_arguments _gtts2(RT=float32, MT=C32);
zgtts2 :: #bake_arguments _gtts2(RT=float64, MT=C64);

#scope_file
_gtts2 :: (
    $$itrans: int,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    // Quick return if possible
    if n == 0 || nrhs == 0 {
        return;
    }

    if itrans == 0 {
        // Solve A * X = B  (No transpose)
        if nrhs <= 1 {
            j : int = 1;
            for i:1..n-1 {
                ip := fget(ipiv, i);
                // IP is either i or i+1. Index selection for B(I+1-IP+I,J):
                // IP==i   => B(i+1,J)
                // IP==i+1 => B(i,J)
                idx := if ip == i { i + 1 } else { i };
                temp : MT = fget(b.*, idx, j) - fget(dl, i)*fget(b.*, ip, j);
                fset(b, i, j, fget(b.*, ip, j));
                fset(b, i + 1, j, temp);
            }

            // Solve U*x = b.
            fset(b, n, j, fget(b.*, n, j) / fget(d, n));
            if n > 1 {
                fset(b, n - 1, j,
                    (fget(b.*, n - 1, j) - fget(du, n - 1)*fget(b.*, n, j)) / fget(d, n - 1));
            }
            if n > 2 {
                for <i:n - 2..1 {
                    fset(b, i, j,
                        (fget(b.*, i, j)
                         - fget(du, i)*fget(b.*, i + 1, j)
                         - fget(du2, i)*fget(b.*, i + 2, j)) / fget(d, i));
                }
            }
        } else {
            // NRHS > 1 path
            for j:1..nrhs {
                // Solve L*x=b.
                for i:1..n-1 {
                    if fget(ipiv, i) == i {
                        fset(b, i + 1, j,
                            fget(b.*, i + 1, j) - fget(dl, i)*fget(b.*, i, j));
                    } else {
                        // Swap B(i,J) and B(i+1,J) under pivot, then update.
                        old_i   : MT = fget(b.*, i, j);
                        old_ip1 : MT = fget(b.*, i + 1, j);
                        fset(b, i, j, old_ip1);
                        fset(b, i + 1, j, old_i - fget(dl, i)*old_ip1);
                    }
                }

                // Solve U*x=b.
                fset(b, n, j, fget(b.*, n, j) / fget(d, n));
                if n > 1 {
                    fset(b, n - 1, j,
                        (fget(b.*, n - 1, j) - fget(du, n - 1)*fget(b.*, n, j)) / fget(d, n - 1));
                }
                if n > 2 {
                    for <i:n - 2..1 {
                        fset(b, i, j,
                            (fget(b.*, i, j)
                             - fget(du, i)*fget(b.*, i + 1, j)
                             - fget(du2, i)*fget(b.*, i + 2, j)) / fget(d, i));
                    }
                }
            }
        }
    } else if itrans == 1 {
        // Solve A^T * X = B  (Transpose)
        if nrhs <= 1 {
            j : int = 1;

            // Solve U^T*x = b.
            fset(b, 1, j, fget(b.*, 1, j) / fget(d, 1));
            if n > 1 {
                fset(b, 2, j,
                    (fget(b.*, 2, j) - fget(du, 1)*fget(b.*, 1, j)) / fget(d, 2));
            }
            if n > 2 {
                for i:3..n {
                    fset(b, i, j,
                        (fget(b.*, i, j)
                         - fget(du, i - 1)*fget(b.*, i - 1, j)
                         - fget(du2, i - 2)*fget(b.*, i - 2, j)) / fget(d, i));
                }
            }

            // Solve L^T*x = b.
            for <i:n - 1..1 {
                ip := fget(ipiv, i);
                temp : MT = fget(b.*, i, j) - fget(dl, i)*fget(b.*, i + 1, j);
                fset(b, i, j, fget(b.*, ip, j));
                fset(b, ip, j, temp);
            }
        } else {
            for j:1..nrhs {
                // Solve U^T*x = b.
                fset(b, 1, j, fget(b.*, 1, j) / fget(d, 1));
                if n > 1 {
                    fset(b, 2, j,
                        (fget(b.*, 2, j) - fget(du, 1)*fget(b.*, 1, j)) / fget(d, 2));
                }
                if n > 2 {
                    for i:3..n {
                        fset(b, i, j,
                            (fget(b.*, i, j)
                             - fget(du, i - 1)*fget(b.*, i - 1, j)
                             - fget(du2, i - 2)*fget(b.*, i - 2, j)) / fget(d, i));
                    }
                }

                // Solve L^T*x = b.
                for <i:n - 1..1 {
                    ip := fget(ipiv, i);
                    temp : MT = fget(b.*, i, j) - fget(dl, i)*fget(b.*, i + 1, j);
                    fset(b, i, j, fget(b.*, ip, j));
                    fset(b, ip, j, temp);
                }
            }
        }
    } else {
        // itrans == 2: Solve A^H * X = B  (Conjugate transpose), complex only
        if nrhs <= 1 {
            j : int = 1;

            // Solve U^H * x = b.
            fset(b, 1, j, fget(b.*, 1, j) / conjugate(fget(d, 1)));
            if n > 1 {
                fset(b, 2, j,
                    (fget(b.*, 2, j) - conjugate(fget(du, 1))*fget(b.*, 1, j)) / conjugate(fget(d, 2)));
            }
            if n > 2 {
                for i:3..n {
                    fset(b, i, j,
                        (fget(b.*, i, j)
                         - conjugate(fget(du, i - 1))*fget(b.*, i - 1, j)
                         - conjugate(fget(du2, i - 2))*fget(b.*, i - 2, j)) / conjugate(fget(d, i)));
                }
            }

            // Solve L^H * x = b.
            for <i:n - 1..1 {
                ip := fget(ipiv, i);
                if ip == i {
                    fset(b, i, j,
                        fget(b.*, i, j) - conjugate(fget(dl, i))*fget(b.*, i + 1, j));
                } else {
                    temp : MT = fget(b.*, i + 1, j);
                    fset(b, i + 1, j,
                        fget(b.*, i, j) - conjugate(fget(dl, i))*temp);
                    fset(b, i, j, temp);
                }
            }
        } else {
            for j:1..nrhs {
                // Solve U^H * x = b.
                fset(b, 1, j, fget(b.*, 1, j) / conjugate(fget(d, 1)));
                if n > 1 {
                    fset(b, 2, j,
                        (fget(b.*, 2, j) - conjugate(fget(du, 1))*fget(b.*, 1, j)) / conjugate(fget(d, 2)));
                }
                if n > 2 {
                    for i:3..n {
                        fset(b, i, j,
                            (fget(b.*, i, j)
                             - conjugate(fget(du, i - 1))*fget(b.*, i - 1, j)
                             - conjugate(fget(du2, i - 2))*fget(b.*, i - 2, j)) / conjugate(fget(d, i)));
                    }
                }

                // Solve L^H * x = b.
                for <i:n - 1..1 {
                    ip := fget(ipiv, i);
                    if ip == i {
                        fset(b, i, j,
                            fget(b.*, i, j) - conjugate(fget(dl, i))*fget(b.*, i + 1, j));
                    } else {
                        temp : MT = fget(b.*, i + 1, j);
                        fset(b, i + 1, j,
                            fget(b.*, i, j) - conjugate(fget(dl, i))*temp);
                        fset(b, i, j, temp);
                    }
                }
            }
        }
    }
}

