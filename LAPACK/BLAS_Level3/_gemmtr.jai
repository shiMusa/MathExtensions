// ? LAPACK *gemmtr()
// ? https://www.netlib.org/lapack/explore-html/

// ? gemmtr = generalized matrix * matrix, triangular result update

sgemmtr :: _gemmtr;
dgemmtr :: _gemmtr;
cgemmtr :: _gemmtr;
zgemmtr :: _gemmtr;

_gemmtr :: (
    $$uplo: u8,
    $$transa: u8,
    $$transb: u8,
    $$n: int,
    $$k: int,
    alpha: $T,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: T,
    c: *$C/ArrayType2D,
    $$ldc: int
) {
    zero_t :: #run zero(T);
    one_t  :: #run one(T);

    nota  := char_eq(transa, #char "n");
    notb  := char_eq(transb, #char "n");
    conja := char_eq(transa, #char "c");
    conjb := char_eq(transb, #char "c");
    nrowa := ifx nota then n else k;
    nrowb := ifx notb then k else n;
    upper := char_eq(uplo, #char "u");

    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == B.data_type);
            assert(T == C.data_type);
        }
        assert(upper || char_eq(uplo, #char "l"));
        assert(nota || conja || char_eq(transa, #char "t"));
        assert(notb || conjb || char_eq(transb, #char "t"));
        assert(n >= 0);
        assert(k >= 0);
        assert(lda >= max(1, nrowa));
        assert(ldb >= max(1, nrowb));
        assert(ldc >= max(1, n));
    }

    if n == 0 then return;

    istart, istop : int = ---;

    if alpha == zero_t {
        if beta == zero_t {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    fset(c, i, j, zero_t);
                }
            }
        } else {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    fset(c, i, j, beta*fget(c.*, i, j));
                }
            }
        }
        return;
    }

    temp : T = ---;

    if notb {
        if nota {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                if beta == zero_t {
                    for i:istart..istop {
                        fset(c, i, j, zero_t);
                    }
                } else if beta != one_t {
                    for i:istart..istop {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                }
                for l:1..k {
                    temp = alpha*fget(b, l, j);
                    for i:istart..istop {
                        fset(c, i, j, fget(c.*, i, j) + temp*fget(a, i, l));
                    }
                }
            }
        } else if conja {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    temp = zero_t;
                    for l:1..k {
                        temp += conjugate(fget(a, l, i))*fget(b, l, j);
                    }
                    if beta == zero_t {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
            }
        } else {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    temp = zero_t;
                    for l:1..k {
                        temp += fget(a, l, i)*fget(b, l, j);
                    }
                    if beta == zero_t {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
            }
        }
    } else if nota {
        if conjb {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                if beta == zero_t {
                    for i:istart..istop {
                        fset(c, i, j, zero_t);
                    }
                } else if beta != one_t {
                    for i:istart..istop {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                }
                for l:1..k {
                    temp = alpha*conjugate(fget(b, j, l));
                    for i:istart..istop {
                        fset(c, i, j, fget(c.*, i, j) + temp*fget(a, i, l));
                    }
                }
            }
        } else {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                if beta == zero_t {
                    for i:istart..istop {
                        fset(c, i, j, zero_t);
                    }
                } else if beta != one_t {
                    for i:istart..istop {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                }
                for l:1..k {
                    temp = alpha*fget(b, j, l);
                    for i:istart..istop {
                        fset(c, i, j, fget(c.*, i, j) + temp*fget(a, i, l));
                    }
                }
            }
        }
    } else if conja {
        if conjb {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    temp = zero_t;
                    for l:1..k {
                        temp += conjugate(fget(a, l, i))*conjugate(fget(b, j, l));
                    }
                    if beta == zero_t {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
            }
        } else {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    temp = zero_t;
                    for l:1..k {
                        temp += conjugate(fget(a, l, i))*fget(b, j, l);
                    }
                    if beta == zero_t {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
            }
        }
    } else {
        if conjb {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    temp = zero_t;
                    for l:1..k {
                        temp += fget(a, l, i)*conjugate(fget(b, j, l));
                    }
                    if beta == zero_t {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
            }
        } else {
            for j:1..n {
                if upper {
                    istart = 1;
                    istop  = j;
                } else {
                    istart = j;
                    istop  = n;
                }
                for i:istart..istop {
                    temp = zero_t;
                    for l:1..k {
                        temp += fget(a, l, i)*fget(b, j, l);
                    }
                    if beta == zero_t {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
            }
        }
    }
}
