// ? LAPACK *larz()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarz/dlarz/clarz/zlarz
// *LARZ applies an elementary reflector produced by *TZRZF.

slarz :: #bake_arguments _larz(RT=float32, MT=float32);
dlarz :: #bake_arguments _larz(RT=float64, MT=float64);
clarz :: #bake_arguments _larz(RT=float32, MT=C32);
zlarz :: #bake_arguments _larz(RT=float64, MT=C64);

_larz :: (
    $$side: u8, $$m: int, $$n: int, $$l: int, v: $V/ArrayType1D, $$incv: int, tau: $MT,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);
    sum : MT = ---;

    if lsame(side, #char "L") {
        if tau != zero_t {
            for j:1..n {
                #if MT == C32 || MT == C64 {
                    sum = conjugate(fget(c.*, 1, j));
                    for i:1..l {
                        sum += conjugate(fget(c.*, m - l + i, j)) * fget(v, i);
                    }
                    fset(work, j, conjugate(sum));
                } else {
                    sum = fget(c.*, 1, j);
                    for i:1..l {
                        sum += fget(c.*, m - l + i, j) * fget(v, i);
                    }
                    fset(work, j, sum);
                }
            }

            for j:1..n {
                fset(c, 1, j, fget(c.*, 1, j) - tau * fget(work.*, j));
            }

            for j:1..n {
                for i:1..l {
                    fset(c, m - l + i, j, fget(c.*, m - l + i, j) - tau * fget(v, i) * fget(work.*, j));
                }
            }
        }
    } else {
        if tau != zero_t {
            for i:1..m {
                sum = fget(c.*, i, 1);
                for j:1..l {
                    sum += fget(c.*, i, n - l + j) * fget(v, j);
                }
                fset(work, i, sum);
            }

            for i:1..m {
                fset(c, i, 1, fget(c.*, i, 1) - tau * fget(work.*, i));
            }

            for i:1..m {
                for j:1..l {
                    #if MT == C32 || MT == C64 {
                        fset(c, i, n - l + j, fget(c.*, i, n - l + j) - tau * fget(work.*, i) * conjugate(fget(v, j)));
                    } else {
                        fset(c, i, n - l + j, fget(c.*, i, n - l + j) - tau * fget(work.*, i) * fget(v, j));
                    }
                }
            }
        }
    }
}
