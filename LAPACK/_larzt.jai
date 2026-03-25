// ? LAPACK *larzt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarzt/dlarzt/clarzt/zlarzt
// *LARZT forms the triangular factor of a block reflector in rowwise storage.

slarzt :: #bake_arguments _larzt(RT=float32, MT=float32);
dlarzt :: #bake_arguments _larzt(RT=float64, MT=float64);
clarzt :: #bake_arguments _larzt(RT=float32, MT=C32);
zlarzt :: #bake_arguments _larzt(RT=float64, MT=C64);

_larzt :: (
    $$direct: u8, $$storev: u8, $$n: int, $$k: int, v: $V/ArrayType2D, $$ldv: int,
    tau: $TAU/ArrayType1D, t: *$T/ArrayType2D, $$ldt: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && TAU.data_type == MT && T.data_type == MT;
} {
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SLARZT"
    } else if MT == float64 {
        "DLARZT"
    } else if MT == C32 {
        "CLARZT"
    } else {
        "ZLARZT"
    };

    info : int = 0;
    sum : MT = ---;
    len : int = ---;

    if !lsame(direct, #char "B") {
        info = -1;
    } else if !lsame(storev, #char "R") {
        info = -2;
    }
    if info != 0 {
        xerbla(routine_name, -info);
        return;
    }

    for <i:1..k {
        if fget(tau, i) == zero_t {
            for j:i..k {
                fset(t, j, i, zero_t);
            }
        } else {
            if i < k {
                for row:i + 1..k {
                    sum = zero_t;
                    for col:1..n {
                        #if MT == C32 || MT == C64 {
                            sum += fget(v, row, col) * conjugate(fget(v, i, col));
                        } else {
                            sum += fget(v, row, col) * fget(v, i, col);
                        }
                    }
                    fset(t, row, i, -fget(tau, i) * sum);
                }

                len = k - i;
                for <row:1..len {
                    sum = zero_t;
                    for col:1..row {
                        sum += fget(t, i + row, i + col) * fget(t.*, i + col, i);
                    }
                    fset(t, i + row, i, sum);
                }
            }
            fset(t, i, i, fget(tau, i));
        }
    }
}
