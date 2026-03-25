// ? LAPACK *larfx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarfx/dlarfx/clarfx/zlarfx
// *LARFX applies an elementary reflector using tuned small-size code paths.

slarfx :: #bake_arguments _larfx(RT=float32, MT=float32);
dlarfx :: #bake_arguments _larfx(RT=float64, MT=float64);
clarfx :: #bake_arguments _larfx(RT=float32, MT=C32);
zlarfx :: #bake_arguments _larfx(RT=float64, MT=C64);

#scope_file
_larfx_apply_left_small :: (
    $$order: int, $$n: int, v: $V/ArrayType1D, tau: $MT, c: *$C/ArrayType2D, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT && C.data_type == MT;
} {
    zero_t :: #run zero(MT);
    sum, vi, ti : MT = ---;

    for j:1..n {
        sum = zero_t;
        for i:1..order {
            vi = fget(v, i);
            #if MT == C32 || MT == C64 {
                sum += conjugate(vi)*fget(c.*, i, j);
            } else {
                sum += vi*fget(c.*, i, j);
            }
        }

        for i:1..order {
            vi = fget(v, i);
            ti = tau*vi;
            fset(c, i, j, fget(c.*, i, j) - sum*ti);
        }
    }
}

_larfx_apply_right_small :: (
    $$m: int, $$order: int, v: $V/ArrayType1D, tau: $MT, c: *$C/ArrayType2D, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT && C.data_type == MT;
} {
    zero_t :: #run zero(MT);
    sum, vi, ti : MT = ---;

    for j:1..m {
        sum = zero_t;
        for i:1..order {
            vi = fget(v, i);
            sum += vi*fget(c.*, j, i);
        }

        for i:1..order {
            vi = fget(v, i);
            #if MT == C32 || MT == C64 {
                ti = tau*conjugate(vi);
            } else {
                ti = tau*vi;
            }
            fset(c, j, i, fget(c.*, j, i) - sum*ti);
        }
    }
}

_larfx :: (
    $$side: u8, $$m: int, $$n: int, v: $V/ArrayType1D, tau: $MT, c: *$C/ArrayType2D, $$ldc: int,
    work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);

    if tau == zero_t then return;

    if lsame(side, #char "L") {
        if m >= 1 && m <= 10 {
            _larfx_apply_left_small(m, n, v, tau, c, MT);
        } else {
            #if MT == float32 {
                slarf(side, m, n, v, 1, tau, c, ldc, work);
            } else if MT == float64 {
                dlarf(side, m, n, v, 1, tau, c, ldc, work);
            } else if MT == C32 {
                clarf(side, m, n, v, 1, tau, c, ldc, work);
            } else {
                zlarf(side, m, n, v, 1, tau, c, ldc, work);
            }
        }
    } else {
        if n >= 1 && n <= 10 {
            _larfx_apply_right_small(m, n, v, tau, c, MT);
        } else {
            #if MT == float32 {
                slarf(side, m, n, v, 1, tau, c, ldc, work);
            } else if MT == float64 {
                dlarf(side, m, n, v, 1, tau, c, ldc, work);
            } else if MT == C32 {
                clarf(side, m, n, v, 1, tau, c, ldc, work);
            } else {
                zlarf(side, m, n, v, 1, tau, c, ldc, work);
            }
        }
    }
}
