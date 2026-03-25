// ? LAPACK *laqr2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqr2/dlaqr2/claqr2/zlaqr2
// *LAQR2 performs aggressive early deflation for the QR algorithm.

slaqr2 :: #bake_arguments _laqr2_real(RT=float32);
dlaqr2 :: #bake_arguments _laqr2_real(RT=float64);
claqr2 :: #bake_arguments _laqr2_complex(RT=float32, MT=C32);
zlaqr2 :: #bake_arguments _laqr2_complex(RT=float64, MT=C64);

#scope_file
_laqr2_set_work1_real :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} { fset(work, 1, x); }

_laqr2_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqr2_real :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ktop: int, $$kbot: int, $$nw: int,
    h: *$H/ArrayType2D, $$ldh: int, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, sr: *$SR/ArrayType1D, si: *$SI/ArrayType1D,
    v: *$V/ArrayType2D, $$ldv: int, $$nh: int, t: *$T/ArrayType2D, $$ldt: int,
    $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int, work: *$WORK/ArrayType1D, $$lwork: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT && Z.data_type == RT && SR.data_type == RT && SI.data_type == RT &&
           V.data_type == RT && T.data_type == RT && WV.data_type == RT && WORK.data_type == RT;
} {
    ns.* = 0;
    nd.* = 0;
    _laqr2_set_work1_real(work, ncast(RT, max(1, nw)), RT);
    _ = wantt, wantz, n, ktop, kbot, h, ldh, iloz, ihiz, z, ldz, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, lwork;
    // Staged port note:
    // The aggressive early deflation body from *laqr2.f is still pending.
}

_laqr2_complex :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ktop: int, $$kbot: int, $$nw: int,
    h: *$H/ArrayType2D, $$ldh: int, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, sh: *$SH/ArrayType1D,
    v: *$V/ArrayType2D, $$ldv: int, $$nh: int, t: *$T/ArrayType2D, $$ldt: int,
    $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int, work: *$WORK/ArrayType1D, $$lwork: int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT && Z.data_type == MT && SH.data_type == MT &&
           V.data_type == MT && T.data_type == MT && WV.data_type == MT && WORK.data_type == MT;
} {
    ns.* = 0;
    nd.* = 0;
    _laqr2_set_work1_complex(work, ncast(RT, max(1, nw)), RT, MT);
    _ = wantt, wantz, n, ktop, kbot, h, ldh, iloz, ihiz, z, ldz, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, lwork;
    // Staged port note:
    // The aggressive early deflation body from *laqr2.f is still pending.
}
