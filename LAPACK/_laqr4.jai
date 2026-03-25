// ? LAPACK *laqr4()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqr4/dlaqr4/claqr4/zlaqr4
// *LAQR4 is the recursive multishift QR driver.

slaqr4 :: #bake_arguments _laqr4_real(RT=float32);
dlaqr4 :: #bake_arguments _laqr4_real(RT=float64);
claqr4 :: #bake_arguments _laqr4_complex(RT=float32, MT=C32);
zlaqr4 :: #bake_arguments _laqr4_complex(RT=float64, MT=C64);

#scope_file
_laqr4_set_work1_real :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} { fset(work, 1, x); }

_laqr4_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqr4_real :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int,
    h: *$H/ArrayType2D, $$ldh: int, wr: *$WR/ArrayType1D, wi: *$WI/ArrayType1D,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT && WR.data_type == RT && WI.data_type == RT && Z.data_type == RT && WORK.data_type == RT;
} {
    info.* = 0;
    _laqr4_set_work1_real(work, ncast(RT, max(1, n)), RT);
    if lwork == -1 || n <= 0 then return;
    _ = wantt, wantz, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz;
    // Staged port note:
    // The recursive multishift QR driver body from *laqr4.f is still pending.
}

_laqr4_complex :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int,
    h: *$H/ArrayType2D, $$ldh: int, w: *$W/ArrayType1D,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT && W.data_type == MT && Z.data_type == MT && WORK.data_type == MT;
} {
    info.* = 0;
    _laqr4_set_work1_complex(work, ncast(RT, max(1, n)), RT, MT);
    if lwork == -1 || n <= 0 then return;
    _ = wantt, wantz, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz;
    // Staged port note:
    // The recursive multishift QR driver body from *laqr4.f is still pending.
}
