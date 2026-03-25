// ? LAPACK *laqr5()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqr5/dlaqr5/claqr5/zlaqr5
// *LAQR5 carries out a blocked small-bulge multishift QR sweep.

slaqr5 :: #bake_arguments _laqr5_real(RT=float32);
dlaqr5 :: #bake_arguments _laqr5_real(RT=float64);
claqr5 :: #bake_arguments _laqr5_complex(RT=float32, MT=C32);
zlaqr5 :: #bake_arguments _laqr5_complex(RT=float64, MT=C64);

#scope_file
_laqr5_real :: (
    $$wantt: bool, $$wantz: bool, $$kacc22: int, $$n: int, $$ktop: int, $$kbot: int, $$nshfts: int,
    sr: $SR/ArrayType1D, si: $SI/ArrayType1D, h: *$H/ArrayType2D, $$ldh: int,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, v: *$V/ArrayType2D, $$ldv: int,
    u: *$U/ArrayType2D, $$ldu: int, $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int,
    $$nh: int, wh: *$WH/ArrayType2D, $$ldwh: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           SR.data_type == RT && SI.data_type == RT && H.data_type == RT && Z.data_type == RT &&
           V.data_type == RT && U.data_type == RT && WV.data_type == RT && WH.data_type == RT;
} {
    _ = wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh;
    // Staged port note:
    // The blocked small-bulge multishift QR sweep from *laqr5.f is still pending.
}

_laqr5_complex :: (
    $$wantt: bool, $$wantz: bool, $$kacc22: int, $$n: int, $$ktop: int, $$kbot: int, $$nshfts: int,
    s: $S/ArrayType1D, h: *$H/ArrayType2D, $$ldh: int,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, v: *$V/ArrayType2D, $$ldv: int,
    u: *$U/ArrayType2D, $$ldu: int, $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int,
    $$nh: int, wh: *$WH/ArrayType2D, $$ldwh: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           S.data_type == MT && H.data_type == MT && Z.data_type == MT &&
           V.data_type == MT && U.data_type == MT && WV.data_type == MT && WH.data_type == MT;
} {
    _ = wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh;
    // Staged port note:
    // The blocked small-bulge multishift QR sweep from *laqr5.f is still pending.
}
