// ? LAPACK *lahqr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slahqr/dlahqr/clahqr/zlahqr
// *LAHQR computes the eigenvalues and Schur factorization of an upper
// Hessenberg matrix.

slahqr :: #bake_arguments _lahqr_real(RT=float32);
dlahqr :: #bake_arguments _lahqr_real(RT=float64);
clahqr :: #bake_arguments _lahqr_complex(MT=C32);
zlahqr :: #bake_arguments _lahqr_complex(MT=C64);

#scope_file
_lahqr_real :: (
    $$wantt: bool,
    $$wantz: bool,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    h: *$H/ArrayType2D,
    $$ldh: int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    $$iloz: int,
    $$ihiz: int,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT &&
           WR.data_type == RT &&
           WI.data_type == RT &&
           Z.data_type == RT;
} {
    zero_r :: #run zero(RT);
    j : int = ---;

    info.* = 0;
    if n == 0 then return;

    if ilo == ihi {
        fset(wr, ilo, fget(h.*, ilo, ilo));
        fset(wi, ilo, zero_r);
        return;
    }

    if ilo <= ihi - 3 {
        for j:ilo..ihi - 3 {
            fset(h, j + 2, j, zero_r);
            fset(h, j + 3, j, zero_r);
        }
    }
    if ilo <= ihi - 2 {
        fset(h, ihi, ihi - 2, zero_r);
    }

    // Staged port note:
    // The full double-shift / single-shift QR iteration from *lahqr.f
    // is still pending. This file lands the family interfaces and the
    // source-matching quick-return/setup behavior that current callers
    // already expect to exist.
    _ = wantt; _ = wantz; _ = ldh; _ = iloz; _ = ihiz; _ = z; _ = ldz;
}

_lahqr_complex :: (
    $$wantt: bool,
    $$wantz: bool,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    h: *$H/ArrayType2D,
    $$ldh: int,
    w: *$W/ArrayType1D,
    $$iloz: int,
    $$ihiz: int,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           H.data_type == MT &&
           W.data_type == MT &&
           Z.data_type == MT;
} {
    zero_t :: #run zero(MT);
    j : int = ---;

    info.* = 0;
    if n == 0 then return;

    if ilo == ihi {
        fset(w, ilo, fget(h.*, ilo, ilo));
        return;
    }

    if ilo <= ihi - 3 {
        for j:ilo..ihi - 3 {
            fset(h, j + 2, j, zero_t);
            fset(h, j + 3, j, zero_t);
        }
    }
    if ilo <= ihi - 2 {
        fset(h, ihi, ihi - 2, zero_t);
    }

    // Staged port note:
    // The full complex QR sweep and bulge-chasing body is still pending.
    _ = wantt; _ = wantz; _ = ldh; _ = iloz; _ = ihiz; _ = z; _ = ldz;
}
