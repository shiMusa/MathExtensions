// ? LAPACK *hetrd_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrd_2stage/zhetrd_2stage
// *HETRD_2STAGE reduces a Hermitian matrix to real symmetric tridiagonal form
// using a two-stage algorithm.

chetrd_2stage :: #bake_arguments _hetrd_2stage(RT=float32, MT=C32);
zhetrd_2stage :: #bake_arguments _hetrd_2stage(RT=float64, MT=C64);

#scope_file
_hetrd_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrd_2stage_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrd_2stage_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type;
    view.data = x.data + (offset - 1);
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrd_2stage :: (
    $$vect: u8,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    hous2: *$HOUS2/ArrayType1D,
    $$lhous2: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           TAU.data_type == MT &&
           HOUS2.data_type == MT &&
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHETRD_2STAGE" } else { "ZHETRD_2STAGE" };

    after_info, kd, ib, lwmin, lhmin, lwrk, ldab, wpos, abpos : int = ---;
    lquery, upper, wantq : bool = ---;

    info.* = 0;
    wantq = lsame(vect, #char "V");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1 || lhous2 == -1;

    kd = ilaenv2stage(1, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", vect, n, -1, -1, -1);
    ib = ilaenv2stage(2, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", vect, n, kd, -1, -1);
    if n == 0 {
        lhmin = 1;
        lwmin = 1;
    } else {
        lhmin = ilaenv2stage(3, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", vect, n, kd, ib, -1);
        lwmin = ilaenv2stage(4, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", vect, n, kd, ib, -1);
    }

    if !lsame(vect, #char "N") {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if lhous2 < lhmin && !lquery {
        info.* = -10;
    } else if lwork < lwmin && !lquery {
        info.* = -12;
    }

    if info.* == 0 {
        fset(hous2, 1, _hetrd_2stage_complex_from_real(_hetrd_2stage_roundup_lwork(ncast(RT, lhmin), RT), RT, MT));
        fset(work, 1, _hetrd_2stage_complex_from_real(_hetrd_2stage_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        fset(work, 1, _hetrd_2stage_complex_from_real(one(RT), RT, MT));
        return;
    }

    ldab = kd + 1;
    lwrk = lwork - ldab*n;
    abpos = 1;
    wpos = abpos + ldab*n;

    ab_view_raw: [] MT;
    ab_view_raw.data = work.data + (abpos - 1);
    ab_view_raw.count = ldab*n;
    ab_view := array2d_view(ab_view_raw, ldab, row_major=false);
    wrk_view := _hetrd_2stage_vector_slice(work, wpos, max(0, lwrk));

    #if MT == C32 {
        chetrd_he2hb(uplo, n, kd, a, lda, *ab_view, ldab, tau, *wrk_view, lwrk, info);
    } else {
        zhetrd_he2hb(uplo, n, kd, a, lda, *ab_view, ldab, tau, *wrk_view, lwrk, info);
    }
    if info.* != 0 {
        after_info = info.*;
        xerbla(ifx MT == C32 then "CHETRD_HE2HB" else "ZHETRD_HE2HB", -after_info);
        return;
    }

    #if MT == C32 {
        chetrd_hb2st(#char "Y", vect, uplo, n, kd, *ab_view, ldab, d, e, hous2, lhous2, *wrk_view, lwrk, info);
    } else {
        zhetrd_hb2st(#char "Y", vect, uplo, n, kd, *ab_view, ldab, d, e, hous2, lhous2, *wrk_view, lwrk, info);
    }
    if info.* != 0 {
        after_info = info.*;
        xerbla(ifx MT == C32 then "CHETRD_HB2ST" else "ZHETRD_HB2ST", -after_info);
        return;
    }

    _ = wantq;
    fset(work, 1, _hetrd_2stage_complex_from_real(_hetrd_2stage_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
}
