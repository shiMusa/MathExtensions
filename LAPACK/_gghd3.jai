// ? LAPACK *gghd3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgghd3/cgghd3/dgghd3/zgghd3
// *GGHD3 reduces a matrix pair (A, B) to generalized upper Hessenberg
// form using the blocked interface and workspace query rules.

sgghd3 :: #bake_arguments _gghd3(RT=float32, MT=float32);
dgghd3 :: #bake_arguments _gghd3(RT=float64, MT=float64);
cgghd3 :: #bake_arguments _gghd3(RT=float32, MT=C32);
zgghd3 :: #bake_arguments _gghd3(RT=float64, MT=C64);

#scope_file
_gghd3_mt_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gghd3_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_gghd3_laset :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    beta: MT,
    a: *$A/ArrayType2D,
    $$lda: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == float64 {
        dlaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_gghd3 :: (
    $$compq: u8,
    $$compz: u8,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           B.data_type == MT &&
           Q.data_type == MT &&
           Z.data_type == MT &&
           WORK.data_type == MT;
} {
    zero_t :: #run _gghd3_mt_from_real(zero(RT), RT, MT);
    one_t  :: #run _gghd3_mt_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == float32 { "SGGHD3" } else if MT == float64 { "DGGHD3" } else if MT == C32 { "CGGHD3" } else { "ZGGHD3" };

    initq, initz, lquery, wantq, wantz : bool = ---;
    ierr, jcol, lwkopt, nb, nbmin, nh, nx : int = ---;

    info.* = 0;
    nb = ilaenv(1, routine_name, " "s, n, ilo, ihi, -1);
    nh = ihi - ilo + 1;
    if nh <= 1 {
        lwkopt = 1;
    } else {
        lwkopt = 6*n*nb;
    }
    fset(work, 1, _gghd3_mt_from_real(_gghd3_roundup_lwork(lwkopt, RT), RT, MT));

    initq = lsame(compq, #char "I");
    wantq = initq || lsame(compq, #char "V");
    initz = lsame(compz, #char "I");
    wantz = initz || lsame(compz, #char "V");
    lquery = lwork == -1;

    if !lsame(compq, #char "N") && !wantq {
        info.* = -1;
    } else if !lsame(compz, #char "N") && !wantz {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ilo < 1 {
        info.* = -4;
    } else if ihi > n || ihi < ilo - 1 {
        info.* = -5;
    } else if lda < max(1, n) {
        info.* = -7;
    } else if ldb < max(1, n) {
        info.* = -9;
    } else if (wantq && ldq < n) || ldq < 1 {
        info.* = -11;
    } else if (wantz && ldz < n) || ldz < 1 {
        info.* = -13;
    } else if lwork < 1 && !lquery {
        info.* = -15;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if initq {
        _gghd3_laset(#char "A", n, n, zero_t, one_t, q, ldq, MT);
    }
    if initz {
        _gghd3_laset(#char "A", n, n, zero_t, one_t, z, ldz, MT);
    }

    if n > 1 {
        for j:1..n - 1 {
            for i:j + 1..n {
                fset(b, i, j, zero_t);
            }
        }
    }

    if nh <= 1 {
        fset(work, 1, one_t);
        return;
    }

    nbmin = ilaenv(2, routine_name, " "s, n, ilo, ihi, -1);
    if nb > 1 && nb < nh {
        nx = max(nb, ilaenv(3, routine_name, " "s, n, ilo, ihi, -1));
        if nx < nh {
            if lwork < lwkopt {
                nbmin = max(2, ilaenv(2, routine_name, " "s, n, ilo, ihi, -1));
                if lwork >= 6*n*nbmin {
                    nb = lwork / (6*n);
                } else {
                    nb = 1;
                }
            }
        }
    }

    if nb < nbmin || nb >= nh {
        jcol = ilo;
    } else {
        // Staged port note: the blocked workspace-accumulation path from
        // LAPACK still needs its own literal follow-up pass, so we fall
        // through to the unblocked kernel starting at ILO.
        jcol = ilo;
    }

    if jcol < ihi {
        #if MT == float32 {
            sgghrd(compq, compz, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        } else if MT == float64 {
            dgghrd(compq, compz, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        } else if MT == C32 {
            cgghrd(compq, compz, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        } else {
            zgghrd(compq, compz, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        }
    }

    fset(work, 1, _gghd3_mt_from_real(_gghd3_roundup_lwork(lwkopt, RT), RT, MT));
}
