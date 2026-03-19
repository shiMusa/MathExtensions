// ? LAPACK *gghrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgghrd/cgghrd/dgghrd/zgghrd
// *GGHRD reduces a matrix pair (A, B) to generalized upper Hessenberg
// form using Givens rotations.

sgghrd :: #bake_arguments _gghrd(RT=float32, MT=float32);
dgghrd :: #bake_arguments _gghrd(RT=float64, MT=float64);
cgghrd :: #bake_arguments _gghrd(RT=float32, MT=C32);
zgghrd :: #bake_arguments _gghrd(RT=float64, MT=C64);

#scope_file
_gghrd_matrix_col_prefix :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = count;
    return view;
}

_gghrd_matrix_row_suffix :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = if count <= 0 then 0 else 1 + (count - 1)*lda;
    return view;
}

_gghrd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gghrd_conj_if_needed :: inline (x: $MT, $$do_conj: bool, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == float32 || MT == float64 {
        return x;
    } else {
        if do_conj then return conj(x);
        return x;
    }
}

_gghrd_laset :: (
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

_gghrd_lartg :: (
    f: $MT,
    g: MT,
    c: *$RT,
    s: *MT,
    r: *MT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 {
        slartg(f, g, c, s, r);
    } else if MT == float64 {
        dlartg(f, g, c, s, r);
    } else if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_gghrd_rot :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    c: $RT,
    s: $MT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else if MT == float64 {
        drot(n, x, incx, y, incy, c, s);
    } else if MT == C32 {
        crot(n, x, incx, y, incy, c, s);
    } else {
        zrot(n, x, incx, y, incy, c, s);
    }
}

_gghrd :: (
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
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           B.data_type == MT &&
           Q.data_type == MT &&
           Z.data_type == MT;
} {
    zero_t :: #run _gghrd_complex_from_real(zero(RT), RT, MT);
    one_t  :: #run _gghrd_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == float32 { "SGGHRD" } else if MT == float64 { "DGGHRD" } else if MT == C32 { "CGGHRD" } else { "ZGGHRD" };

    ilq, ilz : bool = ---;
    icompq, icompz, jcol, jrow : int = ---;
    c : RT = ---;
    s, r, temp : MT = ---;

    if lsame(compq, #char "N") {
        ilq = false;
        icompq = 1;
    } else if lsame(compq, #char "V") {
        ilq = true;
        icompq = 2;
    } else if lsame(compq, #char "I") {
        ilq = true;
        icompq = 3;
    } else {
        icompq = 0;
    }

    if lsame(compz, #char "N") {
        ilz = false;
        icompz = 1;
    } else if lsame(compz, #char "V") {
        ilz = true;
        icompz = 2;
    } else if lsame(compz, #char "I") {
        ilz = true;
        icompz = 3;
    } else {
        icompz = 0;
    }

    info.* = 0;
    if icompq <= 0 {
        info.* = -1;
    } else if icompz <= 0 {
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
    } else if (ilq && ldq < n) || ldq < 1 {
        info.* = -11;
    } else if (ilz && ldz < n) || ldz < 1 {
        info.* = -13;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if icompq == 3 {
        _gghrd_laset(#char "F", n, n, zero_t, one_t, q, ldq, MT);
    }
    if icompz == 3 {
        _gghrd_laset(#char "F", n, n, zero_t, one_t, z, ldz, MT);
    }

    if n <= 1 then return;

    for jcol:1..n - 1 {
        for jrow:jcol + 1..n {
            fset(b, jrow, jcol, zero_t);
        }
    }

    for jcol:ilo..ihi - 2 {
        for <jrow:jcol + 2..ihi {
            temp = fget(a.*, jrow - 1, jcol);
            _gghrd_lartg(temp, fget(a.*, jrow, jcol), *c, *s, *r, RT, MT);
            fset(a, jrow - 1, jcol, r);
            fset(a, jrow, jcol, zero_t);

            if n - jcol > 0 {
                row1 := _gghrd_matrix_row_suffix(a.*, lda, jrow - 1, jcol + 1, n - jcol);
                row2 := _gghrd_matrix_row_suffix(a.*, lda, jrow, jcol + 1, n - jcol);
                _gghrd_rot(n - jcol, row1, lda, row2, lda, c, s, RT, MT);
            }

            if n + 2 - jrow > 0 {
                brow1 := _gghrd_matrix_row_suffix(b.*, ldb, jrow - 1, jrow - 1, n + 2 - jrow);
                brow2 := _gghrd_matrix_row_suffix(b.*, ldb, jrow, jrow - 1, n + 2 - jrow);
                _gghrd_rot(n + 2 - jrow, brow1, ldb, brow2, ldb, c, s, RT, MT);
            }

            if ilq {
                qcol1 := _gghrd_matrix_col_prefix(q.*, ldq, 1, jrow - 1, n);
                qcol2 := _gghrd_matrix_col_prefix(q.*, ldq, 1, jrow, n);
                _gghrd_rot(n, qcol1, 1, qcol2, 1, c, _gghrd_conj_if_needed(s, true, MT), RT, MT);
            }

            temp = fget(b.*, jrow, jrow);
            _gghrd_lartg(temp, fget(b.*, jrow, jrow - 1), *c, *s, *r, RT, MT);
            fset(b, jrow, jrow, r);
            fset(b, jrow, jrow - 1, zero_t);

            acol1 := _gghrd_matrix_col_prefix(a.*, lda, 1, jrow, ihi);
            acol2 := _gghrd_matrix_col_prefix(a.*, lda, 1, jrow - 1, ihi);
            _gghrd_rot(ihi, acol1, 1, acol2, 1, c, s, RT, MT);

            if jrow - 1 > 0 {
                bcol1 := _gghrd_matrix_col_prefix(b.*, ldb, 1, jrow, jrow - 1);
                bcol2 := _gghrd_matrix_col_prefix(b.*, ldb, 1, jrow - 1, jrow - 1);
                _gghrd_rot(jrow - 1, bcol1, 1, bcol2, 1, c, s, RT, MT);
            }

            if ilz {
                zcol1 := _gghrd_matrix_col_prefix(z.*, ldz, 1, jrow, n);
                zcol2 := _gghrd_matrix_col_prefix(z.*, ldz, 1, jrow - 1, n);
                _gghrd_rot(n, zcol1, 1, zcol2, 1, c, s, RT, MT);
            }
        }
    }
}
