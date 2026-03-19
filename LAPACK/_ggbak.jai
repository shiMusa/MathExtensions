// ? LAPACK *ggbak()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggbak/cggbak/dggbak/zggbak
// *GGBAK forms the right or left eigenvectors of a generalized matrix pair
// by backward transformation on the computed eigenvectors of the balanced pair.

sggbak :: #bake_arguments _ggbak(RT=float32, MT=float32);
dggbak :: #bake_arguments _ggbak(RT=float64, MT=float64);
cggbak :: #bake_arguments _ggbak(RT=float32, MT=C32);
zggbak :: #bake_arguments _ggbak(RT=float64, MT=C64);

#scope_file
_ggbak_matrix_row :: inline (
    v: $V/ArrayType2D,
    $$ldv: int,
    $$row: int,
    $$count: int
) -> [] V.data_type {
    view: [] V.data_type;
    view.data = v.data + (row - 1);
    view.count = if count <= 0 then 0 else 1 + (count - 1)*ldv;
    return view;
}

_ggbak_scal :: (
    $$m: int,
    s: $RT,
    v: *$V/ArrayType2D,
    $$ldv: int,
    $$row: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT;
} {
    row_view := _ggbak_matrix_row(v.*, ldv, row, m);

    #if MT == float32 {
        sscal(m, s, *row_view, ldv);
    } else if MT == float64 {
        dscal(m, s, *row_view, ldv);
    } else if MT == C32 {
        csscal(m, s, *row_view, ldv);
    } else {
        zdscal(m, s, *row_view, ldv);
    }
}

_ggbak_swap :: (
    $$m: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    $$row1: int,
    $$row2: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT;
} {
    if row1 == row2 then return;

    row1_view := _ggbak_matrix_row(v.*, ldv, row1, m);
    row2_view := _ggbak_matrix_row(v.*, ldv, row2, m);

    #if MT == float32 {
        sswap(m, *row1_view, ldv, *row2_view, ldv);
    } else if MT == float64 {
        dswap(m, *row1_view, ldv, *row2_view, ldv);
    } else if MT == C32 {
        cswap(m, *row1_view, ldv, *row2_view, ldv);
    } else {
        zswap(m, *row1_view, ldv, *row2_view, ldv);
    }
}

_ggbak :: (
    $$job: u8,
    $$side: u8,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    lscale: $LSCALE/ArrayType1D,
    rscale: $RSCALE/ArrayType1D,
    $$m: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           LSCALE.data_type == RT &&
           RSCALE.data_type == RT &&
           V.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGGBAK"
    } else if MT == float64 {
        "DGGBAK"
    } else if MT == C32 {
        "CGGBAK"
    } else {
        "ZGGBAK"
    };

    leftv, rightv : bool = ---;
    i, k : int = ---;

    rightv = lsame(side, #char "R");
    leftv  = lsame(side, #char "L");

    info.* = 0;
    if !lsame(job, #char "N") &&
       !lsame(job, #char "P") &&
       !lsame(job, #char "S") &&
       !lsame(job, #char "B") {
        info.* = -1;
    } else if !rightv && !leftv {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ilo < 1 {
        info.* = -4;
    } else if n == 0 && ihi == 0 && ilo != 1 {
        info.* = -4;
    } else if n > 0 && (ihi < ilo || ihi > max(1, n)) {
        info.* = -5;
    } else if n == 0 && ilo == 1 && ihi != 0 {
        info.* = -5;
    } else if m < 0 {
        info.* = -8;
    } else if ldv < max(1, n) {
        info.* = -10;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;
    if m == 0 then return;
    if lsame(job, #char "N") then return;

    if ilo != ihi {
        if lsame(job, #char "S") || lsame(job, #char "B") {
            if rightv {
                for i:ilo..ihi {
                    _ggbak_scal(m, fget(rscale, i), v, ldv, i, RT, MT);
                }
            }

            if leftv {
                for i:ilo..ihi {
                    _ggbak_scal(m, fget(lscale, i), v, ldv, i, RT, MT);
                }
            }
        }
    }

    if lsame(job, #char "P") || lsame(job, #char "B") {
        if rightv {
            if ilo != 1 {
                for <i:1..ilo - 1 {
                    k = cast(int) fget(rscale, i);
                    if k != i {
                        _ggbak_swap(m, v, ldv, i, k, MT);
                    }
                }
            }

            if ihi != n {
                for i:ihi + 1..n {
                    k = cast(int) fget(rscale, i);
                    if k != i {
                        _ggbak_swap(m, v, ldv, i, k, MT);
                    }
                }
            }
        }

        if leftv {
            if ilo != 1 {
                for <i:1..ilo - 1 {
                    k = cast(int) fget(lscale, i);
                    if k != i {
                        _ggbak_swap(m, v, ldv, i, k, MT);
                    }
                }
            }

            if ihi != n {
                for i:ihi + 1..n {
                    k = cast(int) fget(lscale, i);
                    if k != i {
                        _ggbak_swap(m, v, ldv, i, k, MT);
                    }
                }
            }
        }
    }
}
