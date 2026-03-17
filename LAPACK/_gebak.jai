// ? LAPACK *gebak()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgebak/cgebak/dgebak/zgebak
// *GEBAK forms the right or left eigenvectors of a general matrix
// by backward transformation on the computed eigenvectors of the
// balanced matrix output by *GEBAL.

sgebak :: #bake_arguments _gebak(RT=float32, MT=float32);
dgebak :: #bake_arguments _gebak(RT=float64, MT=float64);
cgebak :: #bake_arguments _gebak(RT=float32, MT=C32);
zgebak :: #bake_arguments _gebak(RT=float64, MT=C64);

#scope_file
_gebak_matrix_row :: inline (
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

_gebak_scal :: (
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
    row_view := _gebak_matrix_row(v.*, ldv, row, m);

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

_gebak_swap :: (
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

    row1_view := _gebak_matrix_row(v.*, ldv, row1, m);
    row2_view := _gebak_matrix_row(v.*, ldv, row2, m);

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

_gebak :: (
    $$job: u8,
    $$side: u8,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    scale: $SCALE/ArrayType1D,
    $$m: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           SCALE.data_type == RT &&
           V.data_type == MT;
} {
    one_r :: #run one(RT);
    routine_name :: #run if MT == float32 {
        "SGEBAK"
    } else if MT == float64 {
        "DGEBAK"
    } else if MT == C32 {
        "CGEBAK"
    } else {
        "ZGEBAK"
    };

    leftv, rightv : bool = ---;
    i, ii, k : int = ---;
    s : RT = ---;

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
    } else if ilo < 1 || ilo > max(1, n) {
        info.* = -4;
    } else if ihi < min(ilo, n) || ihi > n {
        info.* = -5;
    } else if m < 0 {
        info.* = -7;
    } else if ldv < max(1, n) {
        info.* = -9;
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
                    s = fget(scale, i);
                    _gebak_scal(m, s, v, ldv, i, RT, MT);
                }
            }

            if leftv {
                for i:ilo..ihi {
                    s = one_r / fget(scale, i);
                    _gebak_scal(m, s, v, ldv, i, RT, MT);
                }
            }
        }
    }

    if lsame(job, #char "P") || lsame(job, #char "B") {
        if rightv {
            for ii:1..n {
                i = ii;
                if i >= ilo && i <= ihi then continue;
                if i < ilo then i = ilo - ii;
                k = cast(int) fget(scale, i);
                if k == i then continue;
                _gebak_swap(m, v, ldv, i, k, MT);
            }
        }

        if leftv {
            for ii:1..n {
                i = ii;
                if i >= ilo && i <= ihi then continue;
                if i < ilo then i = ilo - ii;
                k = cast(int) fget(scale, i);
                if k == i then continue;
                _gebak_swap(m, v, ldv, i, k, MT);
            }
        }
    }
}
