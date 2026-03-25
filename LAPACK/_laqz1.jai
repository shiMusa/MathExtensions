// ? LAPACK *laqz1()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqz1/dlaqz1/claqz1/zlaqz1
// *LAQZ1 provides the small starter / bulge-chasing helper used by QZ.

slaqz1 :: #bake_arguments _laqz1_real(RT=float32);
dlaqz1 :: #bake_arguments _laqz1_real(RT=float64);
claqz1 :: #bake_arguments _laqz1_complex(RT=float32, MT=C32);
zlaqz1 :: #bake_arguments _laqz1_complex(RT=float64, MT=C64);

#scope_file
_laqz1_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_laqz1_storage :: inline (data: *$T, $$count: int, $$inc: int) -> [] T {
    storage_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _laqz1_raw_slice(data, storage_count);
}

_laqz1_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_laqz1_lamch :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch("Safe minimum".data);
    } else {
        return dlamch("Safe minimum".data);
    }
}

_laqz1_lartg_complex :: (f: $MT, g: MT, c: *$RT, s: *MT, r: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_laqz1_rot_complex :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $$c: $RT, s: $MT, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        crot(n, x, incx, y, incy, c, s);
    } else {
        zrot(n, x, incx, y, incy, c, s);
    }
}

_laqz1_real :: (
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    $$sr1: $RT,
    $$sr2: RT,
    $$si: RT,
    $$beta1: RT,
    $$beta2: RT,
    v: *$V/ArrayType1D,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT &&
           V.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    w1, w2, safmin, safmax, scale1, scale2 : RT = ---;

    safmin = _laqz1_lamch(RT);
    safmax = one_r / safmin;

    w1 = beta1*fget(a, 1, 1) - sr1*fget(b, 1, 1);
    w2 = beta1*fget(a, 2, 1) - sr1*fget(b, 2, 1);
    scale1 = sqrt(abs(w1))*sqrt(abs(w2));
    if scale1 >= safmin && scale1 <= safmax {
        w1 = w1 / scale1;
        w2 = w2 / scale1;
    }

    w2 = w2 / fget(b, 2, 2);
    w1 = (w1 - fget(b, 1, 2)*w2) / fget(b, 1, 1);
    scale2 = sqrt(abs(w1))*sqrt(abs(w2));
    if scale2 >= safmin && scale2 <= safmax {
        w1 = w1 / scale2;
        w2 = w2 / scale2;
    }

    fset(v, 1, beta2*(fget(a, 1, 1)*w1 + fget(a, 1, 2)*w2) - sr2*(fget(b, 1, 1)*w1 + fget(b, 1, 2)*w2));
    fset(v, 2, beta2*(fget(a, 2, 1)*w1 + fget(a, 2, 2)*w2) - sr2*(fget(b, 2, 1)*w1 + fget(b, 2, 2)*w2));
    fset(v, 3, beta2*(fget(a, 3, 1)*w1 + fget(a, 3, 2)*w2) - sr2*(fget(b, 3, 1)*w1 + fget(b, 3, 2)*w2));

    fset(v, 1, fget(v.*, 1) + si*si*fget(b, 1, 1) / scale1 / scale2);

    if abs(fget(v.*, 1)) > safmax || abs(fget(v.*, 2)) > safmax || abs(fget(v.*, 3)) > safmax ||
       _laqz1_isnan(fget(v.*, 1), RT) || _laqz1_isnan(fget(v.*, 2), RT) || _laqz1_isnan(fget(v.*, 3), RT) {
        fset(v, 1, zero_r);
        fset(v, 2, zero_r);
        fset(v, 3, zero_r);
    }
}

_laqz1_complex :: (
    $$ilq: bool,
    $$ilz: bool,
    $$k: int,
    $$istartm: int,
    $$istopm: int,
    $$ihi: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $$nq: int,
    $$qstart: int,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    $$nz: int,
    $$zstart: int,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           Q.data_type == MT &&
           Z.data_type == MT;
} {
    zero_t :: #run zero(MT);

    c : RT = ---;
    s, temp : MT = ---;

    if k + 1 == ihi {
        _laqz1_lartg_complex(fget(b.*, ihi, ihi), fget(b.*, ihi, ihi - 1), c, s, temp, RT, MT);
        fset(b, ihi, ihi, temp);
        fset(b, ihi, ihi - 1, zero_t);
        _laqz1_rot_complex(ihi - istartm, _laqz1_storage(b.data + (ihi - 1)*ldb + (istartm - 1), ihi - istartm, 1), 1,
                                           _laqz1_storage(b.data + (ihi - 2)*ldb + (istartm - 1), ihi - istartm, 1), 1, c, s, RT, MT);
        _laqz1_rot_complex(ihi - istartm + 1, _laqz1_storage(a.data + (ihi - 1)*lda + (istartm - 1), ihi - istartm + 1, 1), 1,
                                               _laqz1_storage(a.data + (ihi - 2)*lda + (istartm - 1), ihi - istartm + 1, 1), 1, c, s, RT, MT);
        if ilz {
            _laqz1_rot_complex(nz, _laqz1_storage(z.data + (ihi - zstart)*ldz, nz, 1), 1,
                                   _laqz1_storage(z.data + (ihi - 1 - zstart)*ldz, nz, 1), 1, c, s, RT, MT);
        }
    } else {
        _laqz1_lartg_complex(fget(b.*, k + 1, k + 1), fget(b.*, k + 1, k), c, s, temp, RT, MT);
        fset(b, k + 1, k + 1, temp);
        fset(b, k + 1, k, zero_t);
        _laqz1_rot_complex(k + 2 - istartm + 1, _laqz1_storage(a.data + k*lda + (istartm - 1), k + 2 - istartm + 1, 1), 1,
                                               _laqz1_storage(a.data + (k - 1)*lda + (istartm - 1), k + 2 - istartm + 1, 1), 1, c, s, RT, MT);
        _laqz1_rot_complex(k - istartm + 1, _laqz1_storage(b.data + k*ldb + (istartm - 1), k - istartm + 1, 1), 1,
                                           _laqz1_storage(b.data + (k - 1)*ldb + (istartm - 1), k - istartm + 1, 1), 1, c, s, RT, MT);
        if ilz {
            _laqz1_rot_complex(nz, _laqz1_storage(z.data + (k + 1 - zstart)*ldz, nz, 1), 1,
                                   _laqz1_storage(z.data + (k - zstart)*ldz, nz, 1), 1, c, s, RT, MT);
        }

        _laqz1_lartg_complex(fget(a.*, k + 1, k), fget(a.*, k + 2, k), c, s, temp, RT, MT);
        fset(a, k + 1, k, temp);
        fset(a, k + 2, k, zero_t);
        _laqz1_rot_complex(istopm - k, _laqz1_storage(a.data + k*lda + k, istopm - k, lda), lda,
                                       _laqz1_storage(a.data + k*lda + (k + 1), istopm - k, lda), lda, c, s, RT, MT);
        _laqz1_rot_complex(istopm - k, _laqz1_storage(b.data + k*ldb + k, istopm - k, ldb), ldb,
                                       _laqz1_storage(b.data + k*ldb + (k + 1), istopm - k, ldb), ldb, c, s, RT, MT);
        if ilq {
            _laqz1_rot_complex(nq, _laqz1_storage(q.data + (k + 1 - qstart)*ldq, nq, 1), 1,
                                   _laqz1_storage(q.data + (k + 2 - qstart)*ldq, nq, 1), 1, c, conjugate(s), RT, MT);
        }
    }
}
