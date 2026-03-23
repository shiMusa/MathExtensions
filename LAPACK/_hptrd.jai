// ? LAPACK *hptrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chptrd/zhptrd
// *HPTRD reduces a Hermitian packed matrix to real symmetric
// tridiagonal form.

chptrd :: #bake_arguments _hptrd(RT=float32, MT=C32);
zhptrd :: #bake_arguments _hptrd(RT=float64, MT=C64);

#scope_file
_hptrd_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hptrd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hptrd_larfg :: (
    $$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *$MT, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { clarfg(n, alpha, x, incx, tau); } else { zlarfg(n, alpha, x, incx, tau); }
}

_hptrd_hpmv :: (
    $$uplo: u8, $$n: int, alpha: $MT, ap: $AP/ArrayType1D, x: $X/ArrayType1D, $$incx: int,
    beta: MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); } else { zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); }
}

_hptrd_dotc :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type
) -> MT #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { return cdotc(n, x, incx, y, incy); } else { return zdotc(n, x, incx, y, incy); }
}

_hptrd_axpy :: (
    $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { caxpy(n, alpha, x, incx, y, incy); } else { zaxpy(n, alpha, x, incx, y, incy); }
}

_hptrd_hpr2 :: (
    $$uplo: u8, $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int,
    ap: *$AP/ArrayType1D, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT && AP.data_type == MT;
} {
    #if MT == C32 { chpr2(uplo, n, alpha, x, incx, y, incy, ap); } else { zhpr2(uplo, n, alpha, x, incx, y, incy, ap); }
}

_hptrd :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, d: *$D/ArrayType1D, e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && D.data_type == RT && E.data_type == RT && TAU.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_t  :: #run one(MT);
    half_t :: #run _hptrd_complex_from_real(cast(RT) 0.5, RT, MT);
    routine_name :: #run if MT == C32 { "CHPTRD" } else { "ZHPTRD" };

    upper : bool = ---;
    i, i1, i1i1, ii : int = ---;
    alpha, taui : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n <= 0 then return;

    if upper {
        i1 = n*(n - 1)/2 + 1;
        fset(ap, i1 + n - 1, _hptrd_complex_from_real(real(fget(ap.*, i1 + n - 1)), RT, MT));
        for <i:1..n - 1 {
            alpha = fget(ap.*, i1 + i - 1);
            x_tail := _hptrd_subvector(ap.*, i1, max(0, i - 1));
            _hptrd_larfg(i, *alpha, *x_tail, 1, *taui, MT);
            fset(e, i, real(alpha));

            if taui != zero_t {
                fset(ap, i1 + i - 1, one_t);
                v := _hptrd_subvector(ap.*, i1, i);
                tau_head := _hptrd_subvector(tau.*, 1, i);
                ap_head := _hptrd_subvector(ap.*, 1, i*(i + 1)/2);
                _hptrd_hpmv(uplo, i, taui, ap_head, v, 1, zero_t, *tau_head, 1, MT);
                alpha = -half_t*taui*_hptrd_dotc(i, tau_head, 1, v, 1, MT);
                _hptrd_axpy(i, alpha, v, 1, *tau_head, 1, MT);
                _hptrd_hpr2(uplo, i, -one_t, v, 1, tau_head, 1, *ap_head, MT);
            }

            fset(ap, i1 + i - 1, _hptrd_complex_from_real(fget(e.*, i), RT, MT));
            fset(d, i + 1, real(fget(ap.*, i1 + i)));
            fset(tau, i, taui);
            i1 -= i;
        }
        fset(d, 1, real(fget(ap.*, 1)));
    } else {
        ii = 1;
        fset(ap, 1, _hptrd_complex_from_real(real(fget(ap.*, 1)), RT, MT));
        for i:1..n - 1 {
            i1i1 = ii + n - i + 1;
            alpha = fget(ap.*, ii + 1);
            x_tail := _hptrd_subvector(ap.*, ii + 2, max(0, n - i - 1));
            _hptrd_larfg(n - i, *alpha, *x_tail, 1, *taui, MT);
            fset(e, i, real(alpha));

            if taui != zero_t {
                fset(ap, ii + 1, one_t);
                v := _hptrd_subvector(ap.*, ii + 1, n - i);
                tau_tail := _hptrd_subvector(tau.*, i, n - i);
                ap_tail := _hptrd_subvector(ap.*, i1i1, (n - i)*(n - i + 1)/2);
                _hptrd_hpmv(uplo, n - i, taui, ap_tail, v, 1, zero_t, *tau_tail, 1, MT);
                alpha = -half_t*taui*_hptrd_dotc(n - i, tau_tail, 1, v, 1, MT);
                _hptrd_axpy(n - i, alpha, v, 1, *tau_tail, 1, MT);
                _hptrd_hpr2(uplo, n - i, -one_t, v, 1, tau_tail, 1, *ap_tail, MT);
            }

            fset(ap, ii + 1, _hptrd_complex_from_real(fget(e.*, i), RT, MT));
            fset(d, i, real(fget(ap.*, ii)));
            fset(tau, i, taui);
            ii = i1i1;
        }
        fset(d, n, real(fget(ap.*, ii)));
    }
}
