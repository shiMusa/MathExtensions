// ? LAPACK *lapll()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slapll/clapll/dlapll/zlapll
// *LAPLL measures the linear dependence of two vectors.

slapll :: #bake_arguments _lapll_real(RT=float32);
dlapll :: #bake_arguments _lapll_real(RT=float64);
clapll :: #bake_arguments _lapll_complex(RT=float32, MT=C32);
zlapll :: #bake_arguments _lapll_complex(RT=float64, MT=C64);

#scope_file
_lapll_dot :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    $RT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    #if RT == float32 {
        return sdot(n, x, incx, y, incy);
    } else {
        return ddot(n, x, incx, y, incy);
    }
}

_lapll_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        saxpy(n, alpha, x, incx, y, incy);
    } else if MT == float64 {
        daxpy(n, alpha, x, incx, y, incy);
    } else if MT == C32 {
        caxpy(n, alpha, x, incx, y, incy);
    } else {
        zaxpy(n, alpha, x, incx, y, incy);
    }
}

_lapll_larfg_real :: (
    $$n: int,
    alpha: *$RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    tau: *RT,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        slarfg(n, alpha, x, incx, tau);
    } else {
        dlarfg(n, alpha, x, incx, tau);
    }
}

_lapll_larfg_complex :: (
    $$n: int,
    alpha: *$MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    tau: *MT,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_lapll_las2 :: (
    $$f: $RT,
    $$g: RT,
    $$h: RT,
    ssmin: *RT,
    ssmax: *RT,
    $RT: Type
) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slas2(f, g, h, ssmin, ssmax);
    } else {
        dlas2(f, g, h, ssmin, ssmax);
    }
}

_lapll_conjugated_dot :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) -> MT #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        return cdotc(n, x, incx, y, incy);
    } else {
        return zdotc(n, x, incx, y, incy);
    }
}

_lapll_subvector :: inline (
    x: *$X/ArrayType1D,
    $$offset: int,
    $$count: int
) -> [] X.data_type {
    view: [] X.data_type = x.*;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_lapll_real :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    ssmin: *RT,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    a11, a12, a22, c, ssmax, tau : RT = ---;

    if n <= 1 {
        ssmin.* = zero_r;
        return;
    }

    _lapll_larfg_real(n, x.data, *x, incx, *tau, RT);
    a11 = fget(x.*, 1);
    fset(x, 1, one_r);

    c = -tau*_lapll_dot(n, x.*, incx, y.*, incy, RT);
    _lapll_axpy(n, c, x.*, incx, y, incy, RT);

    tail_count := ifx n <= 2 then 0 else 1 + (n - 3)*incy;
    _lapll_larfg_real(n - 1, y.data + incy, _lapll_subvector(y, 1 + 2*incy, tail_count), incy, *tau, RT);

    a12 = fget(y.*, 1);
    a22 = fget(y.*, 1 + incy);

    _lapll_las2(a11, a12, a22, ssmin, *ssmax, RT);
}

_lapll_complex :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    ssmin: *RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    zero_r :: #run zero(RT);
    cone   :: #run cpx(one(RT), zero(RT));

    a11, a12, a22, c, tau : MT = ---;
    ssmax : RT = ---;

    if n <= 1 {
        ssmin.* = zero_r;
        return;
    }

    _lapll_larfg_complex(n, x.data, *x, incx, *tau, MT);
    a11 = fget(x.*, 1);
    fset(x, 1, cone);

    c = -conj(tau)*_lapll_conjugated_dot(n, x.*, incx, y.*, incy, MT);
    _lapll_axpy(n, c, x.*, incx, y, incy, MT);

    tail_count := ifx n <= 2 then 0 else 1 + (n - 3)*incy;
    _lapll_larfg_complex(n - 1, y.data + incy, _lapll_subvector(y, 1 + 2*incy, tail_count), incy, *tau, MT);

    a12 = fget(y.*, 1);
    a22 = fget(y.*, 1 + incy);

    _lapll_las2(abs(a11), abs(a12), abs(a22), ssmin, *ssmax, RT);
}
