// ? LAPACK *larnv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarnv/dlarnv/clarnv/zlarnv
// *LARNV generates random vectors from several distributions.

slarnv :: #bake_arguments _larnv_real(RT=float32);
dlarnv :: #bake_arguments _larnv_real(RT=float64);
clarnv :: #bake_arguments _larnv_complex(RT=float32, MT=C32);
zlarnv :: #bake_arguments _larnv_complex(RT=float64, MT=C64);

_larnv_laruv :: (iseed: *$ISEED/ArrayType1D, $$n: int, u: *$U/ArrayType1D, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && ISEED.data_type == int && U.data_type == RT;
} {
    #if RT == float32 { slaruv(iseed, n, u); } else { dlaruv(iseed, n, u); }
}

_larnv_real :: (
    $$idist: int, iseed: *$ISEED/ArrayType1D, $$n: int, x: *$X/ArrayType1D, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && ISEED.data_type == int && X.data_type == RT;
} {
    one_r   :: #run one(RT);
    two_r   :: #run ncast(RT, 2.0);
    twopi_r :: #run ncast(RT, 6.28318530717958647692528676655900576839);
    LV :: 128;
    u_storage: [LV] RT;
    u_view := u_storage[0..];
    i, il, il2, iv : int = ---;

    for iv:1..n by LV / 2 {
        il = min(LV / 2, n - iv + 1);
        if idist == 3 { il2 = 2*il; } else { il2 = il; }
        _larnv_laruv(iseed, il2, *u_view, RT);
        if idist == 1 {
            for i:1..il { fset(x, iv + i - 1, u_storage[i - 1]); }
        } else if idist == 2 {
            for i:1..il { fset(x, iv + i - 1, two_r*u_storage[i - 1] - one_r); }
        } else if idist == 3 {
            for i:1..il {
                fset(x, iv + i - 1, sqrt(-two_r*log(u_storage[2*i - 2])) * cos(twopi_r*u_storage[2*i - 1]));
            }
        }
    }
}

_larnv_complex :: (
    $$idist: int, iseed: *$ISEED/ArrayType1D, $$n: int, x: *$X/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           ISEED.data_type == int && X.data_type == MT;
} {
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    two_r   :: #run ncast(RT, 2.0);
    twopi_r :: #run ncast(RT, 6.28318530717958647692528676655900576839);
    LV :: 128;
    u_storage: [LV] RT;
    u_view := u_storage[0..];
    i, il, iv : int = ---;

    for iv:1..n by LV / 2 {
        il = min(LV / 2, n - iv + 1);
        _larnv_laruv(iseed, 2*il, *u_view, RT);
        if idist == 1 {
            for i:1..il { fset(x, iv + i - 1, cpx(u_storage[2*i - 2], u_storage[2*i - 1])); }
        } else if idist == 2 {
            for i:1..il { fset(x, iv + i - 1, cpx(two_r*u_storage[2*i - 2] - one_r, two_r*u_storage[2*i - 1] - one_r)); }
        } else if idist == 3 {
            for i:1..il {
                mag := sqrt(-two_r*log(u_storage[2*i - 2]));
                ang := twopi_r*u_storage[2*i - 1];
                fset(x, iv + i - 1, mag*exp(cpx(zero_r, ang)));
            }
        } else if idist == 4 {
            for i:1..il {
                mag := sqrt(u_storage[2*i - 2]);
                ang := twopi_r*u_storage[2*i - 1];
                fset(x, iv + i - 1, mag*exp(cpx(zero_r, ang)));
            }
        } else if idist == 5 {
            for i:1..il {
                ang := twopi_r*u_storage[2*i - 1];
                fset(x, iv + i - 1, exp(cpx(zero_r, ang)));
            }
        }
    }
}
