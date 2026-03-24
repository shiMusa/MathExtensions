// ? LAPACK *laesy()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment claesy/zlaesy
// *LAESY computes the eigenvalues and eigenvectors of a 2-by-2
// complex symmetric matrix.

claesy :: #bake_arguments _laesy(RT=float32, MT=C32);
zlaesy :: #bake_arguments _laesy(RT=float64, MT=C64);

#scope_file
_laesy_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laesy :: ($$a: $MT, $$b: MT, $$c: MT, rt1: *$MT, rt2: *MT, evscal: *MT, cs1: *MT, sn1: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    zero_r   :: #run zero(RT);
    one_r    :: #run one(RT);
    half_r   :: #run ncast(RT, 0.5);
    thresh_r :: #run ncast(RT, 0.1);
    zero_t   :: #run zero(MT);
    one_t    :: #run _laesy_complex_from_real(one_r, RT, MT);

    babs, evnorm, tabs, z : RT = ---;
    s, t, tmp : MT = ---;

    if abs(b) == zero_r {
        rt1.* = a;
        rt2.* = c;
        if abs(rt1.*) < abs(rt2.*) {
            tmp = rt1.*;
            rt1.* = rt2.*;
            rt2.* = tmp;
            cs1.* = zero_t;
            sn1.* = one_t;
        } else {
            cs1.* = one_t;
            sn1.* = zero_t;
        }
        // The reference source leaves EVSCAL untouched in this branch.
        // Keep the diagonal case usable by returning the identity scale.
        evscal.* = one_t;
    } else {
        s = (a + c)*half_r;
        t = (a - c)*half_r;

        babs = abs(b);
        tabs = abs(t);
        z = max(babs, tabs);
        if z > zero_r {
            t = _laesy_complex_from_real(z, RT, MT)*sqrt((t / z)*(t / z) + (b / z)*(b / z));
        }

        rt1.* = s + t;
        rt2.* = s - t;
        if abs(rt1.*) < abs(rt2.*) {
            tmp = rt1.*;
            rt1.* = rt2.*;
            rt2.* = tmp;
        }

        sn1.* = (rt1.* - a) / b;
        tabs = abs(sn1.*);
        if tabs > one_r {
            t = _laesy_complex_from_real(tabs, RT, MT)*
                 sqrt((one_r / tabs)*(one_r / tabs) + (sn1.* / tabs)*(sn1.* / tabs));
        } else {
            t = sqrt(one_t + sn1.* * sn1.*);
        }
        evnorm = abs(t);
        if evnorm >= thresh_r {
            evscal.* = one_t / t;
            cs1.* = evscal.*;
            sn1.* = sn1.* * evscal.*;
        } else {
            evscal.* = zero_t;
        }
    }
}
