// ? LAPACK *laev2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaev2/claev2/dlaev2/zlaev2
// *LAEV2 computes the eigenvalues and eigenvectors of a 2-by-2
// symmetric/Hermitian matrix.

slaev2 :: #bake_arguments _laev2_real(RT=float32);
dlaev2 :: #bake_arguments _laev2_real(RT=float64);
claev2 :: #bake_arguments _laev2_complex(RT=float32, MT=C32);
zlaev2 :: #bake_arguments _laev2_complex(RT=float64, MT=C64);

#scope_file
_laev2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laev2_real :: ($$a: $RT, $$b: RT, $$c: RT, rt1: *RT, rt2: *RT, cs1: *RT, sn1: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    half_r :: #run ncast(RT, 0.5);

    sgn1, sgn2 : int = ---;
    ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm, tb, tn : RT = ---;

    sm = a + c;
    df = a - c;
    adf = abs(df);
    tb = b + b;
    ab = abs(tb);
    if abs(a) > abs(c) {
        acmx = a;
        acmn = c;
    } else {
        acmx = c;
        acmn = a;
    }

    if adf > ab {
        rt = adf*sqrt(one_r + (ab / adf)*(ab / adf));
    } else if adf < ab {
        rt = ab*sqrt(one_r + (adf / ab)*(adf / ab));
    } else {
        rt = ab*sqrt(two_r);
    }

    if sm < zero_r {
        rt1.* = half_r*(sm - rt);
        sgn1 = -1;
        rt2.* = (acmx / rt1.*)*acmn - (b / rt1.*)*b;
    } else if sm > zero_r {
        rt1.* = half_r*(sm + rt);
        sgn1 = 1;
        rt2.* = (acmx / rt1.*)*acmn - (b / rt1.*)*b;
    } else {
        rt1.* = half_r*rt;
        rt2.* = -half_r*rt;
        sgn1 = 1;
    }

    if df >= zero_r {
        cs = df + rt;
        sgn2 = 1;
    } else {
        cs = df - rt;
        sgn2 = -1;
    }
    acs = abs(cs);
    if acs > ab {
        ct = -tb / cs;
        sn1.* = one_r / sqrt(one_r + ct*ct);
        cs1.* = ct*sn1.*;
    } else {
        if ab == zero_r {
            cs1.* = one_r;
            sn1.* = zero_r;
        } else {
            tn = -cs / tb;
            cs1.* = one_r / sqrt(one_r + tn*tn);
            sn1.* = tn*cs1.*;
        }
    }

    if sgn1 == sgn2 {
        tn = cs1.*;
        cs1.* = -sn1.*;
        sn1.* = tn;
    }
}

_laev2_complex :: ($$a: $MT, $$b: MT, $$c: MT, rt1: *$RT, rt2: *RT, cs1: *RT, sn1: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    t : RT = ---;
    w : MT = ---;

    if abs(b) == zero_r {
        w = _laev2_complex_from_real(one_r, RT, MT);
    } else {
        w = conjugate(b) / abs(b);
    }
    _laev2_real(real(a), abs(b), real(c), rt1, rt2, cs1, *t, RT);
    sn1.* = w*_laev2_complex_from_real(t, RT, MT);
}
