// ? LAPACK *ladiv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sladiv/cladiv/dladiv/zladiv
// *LADIV performs complex division in real arithmetic, avoiding
// unnecessary overflow.

sladiv :: #bake_arguments _ladiv_real(RT=float32);
dladiv :: #bake_arguments _ladiv_real(RT=float64);
cladiv :: #bake_arguments _ladiv_complex(RT=float32, MT=C32);
zladiv :: #bake_arguments _ladiv_complex(RT=float64, MT=C64);

#scope_file
_ladiv_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach.data);
    } else {
        return dlamch(cmach.data);
    }
}

_ladiv2 :: ($$a: $RT, $$b: RT, $$c: RT, $$d: RT, $$r: RT, $$t: RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    zero_r :: #run zero(RT);
    br: RT = ---;

    if r != zero_r {
        br = b*r;
        if br != zero_r {
            return (a + br)*t;
        } else {
            return a*t + (b*t)*r;
        }
    } else {
        return (a + d*(b / c))*t;
    }
}

_ladiv1 :: ($$a: $RT, $$b: RT, $$c: RT, $$d: RT, p: *RT, q: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    one_r :: #run one(RT);
    r, t : RT = ---;
    neg_a : RT = ---;

    r = d / c;
    t = one_r / (c + d*r);
    p.* = _ladiv2(a, b, c, d, r, t, RT);
    neg_a = -a;
    q.* = _ladiv2(b, neg_a, c, d, r, t, RT);
}

_ladiv_real :: ($$a: $RT, $$b: RT, $$c: RT, $$d: RT, p: *RT, q: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    bs      :: #run ncast(RT, 2);
    half_r  :: #run ncast(RT, 0.5);
    two_r   :: #run ncast(RT, 2);
    one_r   :: #run one(RT);

    aa, bb, cc, dd, ab, cd, s, ov, un, be, eps : RT = ---;

    aa = a;
    bb = b;
    cc = c;
    dd = d;
    ab = max(abs(a), abs(b));
    cd = max(abs(c), abs(d));
    s = one_r;

    ov = _ladiv_lamch("Overflow threshold", RT);
    un = _ladiv_lamch("Safe minimum", RT);
    eps = _ladiv_lamch("Epsilon", RT);
    be = bs / (eps*eps);

    if ab >= half_r*ov {
        aa = half_r*aa;
        bb = half_r*bb;
        s = two_r*s;
    }
    if cd >= half_r*ov {
        cc = half_r*cc;
        dd = half_r*dd;
        s = half_r*s;
    }
    if ab <= un*bs / eps {
        aa = aa*be;
        bb = bb*be;
        s = s / be;
    }
    if cd <= un*bs / eps {
        cc = cc*be;
        dd = dd*be;
        s = s*be;
    }

    if abs(d) <= abs(c) {
        _ladiv1(aa, bb, cc, dd, p, q, RT);
    } else {
        _ladiv1(bb, aa, dd, cc, p, q, RT);
        q.* = -q.*;
    }

    p.* = p.* * s;
    q.* = q.* * s;
}

_ladiv_complex :: ($$x: $MT, $$y: MT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    zr, zi : RT = ---;

    _ladiv_real(real(x), imag(x), real(y), imag(y), *zr, *zi, RT);
    return cpx(zr, zi);
}
