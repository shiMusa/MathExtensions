
// ? LAPACK *rotg
// ? http://www.netlib.org/lapack/explore-html/df/d28/group__single__blas__level1_gaafa91c51f75df6c3f2182032a221c2db.html#gaafa91c51f75df6c3f2182032a221c2db
// ? http://www.netlib.org/lapack/explore-html/df/d28/group__single__blas__level1_ga97ce56a6808661b36ab62269393ac10e.html#ga97ce56a6808661b36ab62269393ac10e


// ? original comment zrotg
// The computation uses the formulas
// |x| = sqrt( Re(x)**2 + Im(x)**2 )
// sgn(x) = x / |x|  if x /= 0
//         = 1        if x  = 0
// c = |a| / sqrt(|a|**2 + |b|**2)
// s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
// When a and b are real and r /= 0, the formulas simplify to
// r = sgn(a)*sqrt(|a|**2 + |b|**2)
// c = a / r
// s = b / r
// the same as in ZROTG when |a| > |b|.  When |b| >= |a|, the
// sign of c and s will be different from those computed by ZROTG
// if the signs of a and b are not the same.

// Anderson E. (2017)
// Algorithm 978: Safe Scaling in the Level 1 BLAS
// ACM Trans Math Softw 44:1--28
// https://doi.org/10.1145/3061665

zrotg :: (a: *$C, b: C, c: *$R, s: *C) #modify {
    // ? make sure the inner types are of same resolution
    if !is_complex_number(C) || is_non_real_number(R) || get_complex_inner_type(C) != R {
        C = null; return;
    }
} {
    czero := zero(C);
    zero  := zero(R);
    one   := one(R);
    #if R == float64 {
        safmin :: SAFMIN_64;
        safmax :: SAFMAX_64;
        rtmin  :: RTMIN_64;
        rtmax  :: RTMAX_64;
    } else {
        safmin :: SAFMIN_32;
        safmax :: SAFMAX_32;
        rtmin  :: RTMIN_32;
        rtmax  :: RTMAX_32;
    }
    // Local Scalars
    d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w : R = ---;
    f, fs, g, gs, r, t : C = ---;

    f = <<a;
    g = b;
    if g == czero {
        <<c = one;
        <<s = czero;
        r = f;
    } else if f == czero {
        <<c = zero;
        g1 = max(abs(real(g), abs(imag(g))));
        if g1 > rtmin && g1 < rtmax {
            // use unscaled algorithm
            g2 = abs_sq(g);
            d = sqrt(g2);
            <<s = conjugate(g)/d;
            r = d;
        } else {
            // use scaled algorithms
            u = clamp(ncast(R,safmin), g1, ncast(R,safmax));
            uu = one / u;
            gs = g*uu;
            g2 = abs_sq(gs);
            d = sqrt(g2);
            <<s = conjugate(gs)/d;
            r = d*u;
        }
    } else {
        f1 = max( abs(real(f)), abs(imag(f)) );
        g1 = max( abs(real(g)), abs(imag(g)) );
        if f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmin {
            // use unscaled algorithm
            f2 = abs_sq(f);
            g2 = abs_sq(g);
            h2 = f2 + g2;
            if f2 > rtmin && h2 < rtmax {
                d = sqrt(f2*h2);
            } else {
                d = sqrt(f2)*sqrt(h2);
            }
            p = one/d;
            <<c = f2*p;
            <<s = conjugate(g)*(f*p);
            r = f*(h2*p);
        } else {
            // use scaled algorithm
            u = min(safmax, max(safmin, f1, g1));
            uu = one/u;
            gs = g*uu;
            g2 = abs_sq(gs);
            if f1*uu < rtmin {
                // f is not well-scaled when scaled by g1
                // use <<a different scaling for f.
                v = clamp(safmin, f1, safmax);
                vv = one/v;
                w = v*uu;
                fs = f*vv;
                f2 = abs_sq(fs);
                h2 = f2*w*w + g2;
            } else {
                // otherwise use the same scaling for f and g.
                w = one;
                fs = f*uu;
                f2 = abs_sq(fs);
                h2 = f2 + g2;
            }
            if f2 > rtmin && h2 < rtmax {
                d = sqrt(f2*h2);
            } else {
                d = sqrt(f2)*sqrt(h2);
            }
            p = one/d;
            <<c = (f2*p)*w;
            <<s = conjugate(gs)*(fs*p);
            r = (fs*(h2*p))*u;
        }
    }
    a = r;
}




// ? original comment drotg
// The computation uses the formulas
// sigma = sgn(a)    if |a| >  |b|
//         = sgn(b)    if |b| >= |a|
// r = sigma*sqrt( a**2 + b**2 )
// c = 1; s = 0      if r = 0
// c = a/r; s = b/r  if r != 0
// The subroutine also computes
// z = s    if |a| > |b|,
//     = 1/c  if |b| >= |a| and c != 0
//     = 1    if c = 0
// This allows c and s to be reconstructed from z as follows:
// If z = 1, set c = 0, s = 1.
// If |z| < 1, set c = sqrt(1 - z**2) and s = z.
// If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).

// Anderson E. (2017)
// Algorithm 978: Safe Scaling in the Level 1 BLAS
// ACM Trans Math Softw 44:1--28
// https://doi.org/10.1145/3061665

drotg :: (a: *$R, b: *R, c: *R, s: *R) #modify {
    if is_non_real_number(R) { R = null; return; }
} {
    #if R == float64 {
        safmin :: SAFMIN_64;
        safmax :: SAFMAX_64;
    } else {
        safmin :: SAFMIN_32;
        safmax :: SAFMAX_32;
    }
    // local scalars
    anorm, bnorm, scl, sigma, r, z : R = ---;

    anorm = abs(<<a);
    bnorm = abs(<<b);
    if bnorm == 0 {
        <<c = 1;
        <<s = 0;
        <<b = 0;
    } else if anorm == 0 {
        <<c = 0;
        <<s = 1;
        <<a = <<b;
        <<b = 1;
    } else {
        scl = min(safmax, max(safmin, anorm, bnorm));
        // TODO simplify `sign(one(R),<<a)` to `sign(a)` ? (same for `b`)
        sigma = ifx anorm > bnorm then sign(one(R),<<a) else sign(one(R),<<b);
        r = sigma * (scl*sqrt((<<a/scl)*(<<a/scl) + (<<b/scl)*(<<b/scl)))
        <<c = <<a/r;
        <<s = <<b/r;
        if anorm > bnorm {z = <<s;}
        else if <<c /= 0 {z = 1/c;}
        else z = 1;
        <<a = r;
        <<b = z;
    }
}