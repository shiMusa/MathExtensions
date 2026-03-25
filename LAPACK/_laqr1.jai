// ? LAPACK *laqr1()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqr1/dlaqr1/claqr1/zlaqr1
// *LAQR1 forms the first column of a shifted 2x2 or 3x3 QR starter.

slaqr1 :: #bake_arguments _laqr1_real(RT=float32);
dlaqr1 :: #bake_arguments _laqr1_real(RT=float64);
claqr1 :: #bake_arguments _laqr1_complex(RT=float32, MT=C32);
zlaqr1 :: #bake_arguments _laqr1_complex(RT=float64, MT=C64);

#scope_file
_laqr1_cabs1 :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return abs(real(x)) + abs(imag(x));
}

_laqr1_real :: (
    $$n: int,
    h: $H/ArrayType2D,
    $$ldh: int,
    $$sr1: $RT,
    $$si1: RT,
    $$sr2: RT,
    $$si2: RT,
    v: *$V/ArrayType1D,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT &&
           V.data_type == RT;
} {
    zero_r :: #run zero(RT);

    s, h21s, h31s : RT = ---;

    if n != 2 && n != 3 then return;

    if n == 2 {
        s = abs(fget(h, 1, 1) - sr2) + abs(si2) + abs(fget(h, 2, 1));
        if s == zero_r {
            fset(v, 1, zero_r);
            fset(v, 2, zero_r);
        } else {
            h21s = fget(h, 2, 1) / s;
            fset(v, 1, h21s*fget(h, 1, 2) + (fget(h, 1, 1) - sr1)*((fget(h, 1, 1) - sr2) / s) - si1*(si2 / s));
            fset(v, 2, h21s*(fget(h, 1, 1) + fget(h, 2, 2) - sr1 - sr2));
        }
    } else {
        s = abs(fget(h, 1, 1) - sr2) + abs(si2) + abs(fget(h, 2, 1)) + abs(fget(h, 3, 1));
        if s == zero_r {
            fset(v, 1, zero_r);
            fset(v, 2, zero_r);
            fset(v, 3, zero_r);
        } else {
            h21s = fget(h, 2, 1) / s;
            h31s = fget(h, 3, 1) / s;
            fset(v, 1,
                (fget(h, 1, 1) - sr1)*((fget(h, 1, 1) - sr2) / s) -
                si1*(si2 / s) +
                fget(h, 1, 2)*h21s +
                fget(h, 1, 3)*h31s
            );
            fset(v, 2, h21s*(fget(h, 1, 1) + fget(h, 2, 2) - sr1 - sr2) + fget(h, 2, 3)*h31s);
            fset(v, 3, h31s*(fget(h, 1, 1) + fget(h, 3, 3) - sr1 - sr2) + h21s*fget(h, 3, 2));
        }
    }
}

_laqr1_complex :: (
    $$n: int,
    h: $H/ArrayType2D,
    $$ldh: int,
    s1: $MT,
    s2: MT,
    v: *$V/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT &&
           V.data_type == MT;
} {
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);

    s : RT = ---;
    h21s, h31s : MT = ---;

    if n != 2 && n != 3 then return;

    if n == 2 {
        s = _laqr1_cabs1(fget(h, 1, 1) - s2, RT, MT) + _laqr1_cabs1(fget(h, 2, 1), RT, MT);
        if s == zero_r {
            fset(v, 1, zero_t);
            fset(v, 2, zero_t);
        } else {
            h21s = fget(h, 2, 1) / s;
            fset(v, 1, h21s*fget(h, 1, 2) + (fget(h, 1, 1) - s1)*((fget(h, 1, 1) - s2) / s));
            fset(v, 2, h21s*(fget(h, 1, 1) + fget(h, 2, 2) - s1 - s2));
        }
    } else {
        s = _laqr1_cabs1(fget(h, 1, 1) - s2, RT, MT) + _laqr1_cabs1(fget(h, 2, 1), RT, MT) + _laqr1_cabs1(fget(h, 3, 1), RT, MT);
        if s == zero_r {
            fset(v, 1, zero_t);
            fset(v, 2, zero_t);
            fset(v, 3, zero_t);
        } else {
            h21s = fget(h, 2, 1) / s;
            h31s = fget(h, 3, 1) / s;
            fset(v, 1, (fget(h, 1, 1) - s1)*((fget(h, 1, 1) - s2) / s) + fget(h, 1, 2)*h21s + fget(h, 1, 3)*h31s);
            fset(v, 2, h21s*(fget(h, 1, 1) + fget(h, 2, 2) - s1 - s2) + fget(h, 2, 3)*h31s);
            fset(v, 3, h31s*(fget(h, 1, 1) + fget(h, 3, 3) - s1 - s2) + h21s*fget(h, 3, 2));
        }
    }
}
