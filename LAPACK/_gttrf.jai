// ? LAPACK *gttrf()
// ? http://www.netlib.org/lapack/explore-html/
//
// ? original comment sgttrf/cgttrf/dgttrf/zgttrf
// *GTTRF computes an LU factorization of a real tridiagonal matrix A
// using elimination with partial pivoting and row interchanges.
//
sgttrf :: #bake_arguments _gttrf(RT=float32, MT=float32);
dgttrf :: #bake_arguments _gttrf(RT=float64, MT=float64);
cgttrf :: #bake_arguments _gttrf(RT=float32, MT=C32);
zgttrf :: #bake_arguments _gttrf(RT=float64, MT=C64);

#scope_file
_gttrf_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    // Match LAPACK's CABS1 definition for complex:
    // abs(real(z)) + abs(imag(z))
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return abs(real(z)) + abs(imag(z));
    }
}

_gttrf :: (
    $$n: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    du2: *$DU2/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int;
} {
    routine_name :: #run if MT == float32 {
        "SGTTRF"
    } else if MT == float64 {
        "DGTTRF"
    } else if MT == C32 {
        "CGTTRF"
    } else {
        "ZGTTRF"
    };

    zero_t :: #run zero(MT);

    info.* = 0;
    if n < 0 {
        info.* = -1;
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    // Initialize IPIV(i) = i and DU2(i) = 0.
    for i:1..n {
        fset(ipiv, i, i);
    }
    if n > 2 {
        for i:1..n-2 {
            fset(du2, i, zero_t);
        }
    }

    // Main elimination steps for i = 1..n-2
    if n > 2 {
        for i:1..n-2 {
            abs_d := _gttrf_abs1(fget(d.*, i), RT, MT);
            abs_dl := _gttrf_abs1(fget(dl.*, i), RT, MT);

            if abs_d >= abs_dl {
                // No row interchange required, eliminate DL(i)
                if fget(d.*, i) != zero_t {
                    fact : MT = fget(dl.*, i) / fget(d.*, i);
                    fset(dl, i, fact);
                    fset(d, i + 1, fget(d.*, i + 1) - fact*fget(du.*, i));
                }
            } else {
                // Interchange rows i and i+1, eliminate DL(i)
                fact : MT = fget(d.*, i) / fget(dl.*, i);

                dl_i_old   : MT = fget(dl.*, i);
                d_i1_old   : MT = fget(d.*, i + 1);
                du_i_old   : MT = fget(du.*, i);
                du_i1_old  : MT = fget(du.*, i + 1);

                fset(d, i, dl_i_old);
                fset(dl, i, fact);

                fset(du, i, d_i1_old);
                fset(d, i + 1, du_i_old - fact*d_i1_old);

                fset(du2, i, du_i1_old);
                fset(du, i + 1, -fact*du_i1_old);
                fset(ipiv, i, i + 1);
            }
        }
    }

    // Last elimination step i = n-1 (if n > 1)
    if n > 1 {
        i_last := n - 1;

        abs_d := _gttrf_abs1(fget(d.*, i_last), RT, MT);
        abs_dl := _gttrf_abs1(fget(dl.*, i_last), RT, MT);

        if abs_d >= abs_dl {
            if fget(d.*, i_last) != zero_t {
                fact : MT = fget(dl.*, i_last) / fget(d.*, i_last);
                fset(dl, i_last, fact);
                fset(d, i_last + 1, fget(d.*, i_last + 1) - fact*fget(du.*, i_last));
            }
        } else {
            fact : MT = fget(d.*, i_last) / fget(dl.*, i_last);

            dl_i_old  : MT = fget(dl.*, i_last);
            d_i1_old  : MT = fget(d.*, i_last + 1);
            du_i_old  : MT = fget(du.*, i_last);

            fset(d, i_last, dl_i_old);
            fset(dl, i_last, fact);

            fset(du, i_last, d_i1_old);
            fset(d, i_last + 1, du_i_old - fact*d_i1_old);
            fset(ipiv, i_last, i_last + 1);
        }
    }

    // Check for a zero on the diagonal of U.
    for i:1..n {
        if fget(d.*, i) == zero_t {
            info.* = i;
            break;
        }
    }
}

