// ? LAPACK *laein()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaein/claein/dlaein/zlaein
// *LAEIN computes a specified right or left eigenvector of an upper
// Hessenberg matrix by inverse iteration.

slaein :: #bake_arguments _laein_real(RT=float32);
dlaein :: #bake_arguments _laein_real(RT=float64);
claein :: #bake_arguments _laein_complex(RT=float32, MT=C32);
zlaein :: #bake_arguments _laein_complex(RT=float64, MT=C64);

#scope_file
_laein_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laein_asum_real :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        return sasum(n, x, incx);
    } else {
        return dasum(n, x, incx);
    }
}

_laein_nrm2_real :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        return snrm2(n, x, incx);
    } else {
        return dnrm2(n, x, incx);
    }
}

_laein_scal_real :: inline ($$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, incx);
    } else {
        dscal(n, alpha, x, incx);
    }
}

_laein_amax_real_1based :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx) + 1;
    } else {
        return idamax(n, x, incx) + 1;
    }
}

_laein_lapy2 :: inline ($$x: $RT, $$y: RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slapy2(x, y);
    } else {
        return dlapy2(x, y);
    }
}

_laein_latrs_real :: inline (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$normin: u8, $$n: int,
    a: $A/ArrayType2D, $$lda: int, x: *$X/ArrayType1D, scale: *$RT,
    cnorm: *$CNORM/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           X.data_type == RT &&
           CNORM.data_type == RT;
} {
    #if RT == float32 {
        slatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    } else {
        dlatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    }
}

_laein_ladiv_real :: inline ($$a: $RT, $$b: RT, $$c: RT, $$d: RT, p: *RT, q: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        sladiv(a, b, c, d, p, q);
    } else {
        dladiv(a, b, c, d, p, q);
    }
}

_laein_asum_complex :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        return scsum1(n, x, incx);
    } else {
        return dzsum1(n, x, incx);
    }
}

_laein_nrm2_complex :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_laein_scal_complex :: inline ($$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_laein_amax_complex_1based :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        return icamax(n, x, incx) + 1;
    } else {
        return izamax(n, x, incx) + 1;
    }
}

_laein_abs1_complex :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return _abs1(x);
}

_laein_latrs_complex :: inline (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$normin: u8, $$n: int,
    a: $A/ArrayType2D, $$lda: int, x: *$X/ArrayType1D, scale: *$RT,
    cnorm: *$CNORM/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           X.data_type == MT &&
           CNORM.data_type == RT;
} {
    #if MT == C32 {
        clatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    } else {
        zlatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    }
}

_laein_ladiv_complex :: inline ($$x: $MT, $$y: MT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        return cladiv(x, y);
    } else {
        return zladiv(x, y);
    }
}

_laein_real :: (
    $$rightv: bool,
    $$noinit: bool,
    $$n: int,
    h: $H/ArrayType2D,
    $$ldh: int,
    $$wr: $RT,
    $$wi: RT,
    vr: *$VR/ArrayType1D,
    vi: *$VI/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    work: *$WORK/ArrayType1D,
    $$eps3: RT,
    $$smlnum: RT,
    $$bignum: RT,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT &&
           VR.data_type == RT &&
           VI.data_type == RT &&
           B.data_type == RT &&
           WORK.data_type == RT;
} {
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    tenth_r :: #run ncast(RT, 1.0e-1);

    i, i1, i2, i3, ierr, its, j : int = ---;
    absbii, absbjj, ei, ej, growto, norm, nrmsml, rec, rootn, scale, temp, vcrit, vmax, vnorm, w, w1, x, xi, xr, y : RT = ---;
    converged : bool = ---;
    normin, trans : u8 = ---;

    info.* = 0;

    rootn = sqrt(ncast(RT, n));
    growto = tenth_r / rootn;
    nrmsml = max(one_r, eps3*rootn)*smlnum;

    for j:1..n {
        for i:1..j - 1 {
            fset(b, i, j, fget(h, i, j));
        }
        fset(b, j, j, fget(h, j, j) - wr);
    }

    if wi == zero_r {
        if noinit {
            for i:1..n {
                fset(vr, i, eps3);
            }
        } else {
            vnorm = _laein_nrm2_real(n, vr.*, 1, RT);
            _laein_scal_real(n, (eps3*rootn) / max(vnorm, nrmsml), vr, 1, RT);
        }

        if rightv {
            for i:1..n - 1 {
                ei = fget(h, i + 1, i);
                if abs(fget(b.*, i, i)) < abs(ei) {
                    x = fget(b.*, i, i) / ei;
                    fset(b, i, i, ei);
                    for j:i + 1..n {
                        temp = fget(b.*, i + 1, j);
                        fset(b, i + 1, j, fget(b.*, i, j) - x*temp);
                        fset(b, i, j, temp);
                    }
                } else {
                    if fget(b.*, i, i) == zero_r {
                        fset(b, i, i, eps3);
                    }
                    x = ei / fget(b.*, i, i);
                    if x != zero_r {
                        for j:i + 1..n {
                            fset(b, i + 1, j, fget(b.*, i + 1, j) - x*fget(b.*, i, j));
                        }
                    }
                }
            }
            if fget(b.*, n, n) == zero_r {
                fset(b, n, n, eps3);
            }
            trans = #char "N";
        } else {
            for <j:n..2 {
                ej = fget(h, j, j - 1);
                if abs(fget(b.*, j, j)) < abs(ej) {
                    x = fget(b.*, j, j) / ej;
                    fset(b, j, j, ej);
                    for i:1..j - 1 {
                        temp = fget(b.*, i, j - 1);
                        fset(b, i, j - 1, fget(b.*, i, j) - x*temp);
                        fset(b, i, j, temp);
                    }
                } else {
                    if fget(b.*, j, j) == zero_r {
                        fset(b, j, j, eps3);
                    }
                    x = ej / fget(b.*, j, j);
                    if x != zero_r {
                        for i:1..j - 1 {
                            fset(b, i, j - 1, fget(b.*, i, j - 1) - x*fget(b.*, i, j));
                        }
                    }
                }
            }
            if fget(b.*, 1, 1) == zero_r {
                fset(b, 1, 1, eps3);
            }
            trans = #char "T";
        }

        normin = #char "N";
        converged = false;
        for its:1..n {
            _laein_latrs_real(#char "U", trans, #char "N", normin, n, b.*, ldb, vr, *scale, work, *ierr, RT);
            normin = #char "Y";

            vnorm = _laein_asum_real(n, vr.*, 1, RT);
            if vnorm >= growto*scale {
                converged = true;
                break;
            }

            temp = eps3 / (rootn + one_r);
            fset(vr, 1, eps3);
            for i:2..n {
                fset(vr, i, temp);
            }
            fset(vr, n - its + 1, fget(vr.*, n - its + 1) - eps3*rootn);
        }

        if !converged {
            info.* = 1;
        }

        i = _laein_amax_real_1based(n, vr.*, 1, RT);
        _laein_scal_real(n, one_r / abs(fget(vr.*, i)), vr, 1, RT);
    } else {
        if noinit {
            for i:1..n {
                fset(vr, i, eps3);
                fset(vi, i, zero_r);
            }
        } else {
            norm = _laein_lapy2(_laein_nrm2_real(n, vr.*, 1, RT), _laein_nrm2_real(n, vi.*, 1, RT), RT);
            rec = (eps3*rootn) / max(norm, nrmsml);
            _laein_scal_real(n, rec, vr, 1, RT);
            _laein_scal_real(n, rec, vi, 1, RT);
        }

        if rightv {
            fset(b, 2, 1, -wi);
            for i:2..n {
                fset(b, i + 1, 1, zero_r);
            }

            for i:1..n - 1 {
                absbii = _laein_lapy2(fget(b.*, i, i), fget(b.*, i + 1, i), RT);
                ei = fget(h, i + 1, i);
                if absbii < abs(ei) {
                    xr = fget(b.*, i, i) / ei;
                    xi = fget(b.*, i + 1, i) / ei;
                    fset(b, i, i, ei);
                    fset(b, i + 1, i, zero_r);
                    for j:i + 1..n {
                        temp = fget(b.*, i + 1, j);
                        fset(b, i + 1, j, fget(b.*, i, j) - xr*temp);
                        fset(b, j + 1, i + 1, fget(b.*, j + 1, i) - xi*temp);
                        fset(b, i, j, temp);
                        fset(b, j + 1, i, zero_r);
                    }
                    fset(b, i + 2, i, -wi);
                    fset(b, i + 1, i + 1, fget(b.*, i + 1, i + 1) - xi*wi);
                    fset(b, i + 2, i + 1, fget(b.*, i + 2, i + 1) + xr*wi);
                } else {
                    if absbii == zero_r {
                        fset(b, i, i, eps3);
                        fset(b, i + 1, i, zero_r);
                        absbii = eps3;
                    }
                    ei = (ei / absbii) / absbii;
                    xr = fget(b.*, i, i)*ei;
                    xi = -fget(b.*, i + 1, i)*ei;
                    for j:i + 1..n {
                        fset(b, i + 1, j,
                            fget(b.*, i + 1, j) - xr*fget(b.*, i, j) + xi*fget(b.*, j + 1, i));
                        fset(b, j + 1, i + 1,
                            -xr*fget(b.*, j + 1, i) - xi*fget(b.*, i, j));
                    }
                    fset(b, i + 2, i + 1, fget(b.*, i + 2, i + 1) - wi);
                }

                fset(work, i, zero_r);
                for j:i + 1..n {
                    fset(work, i, fget(work.*, i) + abs(fget(b.*, i, j)));
                    fset(work, i, fget(work.*, i) + abs(fget(b.*, j + 1, i)));
                }
            }
            if fget(b.*, n, n) == zero_r && fget(b.*, n + 1, n) == zero_r {
                fset(b, n, n, eps3);
            }
            fset(work, n, zero_r);

            i1 = n;
            i2 = 1;
            i3 = -1;
        } else {
            fset(b, n + 1, n, wi);
            for j:1..n - 1 {
                fset(b, n + 1, j, zero_r);
            }

            for <j:n..2 {
                ej = fget(h, j, j - 1);
                absbjj = _laein_lapy2(fget(b.*, j, j), fget(b.*, j + 1, j), RT);
                if absbjj < abs(ej) {
                    xr = fget(b.*, j, j) / ej;
                    xi = fget(b.*, j + 1, j) / ej;
                    fset(b, j, j, ej);
                    fset(b, j + 1, j, zero_r);
                    for i:1..j - 1 {
                        temp = fget(b.*, i, j - 1);
                        fset(b, i, j - 1, fget(b.*, i, j) - xr*temp);
                        fset(b, j, i, fget(b.*, j + 1, i) - xi*temp);
                        fset(b, i, j, temp);
                        fset(b, j + 1, i, zero_r);
                    }
                    fset(b, j + 1, j - 1, wi);
                    fset(b, j - 1, j - 1, fget(b.*, j - 1, j - 1) + xi*wi);
                    fset(b, j, j - 1, fget(b.*, j, j - 1) - xr*wi);
                } else {
                    if absbjj == zero_r {
                        fset(b, j, j, eps3);
                        fset(b, j + 1, j, zero_r);
                        absbjj = eps3;
                    }
                    ej = (ej / absbjj) / absbjj;
                    xr = fget(b.*, j, j)*ej;
                    xi = -fget(b.*, j + 1, j)*ej;
                    for i:1..j - 1 {
                        fset(b, i, j - 1,
                            fget(b.*, i, j - 1) - xr*fget(b.*, i, j) + xi*fget(b.*, j + 1, i));
                        fset(b, j, i,
                            -xr*fget(b.*, j + 1, i) - xi*fget(b.*, i, j));
                    }
                    fset(b, j, j - 1, fget(b.*, j, j - 1) + wi);
                }

                fset(work, j, zero_r);
                for i:1..j - 1 {
                    fset(work, j, fget(work.*, j) + abs(fget(b.*, i, j)));
                    fset(work, j, fget(work.*, j) + abs(fget(b.*, j + 1, i)));
                }
            }
            if fget(b.*, 1, 1) == zero_r && fget(b.*, 2, 1) == zero_r {
                fset(b, 1, 1, eps3);
            }
            fset(work, 1, zero_r);

            i1 = 1;
            i2 = n;
            i3 = 1;
        }

        converged = false;
        for its:1..n {
            scale = one_r;
            vmax = one_r;
            vcrit = bignum;

            i = i1;
            while true {
                if fget(work.*, i) > vcrit {
                    rec = one_r / vmax;
                    _laein_scal_real(n, rec, vr, 1, RT);
                    _laein_scal_real(n, rec, vi, 1, RT);
                    scale = scale*rec;
                    vmax = one_r;
                    vcrit = bignum;
                }

                xr = fget(vr.*, i);
                xi = fget(vi.*, i);
                if rightv {
                    for j:i + 1..n {
                        xr = xr - fget(b.*, i, j)*fget(vr.*, j) + fget(b.*, j + 1, i)*fget(vi.*, j);
                        xi = xi - fget(b.*, i, j)*fget(vi.*, j) - fget(b.*, j + 1, i)*fget(vr.*, j);
                    }
                } else {
                    for j:1..i - 1 {
                        xr = xr - fget(b.*, j, i)*fget(vr.*, j) + fget(b.*, i + 1, j)*fget(vi.*, j);
                        xi = xi - fget(b.*, j, i)*fget(vi.*, j) - fget(b.*, i + 1, j)*fget(vr.*, j);
                    }
                }

                w = abs(fget(b.*, i, i)) + abs(fget(b.*, i + 1, i));
                if w > smlnum {
                    if w < one_r {
                        w1 = abs(xr) + abs(xi);
                        if w1 > w*bignum {
                            rec = one_r / w1;
                            _laein_scal_real(n, rec, vr, 1, RT);
                            _laein_scal_real(n, rec, vi, 1, RT);
                            xr = fget(vr.*, i);
                            xi = fget(vi.*, i);
                            scale = scale*rec;
                            vmax = vmax*rec;
                        }
                    }

                    _laein_ladiv_real(xr, xi, fget(b.*, i, i), fget(b.*, i + 1, i), *xr, *xi, RT);
                    fset(vr, i, xr);
                    fset(vi, i, xi);
                    vmax = max(abs(xr) + abs(xi), vmax);
                    vcrit = bignum / vmax;
                } else {
                    for j:1..n {
                        fset(vr, j, zero_r);
                        fset(vi, j, zero_r);
                    }
                    fset(vr, i, one_r);
                    fset(vi, i, one_r);
                    scale = zero_r;
                    vmax = one_r;
                    vcrit = bignum;
                }

                if i == i2 break;
                i += i3;
            }

            vnorm = _laein_asum_real(n, vr.*, 1, RT) + _laein_asum_real(n, vi.*, 1, RT);
            if vnorm >= growto*scale {
                converged = true;
                break;
            }

            y = eps3 / (rootn + one_r);
            fset(vr, 1, eps3);
            fset(vi, 1, zero_r);
            for i:2..n {
                fset(vr, i, y);
                fset(vi, i, zero_r);
            }
            fset(vr, n - its + 1, fget(vr.*, n - its + 1) - eps3*rootn);
        }

        if !converged {
            info.* = 1;
        }

        vnorm = zero_r;
        for i:1..n {
            vnorm = max(vnorm, abs(fget(vr.*, i)) + abs(fget(vi.*, i)));
        }
        _laein_scal_real(n, one_r / vnorm, vr, 1, RT);
        _laein_scal_real(n, one_r / vnorm, vi, 1, RT);
    }
}

_laein_complex :: (
    $$rightv: bool,
    $$noinit: bool,
    $$n: int,
    h: $H/ArrayType2D,
    $$ldh: int,
    $$w: $MT,
    v: *$V/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    rwork: *$RWORK/ArrayType1D,
    $$eps3: $RT,
    $$smlnum: RT,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT &&
           V.data_type == MT &&
           B.data_type == MT &&
           RWORK.data_type == RT;
} {
    one_r   :: #run one(RT);
    tenth_r :: #run ncast(RT, 1.0e-1);
    zero_t  :: #run zero(MT);

    i, ierr, its, j : int = ---;
    growto, nrmsml, rootn, rtemp, scale, vnorm : RT = ---;
    converged : bool = ---;
    ei, ej, temp, x : MT = ---;
    normin, trans : u8 = ---;

    info.* = 0;

    rootn = sqrt(ncast(RT, n));
    growto = tenth_r / rootn;
    nrmsml = max(one_r, eps3*rootn)*smlnum;

    for j:1..n {
        for i:1..j - 1 {
            fset(b, i, j, fget(h, i, j));
        }
        fset(b, j, j, fget(h, j, j) - w);
    }

    if noinit {
        for i:1..n {
            fset(v, i, _laein_complex_from_real(eps3, RT, MT));
        }
    } else {
        vnorm = _laein_nrm2_complex(n, v.*, 1, RT, MT);
        _laein_scal_complex(n, (eps3*rootn) / max(vnorm, nrmsml), v, 1, RT, MT);
    }

    if rightv {
        for i:1..n - 1 {
            ei = fget(h, i + 1, i);
            if _laein_abs1_complex(fget(b.*, i, i), RT, MT) < _laein_abs1_complex(ei, RT, MT) {
                x = _laein_ladiv_complex(fget(b.*, i, i), ei, RT, MT);
                fset(b, i, i, ei);
                for j:i + 1..n {
                    temp = fget(b.*, i + 1, j);
                    fset(b, i + 1, j, fget(b.*, i, j) - x*temp);
                    fset(b, i, j, temp);
                }
            } else {
                if fget(b.*, i, i) == zero_t {
                    fset(b, i, i, _laein_complex_from_real(eps3, RT, MT));
                }
                x = _laein_ladiv_complex(ei, fget(b.*, i, i), RT, MT);
                if x != zero_t {
                    for j:i + 1..n {
                        fset(b, i + 1, j, fget(b.*, i + 1, j) - x*fget(b.*, i, j));
                    }
                }
            }
        }
        if fget(b.*, n, n) == zero_t {
            fset(b, n, n, _laein_complex_from_real(eps3, RT, MT));
        }
        trans = #char "N";
    } else {
        for <j:n..2 {
            ej = fget(h, j, j - 1);
            if _laein_abs1_complex(fget(b.*, j, j), RT, MT) < _laein_abs1_complex(ej, RT, MT) {
                x = _laein_ladiv_complex(fget(b.*, j, j), ej, RT, MT);
                fset(b, j, j, ej);
                for i:1..j - 1 {
                    temp = fget(b.*, i, j - 1);
                    fset(b, i, j - 1, fget(b.*, i, j) - x*temp);
                    fset(b, i, j, temp);
                }
            } else {
                if fget(b.*, j, j) == zero_t {
                    fset(b, j, j, _laein_complex_from_real(eps3, RT, MT));
                }
                x = _laein_ladiv_complex(ej, fget(b.*, j, j), RT, MT);
                if x != zero_t {
                    for i:1..j - 1 {
                        fset(b, i, j - 1, fget(b.*, i, j - 1) - x*fget(b.*, i, j));
                    }
                }
            }
        }
        if fget(b.*, 1, 1) == zero_t {
            fset(b, 1, 1, _laein_complex_from_real(eps3, RT, MT));
        }
        trans = #char "C";
    }

    normin = #char "N";
    converged = false;
    for its:1..n {
        _laein_latrs_complex(#char "U", trans, #char "N", normin, n, b.*, ldb, v, *scale, rwork, *ierr, RT, MT);
        normin = #char "Y";

        vnorm = _laein_asum_complex(n, v.*, 1, RT, MT);
        if vnorm >= growto*scale {
            converged = true;
            break;
        }

        rtemp = eps3 / (rootn + one_r);
        fset(v, 1, _laein_complex_from_real(eps3, RT, MT));
        for i:2..n {
            fset(v, i, _laein_complex_from_real(rtemp, RT, MT));
        }
        fset(v, n - its + 1, fget(v.*, n - its + 1) - _laein_complex_from_real(eps3*rootn, RT, MT));
    }

    if !converged {
        info.* = 1;
    }

    i = _laein_amax_complex_1based(n, v.*, 1, RT, MT);
    _laein_scal_complex(n, one_r / _laein_abs1_complex(fget(v.*, i), RT, MT), v, 1, RT, MT);
}
