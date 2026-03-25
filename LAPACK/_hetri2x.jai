// ? LAPACK *hetri2x()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetri2x/zhetri2x
// *HETRI2X computes the inverse of a complex Hermitian indefinite matrix
// using the blocked workspace kernel behind *HETRI2.

chetri2x :: #bake_arguments _hetri2x(RT=float32, MT=C32);
zhetri2x :: #bake_arguments _hetri2x(RT=float64, MT=C64);

#scope_file
_hetri2x_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetri2x_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetri2x_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetri2x_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetri2x_syconv :: (
    $$uplo: u8, $$way: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D, iinfo: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 {
        csyconv(uplo, way, n, a, lda, ipiv, work, iinfo);
    } else {
        zsyconv(uplo, way, n, a, lda, ipiv, work, iinfo);
    }
}

_hetri2x_trtri :: (
    $$uplo: u8, $$diag: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 {
        ctrtri(uplo, diag, n, a, lda, info);
    } else {
        ztrtri(uplo, diag, n, a, lda, info);
    }
}

_hetri2x_trmm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hetri2x_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int, alpha: $MT,
    a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int, beta: $MT,
    c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_hetri2x_swapr :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, $$i1: int, $$i2: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && A.data_type == MT;
} {
    #if MT == C32 {
        cheswapr(uplo, n, a, lda, i1, i2);
    } else {
        zheswapr(uplo, n, a, lda, i1, i2);
    }
}

_hetri2x :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$nb: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    one_r  :: #run one(RT);
    one_t  :: #run _hetri2x_complex_from_real(one_r, RT, MT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == C32 { "CHETRI2X" } else { "ZHETRI2X" };

    upper : bool = ---;
    count, cut, i, iinfo, invd, ip, j, k, nnb, u11 : int = ---;
    ak, akp1, t : RT = ---;
    akkp1, d_c, u01_i_j, u01_ip1_j, u11_i_j, u11_ip1_j : MT = ---;

    ldwork := n + nb + 1;
    work_all := array2d_view(_hetri2x_raw_slice(work.data, ldwork*(nb + 3)), ldwork, row_major=false);

    info.* = 0;
    upper = lsame(uplo, #char "U");

    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if nb < 1 {
        info.* = -7;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if n == 0 then return;

    _hetri2x_syconv(uplo, #char "C", n, a, lda, ipiv, work, *iinfo, MT);

    if upper {
        for <info.*:1..n {
            if fget(ipiv, info.*) > 0 && fget(a, info.*, info.*) == zero_t {
                return;
            }
        }
    } else {
        for info.*:1..n {
            if fget(ipiv, info.*) > 0 && fget(a, info.*, info.*) == zero_t {
                return;
            }
        }
    }
    info.* = 0;

    u11 = n;
    invd = nb + 2;

    if upper {
        _hetri2x_trtri(uplo, #char "U", n, a, lda, info, MT);

        k = 1;
        while k <= n {
            if fget(ipiv, k) > 0 {
                fset(work_all, k, invd, _hetri2x_complex_from_real(one_r / real(fget(a.*, k, k)), RT, MT));
                fset(work_all, k, invd + 1, zero_t);
                k += 1;
            } else {
                t = abs(fget(work_all, k + 1, 1));
                ak = real(fget(a.*, k, k)) / t;
                akp1 = real(fget(a.*, k + 1, k + 1)) / t;
                akkp1 = fget(work_all, k + 1, 1) / t;
                d_c = _hetri2x_complex_from_real(t*(ak*akp1 - one_r), RT, MT);
                fset(work_all, k, invd, _hetri2x_complex_from_real(akp1, RT, MT) / d_c);
                fset(work_all, k + 1, invd + 1, _hetri2x_complex_from_real(ak, RT, MT) / d_c);
                fset(work_all, k, invd + 1, -akkp1 / d_c);
                fset(work_all, k + 1, invd, conjugate(fget(work_all, k, invd + 1)));
                k += 2;
            }
        }

        cut = n;
        while cut > 0 {
            nnb = nb;
            if cut <= nnb {
                nnb = cut;
            } else {
                count = 0;
                for i:cut + 1 - nnb..cut {
                    if fget(ipiv, i) < 0 then count += 1;
                }
                if mod(count, 2) == 1 then nnb += 1;
            }

            cut -= nnb;

            for i:1..cut {
                for j:1..nnb {
                    fset(work_all, i, j, fget(a.*, i, cut + j));
                }
            }

            for i:1..nnb {
                fset(work_all, u11 + i, i, one_t);
                for j:1..i - 1 {
                    fset(work_all, u11 + i, j, zero_t);
                }
                for j:i + 1..nnb {
                    fset(work_all, u11 + i, j, fget(a.*, cut + i, cut + j));
                }
            }

            i = 1;
            while i <= cut {
                if fget(ipiv, i) > 0 {
                    for j:1..nnb {
                        fset(work_all, i, j, fget(work_all, i, invd) * fget(work_all, i, j));
                    }
                    i += 1;
                } else {
                    for j:1..nnb {
                        u01_i_j = fget(work_all, i, j);
                        u01_ip1_j = fget(work_all, i + 1, j);
                        fset(work_all, i, j, fget(work_all, i, invd) * u01_i_j + fget(work_all, i, invd + 1) * u01_ip1_j);
                        fset(work_all, i + 1, j, fget(work_all, i + 1, invd) * u01_i_j + fget(work_all, i + 1, invd + 1) * u01_ip1_j);
                    }
                    i += 2;
                }
            }

            i = 1;
            while i <= nnb {
                if fget(ipiv, cut + i) > 0 {
                    for j:i..nnb {
                        fset(work_all, u11 + i, j, fget(work_all, cut + i, invd) * fget(work_all, u11 + i, j));
                    }
                    i += 1;
                } else {
                    for j:i..nnb {
                        u11_i_j = fget(work_all, u11 + i, j);
                        u11_ip1_j = fget(work_all, u11 + i + 1, j);
                        fset(work_all, u11 + i, j, fget(work_all, cut + i, invd) * fget(work_all, u11 + i, j) + fget(work_all, cut + i, invd + 1) * fget(work_all, u11 + i + 1, j));
                        fset(work_all, u11 + i + 1, j, fget(work_all, cut + i + 1, invd) * u11_i_j + fget(work_all, cut + i + 1, invd + 1) * u11_ip1_j);
                    }
                    i += 2;
                }
            }

            a11_raw := _hetri2x_matrix_storage(a.*, lda, cut + 1, cut + 1, nnb, nnb);
            a11 := array2d_view(a11_raw, lda, row_major=false);
            w11_raw := _hetri2x_matrix_storage(work_all.*, ldwork, u11 + 1, 1, nnb, nnb);
            w11 := array2d_view(w11_raw, ldwork, row_major=false);
            _hetri2x_trmm(#char "L", #char "U", #char "C", #char "U", nnb, nnb, one_t, a11.*, lda, *w11, ldwork, MT);

            for i:1..nnb {
                for j:i..nnb {
                    fset(a, cut + i, cut + j, fget(work_all, u11 + i, j));
                }
            }

            if cut > 0 {
                a01_raw := _hetri2x_matrix_storage(a.*, lda, 1, cut + 1, cut, nnb);
                a01 := array2d_view(a01_raw, lda, row_major=false);
                w01_raw := _hetri2x_matrix_storage(work_all.*, ldwork, 1, 1, cut, nnb);
                w01 := array2d_view(w01_raw, ldwork, row_major=false);
                w11b_raw := _hetri2x_matrix_storage(work_all.*, ldwork, u11 + 1, 1, nnb, nnb);
                w11b := array2d_view(w11b_raw, ldwork, row_major=false);
                _hetri2x_gemm(#char "C", #char "N", nnb, nnb, cut, one_t, a01.*, lda, w01.*, ldwork, zero_t, *w11b, ldwork, MT);

                for i:1..nnb {
                    for j:i..nnb {
                        fset(a, cut + i, cut + j, fget(a.*, cut + i, cut + j) + fget(work_all, u11 + i, j));
                    }
                }

                a00_raw := _hetri2x_matrix_storage(a.*, lda, 1, 1, cut, cut);
                a00 := array2d_view(a00_raw, lda, row_major=false);
                _hetri2x_trmm(#char "L", uplo, #char "C", #char "U", cut, nnb, one_t, a00.*, lda, *w01, ldwork, MT);

                for i:1..cut {
                    for j:1..nnb {
                        fset(a, i, cut + j, fget(work_all, i, j));
                    }
                }
            }
        }

        i = 1;
        while i <= n {
            if fget(ipiv, i) > 0 {
                ip = fget(ipiv, i);
                if i < ip {
                    _hetri2x_swapr(uplo, n, a, lda, i, ip, RT, MT);
                } else if i > ip {
                    _hetri2x_swapr(uplo, n, a, lda, ip, i, RT, MT);
                }
            } else {
                ip = -fget(ipiv, i);
                i += 1;
                if i - 1 < ip {
                    _hetri2x_swapr(uplo, n, a, lda, i - 1, ip, RT, MT);
                } else if i - 1 > ip {
                    _hetri2x_swapr(uplo, n, a, lda, ip, i - 1, RT, MT);
                }
            }
            i += 1;
        }
    } else {
        _hetri2x_trtri(uplo, #char "U", n, a, lda, info, MT);

        k = n;
        while k >= 1 {
            if fget(ipiv, k) > 0 {
                fset(work_all, k, invd, _hetri2x_complex_from_real(one_r / real(fget(a.*, k, k)), RT, MT));
                fset(work_all, k, invd + 1, zero_t);
                k -= 1;
            } else {
                t = abs(fget(work_all, k - 1, 1));
                ak = real(fget(a.*, k - 1, k - 1)) / t;
                akp1 = real(fget(a.*, k, k)) / t;
                akkp1 = fget(work_all, k - 1, 1) / t;
                d_c = _hetri2x_complex_from_real(t*(ak*akp1 - one_r), RT, MT);
                fset(work_all, k - 1, invd, _hetri2x_complex_from_real(akp1, RT, MT) / d_c);
                fset(work_all, k, invd, _hetri2x_complex_from_real(ak, RT, MT) / d_c);
                fset(work_all, k, invd + 1, -akkp1 / d_c);
                fset(work_all, k - 1, invd + 1, conjugate(fget(work_all, k, invd + 1)));
                k -= 2;
            }
        }

        cut = 0;
        while cut < n {
            nnb = nb;
            if cut + nnb >= n {
                nnb = n - cut;
            } else {
                count = 0;
                for i:cut + 1..cut + nnb {
                    if fget(ipiv, i) < 0 then count += 1;
                }
                if mod(count, 2) == 1 then nnb += 1;
            }

            for i:1..n - cut - nnb {
                for j:1..nnb {
                    fset(work_all, i, j, fget(a.*, cut + nnb + i, cut + j));
                }
            }

            for i:1..nnb {
                fset(work_all, u11 + i, i, one_t);
                for j:i + 1..nnb {
                    fset(work_all, u11 + i, j, zero_t);
                }
                for j:1..i - 1 {
                    fset(work_all, u11 + i, j, fget(a.*, cut + i, cut + j));
                }
            }

            i = n - cut - nnb;
            while i >= 1 {
                if fget(ipiv, cut + nnb + i) > 0 {
                    for j:1..nnb {
                        fset(work_all, i, j, fget(work_all, cut + nnb + i, invd) * fget(work_all, i, j));
                    }
                    i -= 1;
                } else {
                    for j:1..nnb {
                        u01_i_j = fget(work_all, i, j);
                        u01_ip1_j = fget(work_all, i - 1, j);
                        fset(work_all, i, j, fget(work_all, cut + nnb + i, invd) * u01_i_j + fget(work_all, cut + nnb + i, invd + 1) * u01_ip1_j);
                        fset(work_all, i - 1, j, fget(work_all, cut + nnb + i - 1, invd + 1) * u01_i_j + fget(work_all, cut + nnb + i - 1, invd) * u01_ip1_j);
                    }
                    i -= 2;
                }
            }

            i = nnb;
            while i >= 1 {
                if fget(ipiv, cut + i) > 0 {
                    for j:1..nnb {
                        fset(work_all, u11 + i, j, fget(work_all, cut + i, invd) * fget(work_all, u11 + i, j));
                    }
                    i -= 1;
                } else {
                    for j:1..nnb {
                        u11_i_j = fget(work_all, u11 + i, j);
                        u11_ip1_j = fget(work_all, u11 + i - 1, j);
                        fset(work_all, u11 + i, j, fget(work_all, cut + i, invd) * u11_i_j + fget(work_all, cut + i, invd + 1) * u11_ip1_j);
                        fset(work_all, u11 + i - 1, j, fget(work_all, cut + i - 1, invd + 1) * u11_i_j + fget(work_all, cut + i - 1, invd) * u11_ip1_j);
                    }
                    i -= 2;
                }
            }

            a11_raw := _hetri2x_matrix_storage(a.*, lda, cut + 1, cut + 1, nnb, nnb);
            a11 := array2d_view(a11_raw, lda, row_major=false);
            w11_raw := _hetri2x_matrix_storage(work_all.*, ldwork, u11 + 1, 1, nnb, nnb);
            w11 := array2d_view(w11_raw, ldwork, row_major=false);
            _hetri2x_trmm(#char "L", uplo, #char "C", #char "U", nnb, nnb, one_t, a11.*, lda, *w11, ldwork, MT);

            for i:1..nnb {
                for j:1..i {
                    fset(a, cut + i, cut + j, fget(work_all, u11 + i, j));
                }
            }

            if cut + nnb < n {
                a21_raw := _hetri2x_matrix_storage(a.*, lda, cut + nnb + 1, cut + 1, n - cut - nnb, nnb);
                a21 := array2d_view(a21_raw, lda, row_major=false);
                w01_raw := _hetri2x_matrix_storage(work_all.*, ldwork, 1, 1, n - cut - nnb, nnb);
                w01 := array2d_view(w01_raw, ldwork, row_major=false);
                w11b_raw := _hetri2x_matrix_storage(work_all.*, ldwork, u11 + 1, 1, nnb, nnb);
                w11b := array2d_view(w11b_raw, ldwork, row_major=false);
                _hetri2x_gemm(#char "C", #char "N", nnb, nnb, n - nnb - cut, one_t, a21.*, lda, w01.*, ldwork, zero_t, *w11b, ldwork, MT);

                for i:1..nnb {
                    for j:1..i {
                        fset(a, cut + i, cut + j, fget(a.*, cut + i, cut + j) + fget(work_all, u11 + i, j));
                    }
                }

                a22_raw := _hetri2x_matrix_storage(a.*, lda, cut + nnb + 1, cut + nnb + 1, n - cut - nnb, n - cut - nnb);
                a22 := array2d_view(a22_raw, lda, row_major=false);
                _hetri2x_trmm(#char "L", uplo, #char "C", #char "U", n - cut - nnb, nnb, one_t, a22.*, lda, *w01, ldwork, MT);

                for i:1..n - cut - nnb {
                    for j:1..nnb {
                        fset(a, cut + nnb + i, cut + j, fget(work_all, i, j));
                    }
                }
            } else {
                for i:1..nnb {
                    for j:1..i {
                        fset(a, cut + i, cut + j, fget(work_all, u11 + i, j));
                    }
                }
            }

            cut += nnb;
        }

        i = n;
        while i >= 1 {
            if fget(ipiv, i) > 0 {
                ip = fget(ipiv, i);
                if i < ip {
                    _hetri2x_swapr(uplo, n, a, lda, i, ip, RT, MT);
                } else if i > ip {
                    _hetri2x_swapr(uplo, n, a, lda, ip, i, RT, MT);
                }
            } else {
                ip = -fget(ipiv, i);
                if i < ip {
                    _hetri2x_swapr(uplo, n, a, lda, i, ip, RT, MT);
                } else if i > ip {
                    _hetri2x_swapr(uplo, n, a, lda, ip, i, RT, MT);
                }
                i -= 1;
            }
            i -= 1;
        }
    }
}
