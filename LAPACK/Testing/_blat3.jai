// ? LAPACK BLAT3 generic driver

#load "_blat3_checks.jai";

_blat3 :: ($T: Type, $U: Type, file_url: string, label: string) -> bool {
    #run {
        assert((T == float32 && U == float32) || (T == float64 && U == float64) || (T == C32 && U == float32) || (T == C64 && U == float64));
    }

    println_push(label, color=.FG_GREEN);

    data := read_file_and_split(file_url);

    println("Input file:");
    print(ansi_color(.FG_WHITE));
    for line, li: data {
        print("%: ", li);
        for line {
            print("% ", it);
        }
        print("\n");
    }
    print(ansi_color());

    nmax   : int : 65;
    nidmax : int : 9;
    nalmax : int : 7;
    nbemax : int : 7;

    complex_suite :: #run is_non_real_number(T);
    #if T == float32 || T == float64 {
        nsubs :: 7;
    } else {
        nsubs :: 10;
    }

    a  : Array2D(nmax, nmax, T);
    b  : Array2D(nmax, nmax, T);
    c  : Array2D(nmax, nmax, T);
    aa : Array1D(nmax*nmax, T);
    as : Array1D(nmax*nmax, T);
    bb : Array1D(nmax*nmax, T);
    bs : Array1D(nmax*nmax, T);
    cc : Array1D(nmax*nmax, T);
    cs : Array1D(nmax*nmax, T);
    ct : Array1D(nmax, T);
    g  : Array1D(nmax, U);
    w  : Array1D(2*nmax, T);

    thresh : U;
    nidim, nalf, nbet : int;
    idim : []int;
    alf  : Array1D(nalmax, T);
    bet  : Array1D(nbemax, T);
    ltest : [nsubs]bool;
    snames := NewArray(nsubs, string, false);
    #if T == float32 {
        snames[0] = "SGEMM";  snames[1] = "SSYMM";  snames[2] = "STRMM";
        snames[3] = "STRSM";  snames[4] = "SSYRK";  snames[5] = "SSYR2K";
        snames[6] = "SGEMMTR";
    }
    #if T == float64 {
        snames[0] = "DGEMM";  snames[1] = "DSYMM";  snames[2] = "DTRMM";
        snames[3] = "DTRSM";  snames[4] = "DSYRK";  snames[5] = "DSYR2K";
        snames[6] = "DGEMMTR";
    }
    #if T == C32 {
        snames[0] = "CGEMM";  snames[1] = "CHEMM";  snames[2] = "CSYMM";
        snames[3] = "CTRMM";  snames[4] = "CTRSM";  snames[5] = "CHERK";
        snames[6] = "CSYRK";  snames[7] = "CHER2K"; snames[8] = "CSYR2K";
        snames[9] = "CGEMMTR";
    }
    #if T == C64 {
        snames[0] = "ZGEMM";  snames[1] = "ZHEMM";  snames[2] = "ZSYMM";
        snames[3] = "ZTRMM";  snames[4] = "ZTRSM";  snames[5] = "ZHERK";
        snames[6] = "ZSYRK";  snames[7] = "ZHER2K"; snames[8] = "ZSYR2K";
        snames[9] = "ZGEMMTR";
    }

    {
        thresh = ncast(U, String.parse_float(*data[7][0]));
        println("thresh = %", thresh);

        nidim = String.parse_int(*data[8][0]);
        assert(nidim > 0 && nidim <= nidmax, "nidim out of range: 0 < % <= %\n", nidim, nidmax);
        println("nidim = %", nidim);

        idim = NewArray(nidim, int, false);
        for 0..nidim-1 {
            idim[it] = String.parse_int(*data[9][it]);
            assert(idim[it] >= 0 && idim[it] <= nmax, "idim[%] out of range: % <= % <= %\n", it, 0, idim[it], nmax);
        }
        println("idim = %", idim);

        nalf = String.parse_int(*data[10][0]);
        assert(nalf > 0 && nalf <= nalmax, "nalf out of range: 0 < % <= %\n", nalf, nalmax);
        println("nalf = %", nalf);

        for 0..nalf-1 {
            #if T == float32 {
                fset(*alf, it+1, String.parse_float(*data[11][it]));
            }
            #if T == float64 {
                fset(*alf, it+1, cast(float64) String.parse_float(*data[11][it]));
            }
            #if T == C32 {
                c_num := String.slice(data[11][it], 1, data[11][it].count-2);
                f_vals := String.split(c_num, ",");
                fset(*alf, it+1, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])));
            }
            #if T == C64 {
                c_num := String.slice(data[11][it], 1, data[11][it].count-2);
                f_vals := String.split(c_num, ",");
                fset(*alf, it+1, c64(cast(float64) String.parse_float(*f_vals[0]), cast(float64) String.parse_float(*f_vals[1])));
            }
        }
        println("alf = %", alf);

        nbet = String.parse_int(*data[12][0]);
        assert(nbet > 0 && nbet <= nbemax, "nbet out of range: 0 < % <= %\n", nbet, nbemax);
        println("nbet = %", nbet);

        for 0..nbet-1 {
            #if T == float32 {
                fset(*bet, it+1, String.parse_float(*data[13][it]));
            }
            #if T == float64 {
                fset(*bet, it+1, cast(float64) String.parse_float(*data[13][it]));
            }
            #if T == C32 {
                c_num := String.slice(data[13][it], 1, data[13][it].count-2);
                f_vals := String.split(c_num, ",");
                fset(*bet, it+1, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])));
            }
            #if T == C64 {
                c_num := String.slice(data[13][it], 1, data[13][it].count-2);
                f_vals := String.split(c_num, ",");
                fset(*bet, it+1, c64(cast(float64) String.parse_float(*f_vals[0]), cast(float64) String.parse_float(*f_vals[1])));
            }
        }
        println("bet = %", bet);

        for i:0..nsubs-1 {
            col := 1;
            do_test := data[14+i][col];
            while do_test == " " || do_test == "" {
                col += 1;
                do_test = data[14+i][col];
            }

            if do_test == "T" {
                for *ltest {
                    if snames[it_index] == data[14+i][0] {
                        it.* = true;
                    }
                }
            }
        }
        println("ltest = %", ltest);
    }

    eps := ifx U == float32 then ncast(U, EPS_32) else ncast(U, EPS_64);
    zero_t := zero(T);
    one_t  := one(T);
    zero_u := zero(U);
    trans_h : u8 = ---;
    if complex_suite {
        trans_h = #char "c";
    } else {
        trans_h = #char "t";
    }

    {
        ab_exact : Array2D(nmax, nmax, T);
        bb_exact : Array2D(nmax, nmax, T);
        c_exact  : Array2D(nmax, nmax, T);
        cc_exact : Array1D(nmax, T);

        n := min(32, nmax);
        for j:1..n {
            for i:1..n {
                fset(*ab_exact, i, j, ncast(T, max(i-j+1, 0)));
            }
            fset(*bb_exact, j, 1, ncast(T, j));
            fset(*bb_exact, 1, j, ncast(T, j));
            fset(*c_exact, j, 1, zero_t);
        }
        for j:1..n {
            fset(*cc_exact, j, ncast(T, j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3));
        }

        cc_exact_view := array2d_view(*cc_exact, nmax, row_major=false);

        errv : U = zero_u;
        pass := cblat3_cmmch(#char "n", #char "n", n, 1, n, one_t, ab_exact, nmax, bb_exact, nmax, zero_t, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != zero_u then err("cmmch result not good: transa = n, transb = n", exit_program=true);

        errv = zero_u;
        pass = cblat3_cmmch(#char "n", trans_h, n, 1, n, one_t, ab_exact, nmax, bb_exact, nmax, zero_t, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != zero_u then err("cmmch result not good: transa = n, transb = %", char_to_string(trans_h), exit_program=true);

        for j:1..n {
            fset(*bb_exact, j, 1, ncast(T, n-j+1));
            fset(*bb_exact, 1, j, ncast(T, n-j+1));
        }
        for j:1..n {
            fset(*cc_exact, n-j+1, ncast(T, j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3));
        }

        errv = zero_u;
        pass = cblat3_cmmch(trans_h, #char "n", n, 1, n, one_t, ab_exact, nmax, bb_exact, nmax, zero_t, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != zero_u then err("cmmch result not good: transa = %, transb = n", char_to_string(trans_h), exit_program=true);

        errv = zero_u;
        pass = cblat3_cmmch(trans_h, trans_h, n, 1, n, one_t, ab_exact, nmax, bb_exact, nmax, zero_t, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != zero_u then err("cmmch result not good: transa = %, transb = %", char_to_string(trans_h), char_to_string(trans_h), exit_program=true);
    }

    pass := true;
    for isnum:0..nsubs-1 {
        if !ltest[isnum] then continue;

        srnamt := snames[isnum];
        println("%: %", isnum, srnamt);

        if String.slice(srnamt, 1, 5) == "GEMM" {
            pass &= cblat3_cchk1(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g
            );
        } else if String.slice(srnamt, 1, 5) == "HEMM" || String.slice(srnamt, 1, 5) == "SYMM" {
            pass &= cblat3_cchk2(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g
            );
        } else if String.slice(srnamt, 1, 5) == "TRMM" || String.slice(srnamt, 1, 5) == "TRSM" {
            pass &= cblat3_cchk3(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *ct, *g, *c
            );
        } else if String.slice(srnamt, 1, 4) == "HERK" || String.slice(srnamt, 1, 4) == "SYRK" {
            pass &= cblat3_cchk4(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g
            );
        } else if String.slice(srnamt, 1, 5) == "HER2K" || String.slice(srnamt, 1, 5) == "SYR2K" {
            pass &= cblat3_cchk5(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g, *w
            );
        } else if String.slice(srnamt, 1, 6) == "GEMMTR" {
            pass &= cblat3_cchk6(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g
            );
        } else {
            err("_blat3 | unexpected routine name %", srnamt, exit_program=true);
        }
    }

    if !pass then err("% did not pass.", label);
    return pass;
}

cblat3 :: () -> bool {
    return _blat3(C32, float32, "LAPACK/Testing/cblat3.in", "cblat3");
}

sblat3 :: () -> bool {
    return _blat3(float32, float32, "LAPACK/Testing/sblat3.in", "sblat3");
}

dblat3 :: () -> bool {
    return _blat3(float64, float64, "LAPACK/Testing/dblat3.in", "dblat3");
}

zblat3 :: () -> bool {
    return _blat3(C64, float64, "LAPACK/Testing/zblat3.in", "zblat3");
}
