// ? LAPACK BLAT2 generic driver

#load "_blat2_checks.jai";

_blat2 :: ($T: Type, $U: Type, file_url: string, label: string) -> bool {
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
    incmax : int : 2;

    ninmax : int : 7;
    nidmax : int : 9;
    nkbmax : int : 7;
    nalmax : int : 7;
    nbemax : int : 7;

    complex_suite :: #run is_non_real_number(T);
    #if T == float32 || T == float64 {
        nsubs :: 16;
    } else {
        nsubs :: 17;
    }

    a : Array2D(nmax, nmax, T);
    aa, as : Array1D(nmax*nmax, T);
    x, y, yt : Array1D(nmax, T);
    z : Array1D(2*nmax, T);
    xs, _xx, ys, yy : Array1D(nmax*incmax, T);
    alf : Array1D(nalmax, T);
    bet : Array1D(nbemax, T);
    g   : Array1D(nmax, U);

    thresh : U;
    nidim, nkb, ninc, nalf, nbet : int;
    idim, kb, inc : []int;
    ltest : [nsubs]bool;
    snames := NewArray(nsubs, string, false);
    #if T == float32 {
        snames[0]  = "SGEMV"; snames[1]  = "SGBMV"; snames[2]  = "SSYMV"; snames[3]  = "SSBMV";
        snames[4]  = "SSPMV"; snames[5]  = "STRMV"; snames[6]  = "STBMV"; snames[7]  = "STPMV";
        snames[8]  = "STRSV"; snames[9]  = "STBSV"; snames[10] = "STPSV"; snames[11] = "SGER";
        snames[12] = "SSYR";  snames[13] = "SSPR";  snames[14] = "SSYR2"; snames[15] = "SSPR2";
    }
    #if T == float64 {
        snames[0]  = "DGEMV"; snames[1]  = "DGBMV"; snames[2]  = "DSYMV"; snames[3]  = "DSBMV";
        snames[4]  = "DSPMV"; snames[5]  = "DTRMV"; snames[6]  = "DTBMV"; snames[7]  = "DTPMV";
        snames[8]  = "DTRSV"; snames[9]  = "DTBSV"; snames[10] = "DTPSV"; snames[11] = "DGER";
        snames[12] = "DSYR";  snames[13] = "DSPR";  snames[14] = "DSYR2"; snames[15] = "DSPR2";
    }
    #if T == C32 {
        snames[0]  = "CGEMV"; snames[1]  = "CGBMV"; snames[2]  = "CHEMV"; snames[3]  = "CHBMV";
        snames[4]  = "CHPMV"; snames[5]  = "CTRMV"; snames[6]  = "CTBMV"; snames[7]  = "CTPMV";
        snames[8]  = "CTRSV"; snames[9]  = "CTBSV"; snames[10] = "CTPSV"; snames[11] = "CGERC";
        snames[12] = "CGERU"; snames[13] = "CHER";  snames[14] = "CHPR";  snames[15] = "CHER2";
        snames[16] = "CHPR2";
    }
    #if T == C64 {
        snames[0]  = "ZGEMV"; snames[1]  = "ZGBMV"; snames[2]  = "ZHEMV"; snames[3]  = "ZHBMV";
        snames[4]  = "ZHPMV"; snames[5]  = "ZTRMV"; snames[6]  = "ZTBMV"; snames[7]  = "ZTPMV";
        snames[8]  = "ZTRSV"; snames[9]  = "ZTBSV"; snames[10] = "ZTPSV"; snames[11] = "ZGERC";
        snames[12] = "ZGERU"; snames[13] = "ZHER";  snames[14] = "ZHPR";  snames[15] = "ZHER2";
        snames[16] = "ZHPR2";
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

        nkb = String.parse_int(*data[10][0]);
        assert(nkb > 0 && nkb <= nkbmax, "nkb out of range: 0 < % <= %\n", nkb, nkbmax);
        println("nkb = %", nkb);

        kb = NewArray(nkb, int, false);
        for 0..nkb-1 {
            kb[it] = String.parse_int(*data[11][it]);
            assert(kb[it] >= 0, "kb[%] out of range: % < 0\n", it, kb[it]);
        }
        println("kb = %", kb);

        ninc = String.parse_int(*data[12][0]);
        assert(ninc >= 1 && ninc <= ninmax, "ninc out of range: 0 < % <= %\n", ninc, ninmax);
        println("ninc = %", ninc);

        inc = NewArray(ninc, int, false);
        for 0..ninc-1 {
            inc[it] = String.parse_int(*data[13][it]);
            assert(inc[it] != 0 && abs(inc[it]) <= incmax, "inc[%] out of range: 0 != %, % <= %", it, inc[it], inc[it], incmax);
        }
        println("inc = %", inc);

        nalf = String.parse_int(*data[14][0]);
        assert(nalf > 0 && nalf <= nalmax, "nalf out of range: 0 < % <= %\n", nalf, nalmax);
        println("nalf = %", nalf);

        for 0..nalf-1 {
            #if T == float32 {
                set(*alf, it, String.parse_float(*data[15][it]));
            }
            #if T == float64 {
                set(*alf, it, cast(float64) String.parse_float(*data[15][it]));
            }
            #if T == C32 {
                c_num := String.slice(data[15][it], 1, data[15][it].count-2);
                f_vals := String.split(c_num, ",");
                set(*alf, it, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])));
            }
            #if T == C64 {
                c_num := String.slice(data[15][it], 1, data[15][it].count-2);
                f_vals := String.split(c_num, ",");
                set(*alf, it, c64(cast(float64) String.parse_float(*f_vals[0]), cast(float64) String.parse_float(*f_vals[1])));
            }
        }
        println("alf = %", alf);

        nbet = String.parse_int(*data[16][0]);
        assert(nbet > 0 && nbet <= nbemax, "nbet out of range: 0 < % <= %\n", nbet, nbemax);
        println("nbet = %", nbet);

        for 0..nbet-1 {
            #if T == float32 {
                set(*bet, it, String.parse_float(*data[17][it]));
            }
            #if T == float64 {
                set(*bet, it, cast(float64) String.parse_float(*data[17][it]));
            }
            #if T == C32 {
                c_num := String.slice(data[17][it], 1, data[17][it].count-2);
                f_vals := String.split(c_num, ",");
                set(*bet, it, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])));
            }
            #if T == C64 {
                c_num := String.slice(data[17][it], 1, data[17][it].count-2);
                f_vals := String.split(c_num, ",");
                set(*bet, it, c64(cast(float64) String.parse_float(*f_vals[0]), cast(float64) String.parse_float(*f_vals[1])));
            }
        }
        println("bet = %", bet);

        for i:0..nsubs-1 {
            col := 1;
            do_test := data[18+i][col];
            while do_test == " " || do_test == "" {
                col += 1;
                do_test = data[18+i][col];
            }

            if do_test == "T" {
                for *ltest {
                    if snames[it_index] == data[18+i][0] {
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

    n := min(32, nmax);
    for j:1..n {
        for i:1..n {
            fset(*a, i, j, ncast(T, max(i-j+1, 0)));
        }
        fset(*x, j, ncast(T, j-1));
        fset(*y, j, zero_t);
    }
    for j:1..n {
        fset(*yy, j, ncast(T, (j-1)*(j*(j-1))/2 - (j*(j-1)*(j-2))/3));
    }

    pass := true;
    trans : u8 = #char "n";
    pass &= cmvch(trans, n, n, one_t, a, nmax, x, 1, zero_t, y, 1, *yt, *g, yy, eps, true);
    pass &= lce(yy, yt, n);
    if !pass then err("cmvch result not good: n = %\n%\n%", n, yy, yt, exit_program=true);

    trans = #char "t";
    pass &= cmvch(trans, n, n, one_t, a, nmax, x, -1, zero_t, y, -1, *yt, *g, yy, eps, true);
    pass &= lce(yy, yt, n);
    if !pass then err("cmvch result not good: n = %\n%\n%", n, yy, yt, exit_program=true);

    for isnum:0..nsubs-1 {
        if !ltest[isnum] then continue;

        srnamt := snames[isnum];
        println("%: %", isnum, srnamt);

        if isnum <= 1 {
            pass &= cchk1(
                srnamt, eps, thresh, true,
                nidim, idim, nkb, kb, nalf, alf.data,
                nbet, bet.data, ninc, inc, nmax, incmax, *a, *aa, *as,
                *x, *_xx, *xs, *y, *yy, *ys, *yt, *g
            );
        } else if isnum <= 4 {
            pass &= cchk2(
                srnamt, eps, thresh, true,
                nidim, idim, nkb, kb, nalf, alf.data,
                nbet, bet.data, ninc, inc, nmax, incmax, *a, *aa, *as,
                *x, *_xx, *xs, *y, *yy, *ys, *yt, *g
            );
        } else if isnum <= 10 {
            pass &= cchk3(
                srnamt, eps, thresh, true,
                nidim, idim, nkb, kb, ninc, inc,
                nmax, incmax, *a, *aa, *as, *y, *yy, *ys, *yt, *g, *z
            );
        } else if complex_suite {
            if isnum <= 12 {
                pass &= cchk4(
                    srnamt, eps, thresh, true,
                    nidim, idim, nalf, alf.data, ninc, inc,
                    nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
                    *yt, *g, *z
                );
            } else if isnum <= 14 {
                pass &= cchk5(
                    srnamt, eps, thresh, true,
                    nidim, idim, nalf, alf.data, ninc, inc,
                    nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
                    *yt, *g, *z
                );
            } else {
                pass &= cchk6(
                    srnamt, eps, thresh, true,
                    nidim, idim, nalf, alf.data, ninc, inc,
                    nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
                    *yt, *g, *z
                );
            }
        } else {
            if isnum == 11 {
                pass &= cchk4(
                    srnamt, eps, thresh, true,
                    nidim, idim, nalf, alf.data, ninc, inc,
                    nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
                    *yt, *g, *z
                );
            } else if isnum <= 13 {
                pass &= cchk5(
                    srnamt, eps, thresh, true,
                    nidim, idim, nalf, alf.data, ninc, inc,
                    nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
                    *yt, *g, *z
                );
            } else {
                pass &= cchk6(
                    srnamt, eps, thresh, true,
                    nidim, idim, nalf, alf.data, ninc, inc,
                    nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
                    *yt, *g, *z
                );
            }
        }
    }

    if !pass then err("% did not pass.", label);
    return pass;
}

cblat2 :: () -> bool {
    return _blat2(C32, float32, "LAPACK/Testing/cblat2.in", "cblat2");
}

sblat2 :: () -> bool {
    return _blat2(float32, float32, "LAPACK/Testing/sblat2.in", "sblat2");
}

dblat2 :: () -> bool {
    return _blat2(float64, float64, "LAPACK/Testing/dblat2.in", "dblat2");
}

zblat2 :: () -> bool {
    return _blat2(C64, float64, "LAPACK/Testing/zblat2.in", "zblat2");
}
