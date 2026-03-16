// ? LAPACK cblat2
// ? http://www.netlib.org/lapack/explore-html/d2/d54/cblat2_8f_source.html


// ? Does is make sense to accurately copy the behavior of the test suite?
// ? Meaning, generating output files etc.?


file_url : string : "LAPACK/Testing/cblat2.in";

//  ######  ########  ##          ###    ########  #######
// ##    ## ##     ## ##         ## ##      ##    ##     ##
// ##       ##     ## ##        ##   ##     ##           ##
// ##       ########  ##       ##     ##    ##     #######
// ##       ##     ## ##       #########    ##    ##
// ##    ## ##     ## ##       ##     ##    ##    ##
//  ######  ########  ######## ##     ##    ##    #########

cblat2 :: () -> bool {
    println_push("cblat2", color=.FG_GREEN);

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


    // constants
    nin    : int : 5;
    nsubs  : int : 17;

    nmax   : int : 65;
    incmax : int : 2;

    ninmax : int : 7;
    nidmax : int : 9;
    nkbmax : int : 7;
    nalmax : int : 7;
    nbemax : int : 7;

    czero := zero(C32);
    cone  := one(C32);
    rzero := zero(float32);

    // variables
    a : Array2D(nmax, nmax, C32); // ! There might be an issue with row vs column major interpretation of `a`
    aa, as : Array1D(nmax*nmax, C32);
    x, y, yt : Array1D(nmax, C32);
    z : Array1D(2*nmax, C32);
    xs, _xx, ys, yy : Array1D(nmax*incmax, C32);
    alf : Array1D(nalmax, C32);
    bet : Array1D(nbemax, C32);
    g   : Array1D(nmax, float32);

    thresh : float32;
    nidim, nkb, ninc, nalf, nbet : int;
    idim, kb, inc : []int;
    snames : [nsubs]string : string.[
        "CGEMV", "CGBMV", "CHEMV", "CHBMV",
        "CHPMB", "CTRMV", "CTBMV", "CTPMV",
        "CTRSV", "CTBSV", "CTPSV", "CGERC",
        "CGERU", "CHER" , "CHPR" , "CHER2",
        "CHPR2"
    ];
    ltest : [nsubs]bool; // default to false
    // tsterr : bool = ---; // not used here


    // parse data
    {
        // tsterr = ifx data[6][0] == "T" then true else false; // not used here

        thresh = String.parse_float(*data[7][0]);
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
        assert(nkb > 0 && nkb <= nkbmax, "nkb out of range: % < % <= %\n", 0, nkb, nkbmax);
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
            // complex number!
            c_num := String.slice(data[15][it], 1, data[15][it].count-2); // remove "(" and ")"
            f_vals := String.split(c_num, ",");
            set(*alf, it, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])) );
        }
        println("alf = %", alf);

        nbet = String.parse_int(*data[16][0]);
        assert(nbet > 0 && nbet <= nalmax, "nbet out of range: 0 < % <= %\n", nbet, nalmax);
        println("nbet = %", nbet);

        for 0..nbet-1 {
            // complex number!
            c_num := String.slice(data[17][it], 1, data[17][it].count-2); // remove "(" and ")"
            f_vals := String.split(c_num, ",");
            set(*bet, it, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])) );
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


    eps := EPS_32; // Complex(float32)

    // Check the riliability of CMVCH using exact data.
    n := min(32, nmax);
    for j:1..n {
        for i:1..n {
            fset(*a, i, j, ncast(C32,max(i-j+1, 0)));
        }
        fset(*x, j, ncast(C32,j-1));
        fset(*y, j, czero);
    }
    for j:1..n {
        fset(*yy, j, ncast(C32, (j-1)*(j*(j-1))/2 - (j*(j-1)*(j-2))/3) );
    }
    // YY holds the exact result. On exit from CMVCH YT holds 
    // the result computed by CMVCH.

    pass := true;

    trans : u8 = #char "n";
    pass &= cmvch(trans, n, n, cone, a, nmax, x,  1, czero, y,  1, *yt, *g, yy, eps, true);
    pass &= lce(yy, yt, n);
    if !pass then err("cmvch result not good: n = %\n%\n%", n, yy, yt, exit_program=true);

    trans = #char "t";
    pass &= cmvch(trans, n, n, cone, a, nmax, x, -1, czero, y, -1, *yt, *g, yy, eps, true);
    pass &= lce(yy, yt, n);
    if !pass then err("cmvch result not good: n = %\n%\n%", n, yy, yt, exit_program=true);


    // Test each subroutine in turn.
    // println("Test each subroutine in turn.");
    for isnum:0..nsubs-1 {
        if !ltest[isnum] then continue;

        srnamt := snames[isnum];
        println("%: %", isnum, srnamt);
        // error exits need not be tested

        // Test computations
        if isnum <= 1 {
            // pass &= cchk1(
            //     srnamt, eps, thresh, true,
            //     nidim, idim, nkb, kb, nalf, alf.data,
            //     nbet, bet.data, ninc, inc, nmax, incmax, *a, *aa, *as,
            //     *x, *_xx, *xs, *y, *yy, *ys, *yt, *g
            // );
        } else if isnum <= 4 {
            // pass &= cchk2(
            //     srnamt, eps, thresh, true,
            //     nidim, idim, nkb, kb, nalf, alf.data,
            //     nbet, bet.data, ninc, inc, nmax, incmax, *a, *aa, *as,
            //     *x, *_xx, *xs, *y, *yy, *ys, *yt, *g
            // );
        } else if isnum <= 10 {
            // pass &= cchk3(
            //     srnamt, eps, thresh, true,
            //     nidim, idim, nkb, kb, ninc, inc,
            //     nmax, incmax, *a, *aa, *as, *y, *yy, *ys, *yt, *g, *z
            // );
        } else if isnum <= 12 {
            // pass &= cchk4(
            //     srnamt, eps, thresh, true,
            //     nidim, idim, nalf, alf.data, ninc, inc,
            //     nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
            //     *yt, *g, *z
            // );
        } else if isnum <= 14 {
            pass &= cchk5(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, ninc, inc,
                nmax, incmax, *a, *aa, *as, *x, *_xx, *xs, *y, *yy, *ys,
                *yt, *g, *z
            );
        } else {
            
        }
    }

    // Clean up.
    if !pass then err("cblat2 did not pass.");
    return pass;
}


//  ######  ##     ## ##     ##  ######  ##     ##
// ##    ## ###   ### ##     ## ##    ## ##     ##
// ##       #### #### ##     ## ##       ##     ##
// ##       ## ### ## ##     ## ##       #########
// ##       ##     ##  ##   ##  ##       ##     ##
// ##    ## ##     ##   ## ##   ##    ## ##     ##
//  ######  ##     ##    ###     ######  ##     ##

cmvch :: (
    trans: u8,
    m: int,
    n: int,
    alpha: $T,
    a: $A/ArrayType2D,
    nmax: int,
    x: $X/ArrayType1D,
    incx: int,
    beta: T,
    y: $Y/ArrayType1D,
    $$incy: int,
    yt: *$YT/ArrayType1D,
    g: *$G/ArrayType1D,
    yy: $YY/ArrayType1D,
    eps: $U,
    $fatal: bool
) -> bool {
    #if CHECKS {
        #run {
            assert(!is_non_real_number(U));
            assert(U == get_number_inner_type(T));
            assert(T ==  A.data_type);
            assert(T ==  X.data_type);
            assert(T ==  Y.data_type);
            assert(T == YY.data_type);
            assert(U ==  G.data_type);
        }
    }
    DEBUG :: false;
    #if DEBUG println("trans = %", char_to_string(trans), color=.FG_CYAN);

    println_push("cmvch");

    czero := zero(T);
    rzero := zero(U);
    rone  := one(U);

    tran  := char_eq(trans, #char "t");
    ctran := char_eq(trans, #char "c");
    #if DEBUG println("tran = %, ctran = %", tran, ctran, color=.FG_CYAN);
    ml := ifx tran || ctran then n else m;
    nl := ifx tran || ctran then m else n;
    kx    := ifx incx < 0 then nl else 1;
    incxl := ifx incx < 0 then -1 else 1;
    ky    := ifx incy < 0 then ml else 1;
    incyl := ifx incy < 0 then -1 else 1;
    iy, jx: int = ---;

    // Compute expected result in YT using data in A, X and Y.
    // Compute gauges in G.
    iy = ky;
    for i:1..ml {
        fset(yt, iy, czero);
        fset(g , iy, rzero);
        jx = kx;
        if tran {
            for j:1..nl {
                fset(yt, iy, fget(yt.*,iy) + fget(a,j,i)*fget(x,jx));
                fset(g , iy, fget(g.*,iy) + abs1(fget(a,j,i))*abs1(fget(x,jx)));
                jx += incxl;
            }
        } else if ctran {
            for j:1..nl {
                fset(yt, iy, fget(yt.*,iy) + conjugate(fget(a,j,i))*fget(x,jx));
                fset( g, iy, fget(g.*,iy) + abs1(fget(a,j,i))*abs1(fget(x,jx)));
                jx += incxl;
            }
        } else {
            for j:1..nl {
                fset(yt, iy, fget(yt.*,iy) + fget(a,i,j)*fget(x,jx));
                fset( g, iy, fget(g.*,iy) + abs1(fget(a,i,j))*abs1(fget(x,jx)));
                jx += incxl;
            }
        }
        fset(yt, iy, alpha*fget(yt.*,iy) + beta*fget(y,iy));
        fset( g, iy, abs1(alpha)*fget(g.*,iy) + abs1(beta)*abs1(fget(y,iy)));
        iy += incyl;
    }

    // Compute the error ration for this result.
    rerr := rzero;
    res := true;
    for i:1..ml {
        yii := 1 + (i-1)*abs(incy);
        erri := abs(fget(yt.*,i) - fget(yy,yii))/eps;
        if fget(g.*,i) != rzero {
            erri /= fget(g.*,i);
        }
        rerr = max(rerr,erri);
        if rerr*sqrt(eps) >= rone {
            _trans := char_to_string(trans);
            defer free(_trans);
            err("cmvch | trans = %", _trans);
            err("cmvch | m     = %", m    );
            err("cmvch | n     = %", n    );
            err("cmvch | alpha = %", alpha);
            err("cmvch | a     = %", a    );
            err("cmvch | nmax  = %", nmax );
            err("cmvch | x     = %", x    );
            err("cmvch | incx  = %", incx );
            err("cmvch | beta  = %", beta );
            err("cmvch | y     = %", y    );
            err("cmvch | incy  = %", incy );
            err("cmvch | yt    = %", yt.*   );
            err("cmvch | g     = %", g.*    );
            err("cmvch | yy    = %", yy   );
            err("cmvch | eps   = %", formatFloat(eps,mode=FormatFloat.Mode.SCIENTIFIC)  );
            err("cmvch | erri % : %", type_of(erri), erri);
            err("cmvch | yt[%] = %, yy[1 + (%-1)*abs(%)] = %", i, fget(yt.*,i), i, incy, fget(yy,yii));
            err("cmvch | abs(get(yt,%) - get(yy,1 + (%-1)*abs(%))) = %", i, i, incy, abs(fget(yt.*,i) - fget(yy,yii)));
            err("cmvch | g[%] = %", i, fget(g.*,i));
            err("cmvch | %*sqrt(%) = % >= 1.0", rerr, formatFloat(eps, mode=FormatFloat.Mode.SCIENTIFIC), rerr*sqrt(eps));
            err("", exit_program=fatal);
            res = false;
        }
    }
    return res;
}


// ##        ######  ########
// ##       ##    ## ##
// ##       ##       ##
// ##       ##       ######
// ##       ##       ##
// ##       ##    ## ##
// ########  ######  ########

lce :: (ri: $A/ArrayType1D, rj: $B/ArrayType1D, $$lr: int) -> bool {
    for 0..lr-1 {
        if get(ri,it) != get(rj,it) then return false;
    }
    return true;
}

// ##        ######  ######## ########  ########  ######
// ##       ##    ## ##       ##     ## ##       ##    ##
// ##       ##       ##       ##     ## ##       ##
// ##       ##       ######   ########  ######    ######
// ##       ##       ##       ##   ##   ##             ##
// ##       ##    ## ##       ##    ##  ##       ##    ##
// ########  ######  ######## ##     ## ########  ######

// Tests if selected elements in two arrays are equal.
//
// TYPE is 'GE', 'HE' or 'HP'.
//
// Auxiliary routine for test program for Level 2 Blas.
//
// -- Written on 10-August-1987.
//    Richard Hanson, Sandia National Labs.
//    Jeremy Du Croz, NAG Central Office.

lceres :: (
    type: string, 
    uplo: u8, 
    m: int, 
    n: int, 
    aa: $AA/ArrayType2D, 
    as: $AS/ArrayType2D, 
    lda: int
) -> bool {
    upper := char_eq(uplo, #char "u");
    if String.equal_nocase(type, "ge") {
        for j:0..n-1 {
            for i:m+1..lda-1 {
                if get(aa, i, j) != get(as, i, j) {
                    err("aa[%, %] = % != as[%, %] = %", i, j, get(aa,i,j), i, j, get(as,i,j));
                    return false;
                }
            }
        }
    } else if String.equal_nocase(type, "he") {
        for j:0..n-1 {
            ibeg := ifx upper then 0 else j;
            iend := ifx upper then j else n-1;

            if ibeg > 0 {
                for i:0..ibeg-1 {
                    if get(aa, i, j) != get(as, i, j) {
                        err("aa[%, %] = % != as[%, %] = %", i, j, get(aa,i,j), i, j, get(as,i,j));
                        return false;
                    }
                }
            }

            if iend+1 <= lda-1 {
                for i:iend+1..lda-1 {
                    if get(aa, i, j) != get(as, i ,j) {
                        err("aa[%, %] = % != as[%, %] = %", i, j, get(aa,i,j), i, j, get(as,i,j));
                        return false;
                    }
                }
            }
        }
    }

    return true;
}


//  ######   ######  ##     ## ##    ##    ##
// ##    ## ##    ## ##     ## ##   ##   ####
// ##       ##       ##     ## ##  ##      ##
// ##       ##       ######### #####       ##
// ##       ##       ##     ## ##  ##      ##
// ##    ## ##    ## ##     ## ##   ##     ##
//  ######   ######  ##     ## ##    ##  ######

cchk1 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nkb: int,
    kb: []int,
    nalf: int,
    alf: []$T,
    nbet: int,
    bet: []T,
    ninc: int,
    inc: []int,
    nmax: int,
    incmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    x: *$X/ArrayType1D,
    _xx: *$XX/ArrayType1D,
    xs: *$XS/ArrayType1D,
    y: *$Y/ArrayType1D,
    yy: *$YY/ArrayType1D,
    ys: *$YS/ArrayType1D,
    yt: *$YT/ArrayType1D,
    g: *$G/ArrayType1D
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == X.data_type);
        assert(T == XX.data_type);
        assert(T == XS.data_type);
        assert(T == Y.data_type);
        assert(T == YY.data_type);
        assert(T == YS.data_type);
        assert(T == YT.data_type);
        assert(U == G.data_type);
    }
    DEBUG :: true;
    println_push("cchk1");

    pass := true;

    ich :: "ntc";

    czero := zero(T);
    cone  := one(T);
    chalf := ncast(T, 0.5);
    rzero := zero(U);
    rone  := one(U);

    full   := char_eq(sname[2], #char "e");
    banded := char_eq(sname[2], #char "b");

    // Define the number of arguments
    nargs := ifx full then 11 else 13;

    nc := 0;
    reset := true;
    errmax := rzero;

    #if DEBUG {
        COUNT := 0;
    }

    for in:0..nidim-1 {
        // println_push("in = %", in);
        n := idim[in];
        nd := n/2 + 1;

        for im:0..1 {
            // println_push("im = %", im);
            m  := ifx im == 0 then max(n-nd,0) else min(n+nd, nmax);
            nk := ifx banded  then nkb         else 1;

            for iku:0..nk-1 {
                // println_push("iku = %", iku);
                ku := ifx banded then kb[iku]     else n-1;
                kl := ifx banded then max(ku-1,0) else m-1;

                // Set lda to 1 more than minimum value if room
                lda := ifx banded then kl+ku+1 else m;
                if lda < nmax then lda += 1;

                // Skip tests if not enough room.
                if lda > nmax then continue;

                laa := lda*n;
                _null := n <= 0 || m <= 0;

                // Generate the matrix a.
                transl := czero;
                _slice := String.slice(sname, 1, 2);
                cmake(_slice, 0, 0, m, n, a, nmax, aa, lda, kl, ku, *reset, transl);

                for ic:0..3-1 {
                    trans := ich[ic];
                    tran := char_eq(trans, #char "t", #char "c");

                    ml := ifx tran then n else m;
                    nl := ifx tran then m else n;

                    for ix:0..ninc-1 {
                        incx := inc[ix];
                        lx := abs(incx)*nl;

                        // Generate the vector X.
                        transl = chalf;
                        x_view := array2d_view(x, 1, row_major=false);
                        cmake("ge", 0, 0, 1, nl, *x_view, 1, _xx, abs(incx), 0, nl-1, *reset, transl);

                        if nl > 1 {
                            set(x, nl/2 - 1, czero);
                            set(_xx, abs(incx)*(nl/2-1), czero);
                        }

                        for iy:0..ninc-1 {
                            incy := inc[iy];
                            ly := abs(incy)*ml;

                            for ia:0..nalf-1 {
                                alpha := alf[ia];

                                for ib:0..nbet-1 {
                                    #if DEBUG {
                                        COUNT += 1;
                                        println("COUNT %", COUNT);
                                    }
                                    println("in %, im %, iku %, ic %, ix %, iy %, ia %, ib %", in, im, iku, ic, ix, iy, ia, ib);
                                    beta := bet[ib];

                                    // Generate the vector Y.

                                    transl := czero;
                                    y_view := array2d_view(y, 1, row_major=false);
                                    cmake("ge", 0, 0, 1, ml, *y_view, 1, yy, abs(incy), 0, ml-1, *reset, transl);

                                    nc += 1;

                                    // Save every datum before calling hte subroutine

                                    transs := trans;
                                    ms := m;
                                    ns := n;
                                    kls := kl;
                                    kus := ku;
                                    als := alpha;
                                    for i:0..laa-1 {
                                        set(as, i, get(aa,i));
                                    }
                                    ldas := lda;
                                    for i:0..lx-1 {
                                        set(xs, i, get(_xx,i));
                                    }
                                    incxs := incx;
                                    bls := beta;
                                    for i:0..ly-1 {
                                        set(ys, i, get(yy,i));
                                    }
                                    incys := incy;

                                    // Call the subroutine.
                                    aa_view := array2d_view(aa, lda, row_major=false);
                                    if full {
                                        println("cgemv");
                                        cgemv(trans, m, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                    } else if banded {
                                        println("cgbmv");
                                        cgbmv(trans, m, n, kl, ku, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                    }

                                    // See what data changed inside subroutines.
                                    // most of this is unnecessary since Jai has immutable variables in the first place. - Fehse, 2022-01-14
                                    if _null {
                                        ok := lce(ys, yy, ly);
                                        pass &= ok;
                                        if !ok {
                                            err(" lce | ly = %", ly);
                                            err(" lce | ys = %", ys);
                                            err(" lce | yy = %", yy, exit_program=fatal);
                                        }
                                    } else {
                                        ys_view := array2d_view(ys, abs(incy), row_major=false);
                                        yy_view := array2d_view(yy, abs(incy), row_major=false);
                                        ok := lceres("GE", 0, 1, ml, ys_view, yy_view, abs(incy));
                                        pass &= ok;
                                        if !ok {
                                            err(" lceres | ml = %, abs(incy) = %", ml, abs(incy));
                                            err(" lceres | ys = %", ys_view);
                                            err(" lceres | yy = %", yy_view, exit_program=fatal);
                                        }
                                    }

                                    if !_null {
                                        // Check result.
                                        pass &= cmvch(
                                            trans, m, n, alpha, a,
                                            nmax, x, incx, beta, y,
                                            incy, yt, g, yy, eps,
                                            fatal
                                        );
                                    } else {
                                        // Avoid repeating tests with m <= 0 or n <= 0
                                        continue im;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return pass;
}


//  ######   ######  ##     ## ##    ##  #######
// ##    ## ##    ## ##     ## ##   ##  ##     ##
// ##       ##       ##     ## ##  ##          ##
// ##       ##       ######### #####     #######
// ##       ##       ##     ## ##  ##   ##
// ##    ## ##    ## ##     ## ##   ##  ##
//  ######   ######  ##     ## ##    ## #########

cchk2 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nkb: int,
    kb: []int,
    nalf: int,
    alf: []$T,
    nbet: int,
    bet: []T,
    ninc: int,
    inc: []int,
    nmax: int,
    incmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    x: *$X/ArrayType1D,
    _xx: *$XX/ArrayType1D,
    xs: *$XS/ArrayType1D,
    y: *$Y/ArrayType1D,
    yy: *$YY/ArrayType1D,
    ys: *$YS/ArrayType1D,
    yt: *$YT/ArrayType1D,
    g: *$G/ArrayType1D
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == X.data_type);
        assert(T == XX.data_type);
        assert(T == XS.data_type);
        assert(T == Y.data_type);
        assert(T == YY.data_type);
        assert(T == YS.data_type);
        assert(T == YT.data_type);
        assert(U == G.data_type);
    }
    DEBUG :: true;
    println_push("cchk2");

    pass := true;

    ichu :: "ul";

    czero := zero(T);
    chalf := ncast(T, 0.5);

    banded := char_eq(sname[2], #char "b");
    packed := char_eq(sname[2], #char "p");
    reset := true;

    #if DEBUG {
        COUNT := 0;
    }

    for in:0..nidim-1 {
        n := idim[in];
        nk := ifx banded then nkb else 1;

        for iku:0..nk-1 {
            k := ifx banded then kb[iku] else n-1;

            lda := ifx banded then k+1 else n;
            if !packed && lda < nmax then lda += 1;
            if !packed && lda > nmax then continue;

            laa := ifx packed then (n*(n+1))/2 else lda*n;
            _null := n <= 0;

            for iu:0..1 {
                uplo := ichu[iu];

                transl := czero;
                _slice := String.slice(sname, 1, 2);
                cmake(_slice, uplo, 0, n, n, a, nmax, aa, lda, k, k, *reset, transl);

                for ix:0..ninc-1 {
                    incx := inc[ix];
                    lx := abs(incx)*n;

                    transl = chalf;
                    x_view := array2d_view(x, 1, row_major=false);
                    cmake("ge", 0, 0, 1, n, *x_view, 1, _xx, abs(incx), 0, n-1, *reset, transl);

                    if n > 1 {
                        set(x, n/2 - 1, czero);
                        set(_xx, abs(incx)*(n/2-1), czero);
                    }

                    for iy:0..ninc-1 {
                        incy := inc[iy];
                        ly := abs(incy)*n;

                        for ia:0..nalf-1 {
                            alpha := alf[ia];

                            for ib:0..nbet-1 {
                                #if DEBUG {
                                    COUNT += 1;
                                    println("COUNT %", COUNT);
                                }
                                println("in %, iku %, iu %, ix %, iy %, ia %, ib %", in, iku, iu, ix, iy, ia, ib);
                                beta := bet[ib];

                                transl := czero;
                                y_view := array2d_view(y, 1, row_major=false);
                                cmake("ge", 0, 0, 1, n, *y_view, 1, yy, abs(incy), 0, n-1, *reset, transl);

                                for i:0..laa-1 {
                                    set(as, i, get(aa,i));
                                }
                                for i:0..lx-1 {
                                    set(xs, i, get(_xx,i));
                                }
                                for i:0..ly-1 {
                                    set(ys, i, get(yy,i));
                                }

                                if packed {
                                    println("chpmv");
                                    chpmv(uplo, n, alpha, aa, _xx, incx, beta, yy, incy);
                                } else {
                                    aa_view := array2d_view(aa, lda, row_major=false);
                                    if banded {
                                        println("chbmv");
                                        chbmv(uplo, n, k, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                    } else {
                                        println("chemv");
                                        chemv(uplo, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                    }
                                }

                                if _null {
                                    ok := lce(ys, yy, ly);
                                    pass &= ok;
                                    if !ok {
                                        err(" lce | ly = %", ly);
                                        err(" lce | ys = %", ys);
                                        err(" lce | yy = %", yy, exit_program=fatal);
                                    }
                                } else {
                                    ys_view := array2d_view(ys, abs(incy), row_major=false);
                                    yy_view := array2d_view(yy, abs(incy), row_major=false);
                                    ok := lceres("GE", 0, 1, n, ys_view, yy_view, abs(incy));
                                    pass &= ok;
                                    if !ok {
                                        err(" lceres | n = %, abs(incy) = %", n, abs(incy));
                                        err(" lceres | ys = %", ys_view);
                                        err(" lceres | yy = %", yy_view, exit_program=fatal);
                                    }
                                }

                                if !_null {
                                    pass &= cmvch(
                                        #char "n", n, n, alpha, a,
                                        nmax, x, incx, beta, y,
                                        incy, yt, g, yy, eps,
                                        fatal
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return pass;
}


//  ######   ######  ##     ## ##    ##  #######
// ##    ## ##    ## ##     ## ##   ##  ##     ##
// ##       ##       ##     ## ##  ##          ##
// ##       ##       ######### #####     #######
// ##       ##       ##     ## ##  ##          ##
// ##    ## ##    ## ##     ## ##   ##  ##     ##
//  ######   ######  ##     ## ##    ##  #######

cchk3 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nkb: int,
    kb: []int,
    ninc: int,
    inc: []int,
    nmax: int,
    incmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    x: *$X/ArrayType1D,
    _xx: *$XX/ArrayType1D,
    xs: *$XS/ArrayType1D,
    xt: *$XT/ArrayType1D,
    g: *$G/ArrayType1D,
    z: *$Z/ArrayType1D,
) -> bool {
    T :: A.data_type;
    #run {
        assert(U == get_number_inner_type(T));
        // assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == X.data_type);
        assert(T == XX.data_type);
        assert(T == XS.data_type);
        assert(U == G.data_type);
        assert(T == Z.data_type);
    }
    DEBUG :: true;
    println_push("cchk3");

    pass := true;

    ichu :: "ul";
    icht :: "ntc";
    ichd :: "un";

    czero :: #run zero(T);
    cone  :: #run one(T);
    chalf :: #run ncast(T, 0.5);

    full   := char_eq(sname[2], #char "r");
    banded := char_eq(sname[2], #char "b");
    packed := char_eq(sname[2], #char "p");

    name_mv := char_eq(sname[3], #char "m") && char_eq(sname[4], #char "v");
    name_sv := char_eq(sname[3], #char "s") && char_eq(sname[4], #char "v");

    reset := true;
    for 0..nmax-1 set(z, it, czero);

    #if DEBUG {
        COUNT := 0;
    }

    for in:0..nidim-1 {
        n := idim[in];
        nk := ifx banded then nkb else 1;

        for iku:0..nk-1 {
            k := ifx banded then kb[iku] else n-1;

            lda := ifx banded then k+1 else n;
            if !packed && lda < nmax then lda += 1;
            if !packed && lda > nmax then continue;

            laa := ifx packed then (n*(n+1))/2 else lda*n;
            _null := n <= 0;

            for icu:0..1 {
                uplo := ichu[icu];

                for ict:0..2 {
                    trans := icht[ict];

                    for icd:0..1 {
                        diag := ichd[icd];

                        // Generate the matrix A.

                        transl := czero;
                        _slice := String.slice(sname, 1, 2);
                        cmake(_slice, uplo, diag, n, n, a, nmax, aa, lda, k, k, *reset, transl);

                        for ix:0..ninc-1 {
                            incx := inc[ix];
                            lx := abs(incx)*n;

                            transl = chalf;
                            x_view := array2d_view(x, 1, row_major=false);
                            cmake("ge", 0, 0, 1, n, *x_view, 1, _xx, abs(incx), 0, n-1, *reset, transl);

                            if n > 1 {
                                set(x, n/2 - 1, czero);
                                set(_xx, abs(incx)*(n/2-1), czero);
                            }

                            // Save every datum  before calling the subroutine.
                            for 0..laa-1 set(as, it, get(aa, it));
                            for 0..lx-1 set(xs, it, get(_xx, it));

                            aa_view := array2d_view(aa, lda, row_major=false);

                            // Call the subroutine
                            if name_mv {
                                if full        _trmv(uplo, trans, diag, n, aa_view, lda, _xx, incx);
                                else if banded _tbmv(uplo, trans, diag, n, k, aa_view, lda, _xx, incx);
                                else if packed _tpmv(uplo, trans, diag, n, aa, _xx, incx);
                            } else if name_sv {
                                if full        _trsv(uplo, trans, diag, n, aa_view, lda, _xx, incx);
                                else if banded _tbsv(uplo, trans, diag, n, k, aa_view, lda, _xx, incx);
                                else if packed _tpsv(uplo, trans, diag, n, aa, _xx, incx);
                            }

                            // See what data chnaged inside subroutines
                            ok := lce(as, aa, laa);
                            pass &= ok;
                            if !ok {
                                err(" lcs | as  = %", as);
                                err(" lcs | aa  = %", aa);
                                err(" lcs | laa = %", laa, exit_program=fatal);
                            }
                            if _null {
                                ok = lce(xs, _xx, lx);
                                pass &= ok;
                                if !ok {
                                    err(" lce | xs = %", xs);
                                    err(" lce | xx = %", _xx);
                                    err(" lce | lx = %", lx, exit_program=fatal);
                                }
                            } else {
                                xs_view := array2d_view(xs, abs(incx), row_major=false);
                                xx_view := array2d_view(_xx, abs(incx), row_major=false);
                                ok = lceres("ge", 0, 1, n, xs_view, xx_view, abs(incx));
                                pass &= ok;
                                if !ok {
                                    err(" lceres | n    = %", n);
                                    err(" lceres | xs   = %", xs);
                                    err(" lceres | xx   = %", _xx);
                                    err(" lceres | incx = %", incx, exit_program=fatal);
                                }
                            }

                            if !_null {
                                if name_mv {
                                    pass &= cmvch(trans, n, n, cone, a, nmax, x, incx, czero, z, incx, xt, g, _xx, eps, fatal);
                                } else if name_sv {
                                    for i:0..n-1 {
                                        set(z, i, get(_xx, i*abs(incx)));
                                        set(_xx, i*abs(incx), get(x, i));
                                    }
                                    pass &= cmvch(trans, n, n, cone, a, nmax, z, incx, czero, x, incx, xt, g, _xx, eps, fatal);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return pass;
}


//  ######   ######  ##     ## ##    ## ##
// ##    ## ##    ## ##     ## ##   ##  ##    ##
// ##       ##       ##     ## ##  ##   ##    ##
// ##       ##       ######### #####    ##    ##
// ##       ##       ##     ## ##  ##   #########
// ##    ## ##    ## ##     ## ##   ##        ##
//  ######   ######  ##     ## ##    ##       ##

cchk4 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nalf: int,
    alf: []$T,
    ninc: int,
    inc: []int,
    nmax: int,
    incmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    x: *$X/ArrayType1D,
    _xx: *$XX/ArrayType1D,
    xs: *$XS/ArrayType1D,
    y: *$Y/ArrayType1D,
    yy: *$YY/ArrayType1D,
    ys: *$YS/ArrayType1D,
    yt: *$YT/ArrayType1D,
    g: *$G/ArrayType1D,
    z: *$Z/ArrayType1D,
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == X.data_type);
        assert(T == XX.data_type);
        assert(T == XS.data_type);
        assert(T == Y.data_type);
        assert(T == YY.data_type);
        assert(T == YS.data_type);
        assert(T == YT.data_type);
        assert(U == G.data_type);
        assert(T == Z.data_type);
    }
    DEBUG :: true;
    println_push("cchk4");

    pass := true;

    czero := zero(T);
    cone  := one(T);
    chalf := ncast(T, 0.5);
    rzero := zero(U);

    conjug := char_eq(sname[4], #char "c");
    nargs := 9;

    nc := 0;
    reset := true;
    errmax := rzero;

    w : Array1D(1, T);

    for in:0..nidim-1 {
        n := idim[in];
        nd := n/2 + 1;

        for im:0..1 {
            m := ifx im == 0 then max(n-nd, 0) else min(n+nd, nmax);

            lda := m;
            if lda < nmax then lda += 1;
            if lda > nmax then continue;

            laa := lda*n;
            _null := n <= 0 || m <= 0;

            for ix:0..ninc-1 {
                incx := inc[ix];
                lx := abs(incx)*m;

                transl := chalf;
                x_view := array2d_view(x, 1, row_major=false);
                cmake("ge", 0, 0, 1, m, *x_view, 1, _xx, abs(incx), 0, m-1, *reset, transl);

                if m > 1 {
                    fset(x, m/2, czero);
                    fset(_xx, 1 + abs(incx)*(m/2-1), czero);
                }

                for iy:0..ninc-1 {
                    incy := inc[iy];
                    ly := abs(incy)*n;

                    transl = czero;
                    y_view := array2d_view(y, 1, row_major=false);
                    cmake("ge", 0, 0, 1, n, *y_view, 1, yy, abs(incy), 0, n-1, *reset, transl);

                    if n > 1 {
                        fset(y, n/2, czero);
                        fset(yy, 1 + abs(incy)*(n/2-1), czero);
                    }

                    for ia:0..nalf-1 {
                        alpha := alf[ia];

                        transl = czero;
                        _slice := String.slice(sname, 1, 2);
                        cmake(_slice, 0, 0, m, n, a, nmax, aa, lda, m-1, n-1, *reset, transl);

                        nc += 1;

                        ms := m;
                        ns := n;
                        als := alpha;
                        for i:0..laa-1 {
                            set(as, i, get(aa, i));
                        }
                        ldas := lda;
                        for i:0..lx-1 {
                            set(xs, i, get(_xx, i));
                        }
                        incxs := incx;
                        for i:0..ly-1 {
                            set(ys, i, get(yy, i));
                        }
                        incys := incy;

                        aa_view := array2d_view(aa, lda, row_major=false);
                        if conjug {
                            cgerc(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                        } else {
                            cgeru(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                        }

                        as_view := array2d_view(as, lda, row_major=false);
                        isame1 := ms == m;
                        isame2 := ns == n;
                        isame3 := als == alpha;
                        isame4 := lce(xs, _xx, lx);
                        isame5 := incxs == incx;
                        isame6 := lce(ys, yy, ly);
                        isame7 := incys == incy;
                        isame8 := ifx _null then lce(as, aa, laa) else lceres("ge", 0, m, n, aa_view, as_view, lda);
                        isame9 := ldas == lda;

                        same := isame1 && isame2 && isame3 && isame4 && isame5 && isame6 && isame7 && isame8 && isame9;
                        pass &= same;
                        if !same {
                            if !isame1 err(" cchk4 | parameter 1 changed incorrectly");
                            if !isame2 err(" cchk4 | parameter 2 changed incorrectly");
                            if !isame3 err(" cchk4 | parameter 3 changed incorrectly");
                            if !isame4 err(" cchk4 | parameter 4 changed incorrectly");
                            if !isame5 err(" cchk4 | parameter 5 changed incorrectly");
                            if !isame6 err(" cchk4 | parameter 6 changed incorrectly");
                            if !isame7 err(" cchk4 | parameter 7 changed incorrectly");
                            if !isame8 err(" cchk4 | parameter 8 changed incorrectly");
                            if !isame9 err(" cchk4 | parameter 9 changed incorrectly");
                            err(" cchk4 | sname = %", sname);
                            err(" cchk4 | m = %, n = %, alpha = %, incx = %, incy = %, lda = %", m, n, alpha, incx, incy, lda, exit_program=fatal);
                        }

                        if !_null {
                            if incx > 0 {
                                for i:1..m {
                                    fset(z, i, fget(x, i));
                                }
                            } else {
                                for i:1..m {
                                    fset(z, i, fget(x, m-i+1));
                                }
                            }

                            z_view := array2d_view(z, nmax, row_major=false);
                            for j:1..n {
                                jy := ifx incy > 0 then 1 + (j-1)*abs(incy) else 1 + (n-j)*abs(incy);
                                fset(*w, 1, fget(yy, jy));
                                if conjug {
                                    fset(*w, 1, conjugate(fget(w, 1)));
                                }

                                for i:1..m {
                                    fset(y,  i, fget(a.*, i, j));
                                    fset(ys, i, fget(aa, i + (j-1)*lda));
                                }

                                ok := cmvch(#char "n", m, 1, alpha, z_view, nmax, w, 1, cone, y, 1, yt, g, ys, eps, fatal);
                                pass &= ok;

                                erri := rzero;
                                for i:1..m {
                                    diff := abs(fget(yt.*, i) - fget(ys, i))/eps;
                                    gi := fget(g.*, i);
                                    if gi != rzero then diff /= gi;
                                    erri = max(erri, diff);
                                }
                                errmax = max(errmax, erri);
                            }
                        }
                    }
                }
            }
        }
    }

    if errmax >= thresh {
        pass = false;
    }

    _ = incmax;
    _ = nargs;
    _ = DEBUG;

    return pass;
}



//  ######   ######  ##     ## ##    ## ########
// ##    ## ##    ## ##     ## ##   ##  ##
// ##       ##       ##     ## ##  ##   ##
// ##       ##       ######### #####    #######
// ##       ##       ##     ## ##  ##         ##
// ##    ## ##    ## ##     ## ##   ##  ##    ##
//  ######   ######  ##     ## ##    ##  ######

cchk5 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nalf: int,
    alf: []$T,
    ninc: int,
    inc: []int,
    nmax: int,
    incmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    x: *$X/ArrayType1D,
    _xx: *$XX/ArrayType1D,
    xs: *$XS/ArrayType1D,
    y: *$Y/ArrayType1D,
    yy: *$YY/ArrayType1D,
    ys: *$YS/ArrayType1D,
    yt: *$YT/ArrayType1D,
    g: *$G/ArrayType1D,
    z: *$Z/ArrayType1D,
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == X.data_type);
        assert(T == XX.data_type);
        assert(T == XS.data_type);
        assert(T == Y.data_type);
        assert(T == YY.data_type);
        assert(T == YS.data_type);
        assert(T == YT.data_type);
        assert(U == G.data_type);
        assert(T == Z.data_type);
    }
    DEBUG :: true;
    println_push("cchk5");

    pass := true;

    ichu :: "ul";

    czero := zero(T);
    cone  := one(T);
    chalf := ncast(T, 0.5);
    rzero := zero(U);
    rone  := one(U);

    full   := char_eq(sname[2], #char "e");
    packed := char_eq(sname[2], #char "p");

    nargs := ifx full then 7 else 6;

    nc := 0;
    reset := true;
    errmax := rzero;

    w : Array1D(1, T);

    for in:0..nidim-1 {
        n := idim[in];
        k := n - 1;

        lda := ifx full then n else 1;
        if full && lda < nmax then lda += 1;
        if lda > nmax then continue;

        laa := ifx full then lda*n else (n*(n+1))/2;

        for iu:0..1 {
            uplo := ichu[iu];
            upper := char_eq(uplo, #char "u");

            for ix:0..ninc-1 {
                incx := inc[ix];
                lx := abs(incx)*n;

                transl := chalf;
                x_view := array2d_view(x, 1, row_major=false);
                cmake("ge", 0, 0, 1, n, *x_view, 1, _xx, abs(incx), 0, n-1, *reset, transl);

                if n > 1 {
                    fset(x, n/2, czero);
                    fset(_xx, 1 + abs(incx)*(n/2-1), czero);
                }

                for ia:0..nalf-1 {
                    ralpha := real(alf[ia]);
                    if imag(alf[ia]) != rzero then continue;

                    alpha := cpx(ralpha, rzero);
                    _null := n <= 0 || ralpha == rzero;

                    transl = czero;
                    _slice := String.slice(sname, 1, 2);
                    cmake(_slice, uplo, 0, n, n, a, nmax, aa, lda, k, k, *reset, transl);

                    nc += 1;

                    for i:0..laa-1 {
                        set(as, i, get(aa, i));
                    }
                    for i:0..lx-1 {
                        set(xs, i, get(_xx, i));
                    }

                    if full {
                        aa_view := array2d_view(aa, lda, row_major=false);
                        cher(uplo, n, ralpha, _xx, incx, *aa_view, lda);
                    } else if packed {
                        chpr(uplo, n, ralpha, _xx, incx, aa);
                    }

                    ok := lce(xs, _xx, lx);
                    pass &= ok;
                    if !ok {
                        err(" lce | xs = %", xs);
                        err(" lce | xx = %", _xx);
                        err(" lce | lx = %", lx, exit_program=fatal);
                    }

                    if _null {
                        ok = lce(as, aa, laa);
                        pass &= ok;
                        if !ok {
                            err(" lce | as  = %", as);
                            err(" lce | aa  = %", aa);
                            err(" lce | laa = %", laa, exit_program=fatal);
                        }
                    } else if full {
                        as_view := array2d_view(as, lda, row_major=false);
                        aa_view := array2d_view(aa, lda, row_major=false);
                        ok = lceres("he", uplo, n, n, aa_view, as_view, lda);
                        pass &= ok;
                        if !ok {
                            err(" lceres | uplo = %", uplo);
                            err(" lceres | n    = %", n);
                            err(" lceres | lda  = %", lda);
                            err(" lceres | as   = %", as_view);
                            err(" lceres | aa   = %", aa_view, exit_program=fatal);
                        }
                    }

                    if !_null {
                        if incx > 0 {
                            for i:1..n {
                                fset(z, i, fget(x, i));
                            }
                        } else {
                            for i:1..n {
                                fset(z, i, fget(x, n-i+1));
                            }
                        }

                        ja := 1;
                        for j:1..n {
                            fset(*w, 1, conjugate(fget(z, j)));
                            jj, lj : int = ---;

                            if upper {
                                jj = 1;
                                lj = j;
                            } else {
                                jj = j;
                                lj = n-j+1;
                            }

                            for i:1..lj {
                                fset(yy, i, fget(z, jj+i-1));
                                fset(y,  i, fget(a.*, jj+i-1, j));
                                fset(ys, i, fget(aa, ja+i-1));
                            }

                            yy_view := array2d_view(yy, lj, row_major=false);
                            ok = cmvch(#char "n", lj, 1, alpha, yy_view, nmax, w, 1, cone, y, 1, yt, g, ys, eps, fatal);
                            pass &= ok;

                            erri := rzero;
                            for i:1..lj {
                                yii := i;
                                gi := fget(g.*, i);
                                diff := abs(fget(yt.*, i) - fget(ys, yii))/eps;
                                if gi != rzero then diff /= gi;
                                erri = max(erri, diff);
                            }
                            errmax = max(errmax, erri);

                            if full {
                                if upper {
                                    ja += lda;
                                } else {
                                    ja += lda + 1;
                                }
                            } else {
                                ja += lj;
                            }
                        }
                    } else {
                        if n <= 0 then continue;
                    }
                }
            }
        }
    }

    if errmax >= thresh {
        pass = false;
    }

    _ = nargs;
    _ = incmax;
    _ = rone;

    return pass;
}



//  ######  ##     ##    ###    ##    ## ########
// ##    ## ###   ###   ## ##   ##   ##  ##
// ##       #### ####  ##   ##  ##  ##   ##
// ##       ## ### ## ##     ## #####    ######
// ##       ##     ## ######### ##  ##   ##
// ##    ## ##     ## ##     ## ##   ##  ##
//  ######  ##     ## ##     ## ##    ## ########

cmake :: (
    type: string,
    uplo: u8,
    diag: u8,
    m: int,
    n: int,
    a: *$A/ArrayType2D,
    nmax: int,
    aa: *$AA/ArrayType1D,
    lda: int,
    kl: int,
    ku: int,
    reset: *bool,
    transl: $T
) {
    #run {
        assert(T == A.data_type);
        assert(T == AA.data_type);
    }
    U :: #run get_number_inner_type(T);

    DEBUG :: false;
    #if DEBUG println_push("cmake: type = %, uplo = %, diag = %, m = %, n = %, kl = %, ku = %, nmax = %, lda = %", type, uplo, diag, m, n, kl, ku, nmax, lda, color=.FG_CYAN);

    czero := zero(T);
    cone  := one(T);
    rzero := zero(U);
    rone  := one(U);

    rogue  := cpx(ncast(U, -1.0e10), ncast(U, 1.0e10));
    rrogue := ncast(U, -1.0e10);

    gen := char_eq(type[0], #char "g");
    sym := char_eq(type[0], #char "h");
    tri := char_eq(type[0], #char "t");
    upper := (sym || tri) && char_eq(uplo, #char "u");
    lower := (sym || tri) && char_eq(uplo, #char "l");
    unit := tri && char_eq(diag, #char "u");

    #if DEBUG println("gen: %, sym: %, tri: %, upper: %, lower: %, unit: %", gen, sym, tri, upper, lower, unit, color=.FG_CYAN);

    // Generate data in array a.
    // println("generate data in array a.");
    for j:1..n {
        for i:1..m {
            if !gen && !(upper && i <= j) && !(lower && i >= j) then continue;

            // FORTRAN: j-i <= ku, i=j <= kl
            if (i <= j && j-i <= ku) || (i >= j && i-j <= kl) {
                fset(a, i, j, cbeg(reset, T)+transl);
            } else {
                fset(a, i, j, czero);
            }

            if i == j then continue;

            if sym {
                fset(a, j, i, conjugate(fget(a.*,i,j)));
            } else {
                fset(a, j, i, czero);
            }
        }
        if sym  then fset(a, j, j, real_nc(fget(a.*,j,j)));
        if tri  then fset(a, j, j, fget(a.*,j,j) + cone);
        if unit then fset(a, j, j, cone);
    }

    // Store elements in array AS in data structure required by routine.
    // println("store elements in array as in data structure required by routine.");
    if String.equal_nocase(type, "ge") {
        #if DEBUG println("ge", color=.FG_CYAN);
        for j:1..n {
            for i:1..m {
                // FORTRAN: i + (j-1)*lda
                // * This assumes a later interpretation of aa as a column-major matrix - Fehse 2022-03-03
                fset(aa, i + (j-1)*lda, fget(a.*,i,j));
            }
            for i:m+1..lda {
                // FORTRAN: i + (j-1)*lda
                // * This assumes a later interpretation of aa as a column-major matrix - Fehse 2022-03-03
                fset(aa, i + (j-1)*lda, rogue);
            }
        }
    } else if String.equal_nocase(type, "gb") {
        #if DEBUG println("gb", color=.FG_CYAN);
        for j:1..n {
            for band_row:1..lda {
                fset(aa, band_row + (j-1)*lda, rogue);
            }

            ibeg := max(1, j-ku);
            iend := min(m, j+kl);
            for i:ibeg..iend {
                band_row := ku + 1 + i - j;
                fset(aa, band_row + (j-1)*lda, fget(a.*, i, j));
            }
        }
    } else if String.equal_nocase(type, "he") || String.equal_nocase(type, "tr") {
        #if DEBUG println("he || tr", color=.FG_CYAN);
        for j:1..n {
            ibeg, iend : int = ---;
            if upper {
                ibeg = 1;
                iend = ifx unit then j-1 else j;
            } else {
                ibeg = ifx unit then j+1 else j;
                iend = n;
            }
            for i:1..ibeg-1 {
                fset(aa, i+(j-1)*lda, rogue);
            }
            for i:ibeg..iend {
                fset(aa, i+(j-1)*lda, fget(a.*,i,j));
            }
            for i:iend+1..lda {
                fset(aa, i+(j-1)*lda, rogue);
            }
            if sym {
                jj := j + (j-1)*lda;
                fset(aa,jj, cpx(real(fget(aa.*,jj)), rrogue));
            }
        }
    } else if String.equal_nocase(type, "hb") || String.equal_nocase(type, "tb") {
        #if DEBUG println("hb || tb", color=.FG_CYAN);
        for j:1..n {
            kk, ibeg, iend : int = ---;
            if upper {
                kk = kl + 1;
                ibeg = max(1, kk-j+1);
                iend = ifx unit then kk-1 else kk;
            } else {
                kk = 1;
                ibeg = ifx unit then 2 else 1;
                iend = min(kl+1, m-j+1);
            }
            for i:1..ibeg-1 {
                fset(aa, i+(j-1)*lda, rogue);
            }
            for i:ibeg..iend {
                fset(aa, i+(j-1)*lda, fget(a.*,i+j-kk,j));
            }
            for i:iend+1..lda {
                fset(aa, i+(j-1)*lda, rogue);
            }
            if sym {
                jj := kk + (j-1)*lda;
                fset(aa,jj, cpx(real(fget(aa.*,jj)), rrogue));
            }
        }
    } else if String.equal_nocase(type, "hp") || String.equal_nocase(type, "tp") {
        #if DEBUG println("hp || tp", color=.FG_CYAN);
        ioff := 1;
        for j:1..n {
            ibeg, iend : int = ---;
            if upper {
                ibeg = 1;
                iend = j;
            } else {
                ibeg = j;
                iend = n;
            }
            for i:ibeg..iend {
                fset(aa, ioff, fget(a.*,i,j));
                if i == j {
                    if sym {
                        fset(aa,ioff, cpx(real(fget(aa.*,ioff)), rrogue));
                    } else {
                        fset(aa,ioff, fget(a.*,i,j));
                    }
                }
                ioff += 1;
            }
        }
    }
}


//  ######  ########  ########  ######
// ##    ## ##     ## ##       ##    ##
// ##       ##     ## ##       ##
// ##       ########  ######   ##   ####
// ##       ##     ## ##       ##    ##
// ##    ## ##     ## ##       ##    ##
//  ######  ########  ########  ######


CBEG_Save :: struct {
    mi : int = 891;
    mj : int = 457;
    i  : int = 7;
    j  : int = 7;
    ic : int = 0;
}
cbeg_save : CBEG_Save;


// Generates complex number as pairs of random numbers uniformly
// distributed between -0.5 and 0.5.

// Written on 10-August-1987.
// Richard Hanson, Sandia National Labs.
// Jeremy Du Croz, NAG Central Office.
cbeg :: (
    reset: *bool,
    $T: Type
) ->  T {
    using cbeg_save;

    DEBUG :: false;
    #if DEBUG println_push("cbeg(%,%)", reset.*, T, color=.FG_CYAN);

    if reset.* {
        // Initialize local variables
        mi = 891;
        mj = 457;
        i = 7;
        j = 7;
        ic = 0;
        reset.* = false;
    }

    // The sequence of values of I or J is bounded between 1 and 999.
    // If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
    // If initial I or J = 4 or 8, the period will be 25.
    // If initial I or J = 5, the period will be 10.
    // IC is used to break up the period by skipping 1 value of I or J in 6.

    ic += 1;
    while true {
        i *= mi;
        j *= mj;
        i -= 1000*(i/1000);
        j -= 1000*(j/1000);
    
        if ic < 5 then break;

        ic = 0;
    }
    X :: #run get_number_inner_type(T);
    return cpx( ncast(X, cast(float64)(i-500)/1001.0), ncast(X, cast(float64)(j-500)/1001.0) );
}
