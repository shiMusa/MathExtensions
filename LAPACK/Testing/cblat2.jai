// ? LAPACK cblat2
// ? http://www.netlib.org/lapack/explore-html/d2/d54/cblat2_8f_source.html


// ? Does is make sense to accurately copy the behavior of the test suite?
// ? Meaning, generating output files etc.?


file_url : string : "../LAPACK/Testing/cblat2.in";

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
    a : Array2D(nmax, nmax, C32) = ---;
    aa, as : Array1D(nmax*nmax, C32) = ---;
    x, y, yt : Array1D(nmax, C32) = ---;
    z : Array1D(2*nmax, C32) = ---;
    xs, _xx, ys, yy : Array1D(nmax*incmax, C32) = ---;
    alf : Array1D(nalmax, C32) = ---;
    bet : Array1D(nbemax, C32) = ---;
    g   : Array1D(nmax, float32) = ---;

    thresh : float32 = ---;
    nidim, nkb, ninc, nalf, nbet : int = ---;
    idim, kb, inc : []int = ---;
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
                        <<it = true;
                    }
                }
            }
        }
        println("ltest = %", ltest);
    }


    eps := EPS_32; // Complex(float32)

    // Check the riliability of CMVCH using exact data.
    n := min(32, nmax);
    for j:0..n-1 {
        for i:0..n-1 {
            set(*a, i,j, ncast(C32,max(i-j+1, 0)));
        }
        set(*x, j, ncast(C32,j));
        set(*y, j, czero);
    }
    for j:0..n-1 {
        set(*yy, j, ncast(C32, j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3) );
    }
    // YY holds the exact result. On exit from CMVCH YT holds 
    // the result computed by CMVCH.

    pass := true;

    trans : u8 = #char "N";
    pass &= cmvch(trans, n, n, cone, a, nmax, x,  1, czero, y,  1, *yt, *g, yy, eps, true);
    pass &= lce(yy, yt, n);
    if !pass then err("cmvch result not good: n = %\n%\n%", n, yy, yt, exit_program=true);

    trans = #char "T";
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
            pass &= cchk1(
                srnamt, eps, thresh,
                nidim, idim, nkb, kb, nalf, alf.data,
                nbet, bet.data, ninc, inc, nmax, incmax, *a, *aa, as,
                x, _xx, xs, y, yy, ys, yt, g
            );
        } else if isnum <= 4 {

        } else if isnum <= 10 {

        } else if isnum <= 12 {

        } else if isnum <= 14 {

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

    println_push("cmvch");

    czero := zero(T);
    rzero := zero(U);
    rone  := one(U);

    tran  := char_eq(trans, #char "t");
    ctran := char_eq(trans, #char "c");
    ml := ifx tran || ctran then n else m;
    nl := ifx tran || ctran then m else n;
    kx    := ifx incx < 0 then nl-1 else 0;
    incxl := ifx incx < 0 then -1 else 1;
    ky    := ifx incy < 0 then ml-1 else 0;
    incyl := ifx incy < 0 then -1 else 1;
    iy, jx: int = ---;
    erri: U = ---;

    // Compute expected result in YT using data in A, X and Y.
    // Compute gauges in G.
    iy = ky;
    for i:0..ml-1 {
        set(yt, iy, czero);
        set(g,iy, rzero);
        jx = kx;
        if tran {
            for j:0..nl-1 {
                set(yt, iy, get(yt,iy) + get(a,j,i)*get(x,jx));
                set(g,iy, get(g,iy) + abs1(get(a,j,i))*abs1(get(x,jx)));
                jx += incxl;
            }
        } else if ctran {
            for j:0..nl-1 {
                set(yt, iy, get(yt,iy) + conjugate(get(a,j,i))*get(x,jx));
                set(g,iy, get(g,iy) + abs1(get(a,j,i))*abs1(get(x,jx)));
                jx += incxl;
            }
        } else {
            for j:0..nl-1 {
                set(yt, iy, get(yt,iy) + get(a,i,j)*get(x,jx));
                set(g,iy, get(g,iy) + abs1(get(a,i,j))*abs1(get(x,jx)));
                jx += incxl;
            }
        }
        set(yt,iy, alpha*get(yt,iy) + beta*get(y,iy));
        set(g,iy, abs1(alpha)*get(g,iy) + abs1(beta)*abs1(get(y,iy)));
        iy += incyl;
    }

    // Compute the error ration for this result.
    rerr := rzero;
    res := true;
    for i:0..ml-1 {
        erri = abs(get(yt,i) - get(yy,i*abs(incy)))/eps;
        if get(g,i) != rzero {
            erri /= get(g,i);
        }
        rerr = max(rerr,erri);
        if rerr*sqrt(eps) >= rone {
            err("%*sqrt(%) = % >= 1.0", rerr, eps, rerr*sqrt(eps), exit_program=fatal);
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
    for 0..lr {
        if get(ri,it) != get(rj,it) then return false;
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
    as: $AS/ArrayType1D,
    x: $X/ArrayType1D,
    _xx: $XX/ArrayType1D,
    xs: $XS/ArrayType1D,
    y: $Y/ArrayType1D,
    yy: $YY/ArrayType1D,
    ys: $YS/ArrayType1D,
    yt: $YT/ArrayType1D,
    g: $G/ArrayType1D
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
    println_push("cchk1");

    czero := zero(C32);
    cone  := one(C32);
    rzero := zero(U);
    rone  := one(U);

    full   := char_eq(sname[2], #char "e");
    banded := char_eq(sname[2], #char "b");

    // Define the number of arguments
    nargs := ifx full then 11 else 13;

    nc := 0;
    reset := true;
    errmax := rzero;

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
                _null := n < 0 || m < 0;

                // Generate the matrix a.
                transl := czero;
                _slice := String.slice(sname, 1, 2);
                cmake(_slice, 0, 0, m, n, a, nmax, aa, lda, kl, ku, *reset, transl);

                // TODO finish
            }
        }
    }

    err("not yet fully implemented");

    return false; // TODO
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

    // println_push("cmake: type = %, uplo = %, diag = %, m = %, n = %, kl = %, ku = %, nmax = %, lda = %", type, uplo, diag, m, n, kl, ku, nmax, lda);

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

    // Generate data in array a.
    // println("generate data in array a.");
    for j:0..n-1 {
        for i:0..m-1 {
            if !gen && !(upper && i <= j) && !(lower && i >= j) then continue;

            if (i <= j && j-i <= ku) || (i >= j && i-j <= kl) {
                set(a,i,j, cbeg(reset, T)+transl);
            } else {
                set(a,i,j, czero);
            }

            if i == j then continue;

            if sym {
                set(a,j,i, conjugate(get(a,i,j)));
            } else {
                set(a,j,i, czero);
            }
        }
        if sym  then set(a,j,j, real_nc(get(a,j,j)));
        if tri  then set(a,j,j, get(a,j,j) + cone);
        if unit then set(a,j,j, cone);
    }

    // Store elements in array AS in data structure required by routine.
    // println("store elements in array as in data structure required by routine.");
    if String.equal_nocase(type, "ge") {
        for j:0..n-1 {
            for i:0..m-1 {
                set(aa, i+j*lda, get(a,i,j));
            }
            for i:m..lda-1 {
                set(aa, i+j*lda, rogue);
            }
        }
    } else if String.equal_nocase(type, "gb") {
        for j:0..n-1 {
            i1, i2, i3 := 0;
            for 0..ku-j {
                set(aa, i1+j*lda, rogue);
                i1 += 1;
            }
            i2 = i1;
            for i1..min(kl+ku, ku+m-j) {
                set(aa, i2+j*lda, get(a, i2+j-ku-1, j));
                i2 += 1;
            }
            for i3:i2..lda-1 {
                set(aa, i3+j*lda, rogue);
            }
        }
    } else if String.equal_nocase(type, "he") || String.equal_nocase(type, "tr") {
        for j:0..n-1 {
            ibeg, iend : int = ---;
            if upper {
                ibeg = 0;
                iend = ifx unit then j-1 else j;
            } else {
                ibeg = ifx unit then j+1 else j;
                iend = n;
            }
            for i:0..ibeg-1 {
                set(aa, i+j*lda, rogue);
            }
            for i:ibeg..iend {
                set(aa, i+j*lda, get(a,i,j));
            }
            for i:iend+1..lda-1 {
                set(aa, i+j*lda, rogue);
            }
            if sym {
                jj := j + j*lda;
                set(aa,jj, cpx(real(get(aa,jj)), rrogue));
            }
        }
    } else if String.equal_nocase(type, "hb") || String.equal_nocase(type, "tb") {
        for j:0..n-1 {
            kk, ibeg, iend : int = ---;
            if upper {
                kk = kl; // correct Fortran indexing
                ibeg = max(0, kl-j); // correct Fortran indexing
                iend = ifx unit then kl-1 else kl; // correct Fortran indexing
            } else {
                kk = 0; // correct Fortran indexing
                ibeg = ifx unit then 1 else 0; // correct Fortran indexing
                iend = min(kl, m-j); // correct Fortran indexing
            }
            for i:0..ibeg-1 {
                set(aa, i+j*lda, rogue);
            }
            for i:ibeg..iend {
                set(aa, i+j*lda, get(a,i+j-kk,j));
            }
            for i:iend+1..lda-1 {
                set(aa, i+j*lda, rogue);
            }
            if sym {
                jj := kk + j*lda;
                set(aa,jj, cpx(real(get(aa,jj)), rrogue));
            }
        }
    } else if String.equal_nocase(type, "hp") || String.equal_nocase(type, "tp") {
        ioff := 0;
        ibeg, iend : int = ---;
        for j:0..n-1 {
            if upper {
                ibeg = 0;
                iend = j-1;
            } else {
                ibeg = j-1;
                iend = n-1;
            }
            for i:ibeg..iend {
                ioff += 1;
                set(aa, ioff, get(a,i,j));
                if i == j {
                    if unit then 
                        set(aa,ioff, get(a,i,j));
                    else
                        set(aa,ioff, cpx(real(get(aa,ioff)), rrogue));
                }
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
    i  : int = 6;
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

    if <<reset {
        // Initialize local variables
        mi = 891;
        mj = 457;
        i = 7;
        j = 7;
        ic = 0;
        <<reset = false;
    }

    // The sequence of values of I or J is bounded between 1 and 999.
    // If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
    // If initial I or J = 4 or 8, the period will be 25.
    // If initial I or J = 5, the period will be 10.
    // IC is used to break up the period by skipping 1 value of I or J in 6.

    ic += 1;
    while true {
        using cbeg_save;
        i *= mi;
        j *= mj;
        i -= 1000*(i/1000);
        j -= 1000*(j/1000);
    
        if ic < 5 then break;

        ic = 0;
    }
    X :: #run get_number_inner_type(T);
    return cpx( ncast(X, (i-500)/1001.0), ncast(X, (j-500)/1001.0) );
}