// ? LAPACK cblat2
// ? http://www.netlib.org/lapack/explore-html/d2/d54/cblat2_8f_source.html


// ? Does is make sense to accurately copy the behavior of the test suite?
// ? Meaning, generating output files etc.?


file_url : string : "../LAPACK/Testing/cblat2.in";

cblat2 :: () -> bool {
    println_push("cblat2", color=.FG_GREEN);

    data := read_file_and_split(file_url);
    defer {
        for data { array_free(it); }
        array_free(data);
    };

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
    xs, _xx, yx, yy : Array1D(nmax*incmax, C32) = ---;
    alf : Array1D(nalmax, C32) = ---;
    bet : Array1D(nbemax, C32) = ---;
    g : Array1D(nmax, float32) = ---;

    sfac : float32 = ---;
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


    // parse data
    {
        sfac = String.parse_float(*data[7][0]);
        println("sfac = %", sfac);
        
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
    if !pass then err("ERROR: cmvch result not good: n = %\n%\n%", n, yy, yt, exit_program=true);

    trans = #char "T";
    pass &= cmvch(trans, n, n, cone, a, nmax, x, -1, czero, y, -1, *yt, *g, yy, eps, true);
    pass &= lce(yy, yt, n);
    if !pass then err("ERROR: cmvch result not good: n = %\n%\n%", n, yy, yt, exit_program=true);


    // Test each subroutine in turn


    array_free(idim);
    array_free(kb);
    array_free(inc);
    return false; // TODO
}



cmvch :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$alpha: $T,
    a: $A/ArrayType2D,
    $$nmax: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $$beta: T,
    y: $Y/ArrayType1D,
    $$incy: int,
    yt: *$YT/ArrayType1D,
    g: *$G/ArrayType1D,
    yy: $YY/ArrayType1D,
    $$eps: $U,
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