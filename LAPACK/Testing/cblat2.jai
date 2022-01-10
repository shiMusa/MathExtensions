// ? LAPACK cblat2
// ? http://www.netlib.org/lapack/explore-html/d2/d54/cblat2_8f_source.html


// ? Does is make sense to accurately copy the behavior of the test suite?
// ? Meaning, generating output files etc.?


file_url : string : "LAPACK/Testing/cblat2.in";

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

    nin    : int : 5;
    nsubs  : int : 17;

    nmax   : int : 65;
    incmax : int : 2;

    ninmax : int : 7;
    nidmax : int : 9;
    nkbmax : int : 7;
    nalmax : int : 7;
    nbemax : int : 7;

    zero :: #run zero(C32);
    one  :: #run one(C32);

    sfac : float32 = ---;
    nidim, nkb, ninc, nalf, nbet : int = ---;
    idim, kb, inc : []int = ---;
    alf, bet : []C32 = ---;
    snames : [nsubs]string : string.[
        "CGEMV", "CGBMV", "CHEMV", "CHBMV",
        "CHPMB", "CTRMV", "CTBMV", "CTPMV",
        "CTRSV", "CTBSV", "CTPSV", "CGERC",
        "CGERU", "CHER" , "CHPR" , "CHER2",
        "CHPR2"
    ];
    ltest : [nsubs]bool; // default to false

    // parse data
    sfac = String.parse_float(*data[7][0]);
    println("sfac = %", sfac);
    
    nidim = String.parse_int(*data[8][0]);
    assert(nidim > 0 && nidim <= nidmax, "nidim out of range: 0 < % <= %\n", nidim, nidmax);
    println("nidim = %", nidim);
    
    idim = NewArray(nidim, int, false);
    defer array_free(idim);
    for 0..nidim-1 {
        idim[it] = String.parse_int(*data[9][it]);
        assert(idim[it] >= 0 && idim[it] <= nmax, "idim[%] out of range: % <= % <= %\n", it, 0, idim[it], nmax);
    }
    println("idim = %", idim);

    nkb = String.parse_int(*data[10][0]);
    assert(nkb > 0 && nkb <= nkbmax, "nkb out of range: % < % <= %\n", 0, nkb, nkbmax);
    println("nkb = %", nkb);

    kb = NewArray(nkb, int, false);
    defer array_free(kb);
    for 0..nkb-1 {
        kb[it] = String.parse_int(*data[11][it]);
        assert(kb[it] >= 0, "kb[%] out of range: % < 0\n", it, kb[it]);
    }
    println("kb = %", kb);

    ninc = String.parse_int(*data[12][0]);
    assert(ninc >= 1 && ninc <= ninmax, "ninc out of range: 0 < % <= %\n", ninc, ninmax);
    println("ninc = %", ninc);

    inc = NewArray(ninc, int, false);
    defer array_free(inc);
    for 0..ninc-1 {
        inc[it] = String.parse_int(*data[13][it]);
        assert(inc[it] != 0 && abs(inc[it]) <= incmax, "inc[%] out of range: 0 != %, % <= %", it, inc[it], inc[it], incmax);
    }
    println("inc = %", inc);

    nalf = String.parse_int(*data[14][0]);
    assert(nalf > 0 && nalf <= nalmax, "nalf out of range: 0 < % <= %\n", nalf, nalmax);
    println("nalf = %", nalf);

    alf = NewArray(nalf, C32, false);
    defer array_free(alf);
    for 0..nalf-1 {
        // complex number!
        c_num := String.slice(data[15][it], 1, data[15][it].count-2); // remove "(" and ")"
        f_vals := String.split(c_num, ",");
        alf[it] = c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1]));
    }
    println("alf = %", alf);

    nbet = String.parse_int(*data[16][0]);
    assert(nbet > 0 && nbet <= nalmax, "nbet out of range: 0 < % <= %\n", nbet, nalmax);
    println("nbet = %", nbet);

    bet = NewArray(nbet, C32, false);
    defer array_free(bet);
    for 0..nbet-1 {
        // complex number!
        c_num := String.slice(data[17][it], 1, data[17][it].count-2); // remove "(" and ")"
        f_vals := String.split(c_num, ",");
        bet[it] = c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1]));
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





    eps := EPS_32; // Complex(float32)

    // Check the riliability of CMVCH using exact data.
    n = min(32, nmax);
    for j:0..n-1 {
        for i:0..n-1 {
            set(*a, i,j, max(i-j+1, 0));
        }
        set(*x, j, j);
        set(*y, j, zero);
    }
    for j:0..n-1 {
        set(*yy, j, j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3 );
    }
    // YY holds the exact result. On exit from CMVCH YT holds 
    // the result computed by CMVCH.





    return false; // TODO
}