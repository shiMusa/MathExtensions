// ? LAPACK cblat2
// ? 

file_url : string : "LAPACK/Testing/cblat2.in";

cblat2 :: () -> bool {

    data := read_file_and_split(file_url);

    println_pause();

    for line, li: data {
        print("%: ", li);
        for line {
            print("% ", it);
        }
        print("\n");
    }

    nin    : int : 5;
    nsubs  : int : 17;

    nmax   : int : 65;
    incmax : int : 2;

    ninmax : int : 7;
    nidmax : int : 9;
    nkbmax : int : 7;
    nalmax : int : 7;
    nbemax : int : 7;

    sfac : float32 = ---;
    nidim, nkb, ninc, nalf, nbet : int = ---;
    idim, kb, inc : []int = ---;
    alf, bet : []C32 = ---;

    // read data
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

    println_unpause();

    return false; // TODO
}