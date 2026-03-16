// ? LAPACK cblat3
// ? http://www.netlib.org/lapack/explore-html/d9/ddb/cblat3_8f_source.html


cblat3_file_url : string : "LAPACK/Testing/cblat3.in";


cblat3 :: () -> bool {
    println_push("cblat3", color=.FG_GREEN);

    data := read_file_and_split(cblat3_file_url);

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
    nsubs  : int : 9;

    nmax   : int : 65;
    nidmax : int : 9;
    nalmax : int : 7;
    nbemax : int : 7;

    a  : Array2D(nmax, nmax, C32);
    b  : Array2D(nmax, nmax, C32);
    c  : Array2D(nmax, nmax, C32);
    aa : Array1D(nmax*nmax, C32);
    as : Array1D(nmax*nmax, C32);
    bb : Array1D(nmax*nmax, C32);
    bs : Array1D(nmax*nmax, C32);
    cc : Array1D(nmax*nmax, C32);
    cs : Array1D(nmax*nmax, C32);
    ct : Array1D(nmax, C32);
    g  : Array1D(nmax, float32);
    w  : Array1D(2*nmax, C32);

    thresh : float32;
    nidim, nalf, nbet : int;
    idim : []int;
    alf  : Array1D(nalmax, C32);
    bet  : Array1D(nbemax, C32);
    ltest : [nsubs]bool;
    snames : [nsubs]string : string.[
        "CGEMM", "CHEMM", "CSYMM",
        "CTRMM", "CTRSM", "CHERK",
        "CSYRK", "CHER2K", "CSYR2K"
    ];

    {
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

        nalf = String.parse_int(*data[10][0]);
        assert(nalf > 0 && nalf <= nalmax, "nalf out of range: 0 < % <= %\n", nalf, nalmax);
        println("nalf = %", nalf);

        for 0..nalf-1 {
            c_num := String.slice(data[11][it], 1, data[11][it].count-2);
            f_vals := String.split(c_num, ",");
            fset(*alf, it+1, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])));
        }
        println("alf = %", alf);

        nbet = String.parse_int(*data[12][0]);
        assert(nbet > 0 && nbet <= nbemax, "nbet out of range: 0 < % <= %\n", nbet, nbemax);
        println("nbet = %", nbet);

        for 0..nbet-1 {
            c_num := String.slice(data[13][it], 1, data[13][it].count-2);
            f_vals := String.split(c_num, ",");
            fset(*bet, it+1, c32(String.parse_float(*f_vals[0]), String.parse_float(*f_vals[1])));
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

    eps := EPS_32;
    czero := zero(C32);
    cone  := one(C32);
    rzero := zero(float32);

    {
        ab_exact : Array2D(nmax, nmax, C32);
        bb_exact : Array2D(nmax, nmax, C32);
        c_exact  : Array2D(nmax, nmax, C32);
        cc_exact : Array1D(nmax, C32);

        n := min(32, nmax);
        for j:1..n {
            for i:1..n {
                fset(*ab_exact, i, j, ncast(C32, max(i-j+1, 0)));
            }
            fset(*bb_exact, j, 1, ncast(C32, j));
            fset(*bb_exact, 1, j, ncast(C32, j));
            fset(*c_exact, j, 1, czero);
        }
        for j:1..n {
            fset(*cc_exact, j, ncast(C32, j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3));
        }

        cc_exact_view := array2d_view(*cc_exact, nmax, row_major=false);

        errv : float32 = rzero;
        pass := cblat3_cmmch(#char "n", #char "n", n, 1, n, cone, ab_exact, nmax, bb_exact, nmax, czero, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != rzero then err("cmmch result not good: transa = n, transb = n", exit_program=true);

        errv = rzero;
        pass = cblat3_cmmch(#char "n", #char "c", n, 1, n, cone, ab_exact, nmax, bb_exact, nmax, czero, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != rzero then err("cmmch result not good: transa = n, transb = c", exit_program=true);

        for j:1..n {
            fset(*bb_exact, j, 1, ncast(C32, n-j+1));
            fset(*bb_exact, 1, j, ncast(C32, n-j+1));
        }
        for j:1..n {
            fset(*cc_exact, n-j+1, ncast(C32, j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3));
        }

        errv = rzero;
        pass = cblat3_cmmch(#char "c", #char "n", n, 1, n, cone, ab_exact, nmax, bb_exact, nmax, czero, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != rzero then err("cmmch result not good: transa = c, transb = n", exit_program=true);

        errv = rzero;
        pass = cblat3_cmmch(#char "c", #char "c", n, 1, n, cone, ab_exact, nmax, bb_exact, nmax, czero, c_exact, nmax, *ct, *g, cc_exact_view, nmax, eps, *errv, true, true);
        pass &= lce(cc_exact, ct, n);
        if !pass || errv != rzero then err("cmmch result not good: transa = c, transb = c", exit_program=true);
    }

    pass := true;
    for isnum:0..nsubs-1 {
        if !ltest[isnum] then continue;

        srnamt := snames[isnum];
        println("%: %", isnum, srnamt);

        if isnum == 0 {
            pass &= cblat3_cchk1(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g
            );
        } else if isnum == 1 || isnum == 2 {
            pass &= cblat3_cchk2(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g
            );
        } else if isnum == 3 || isnum == 4 {
            pass &= cblat3_cchk3(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *ct, *g, *c
            );
        } else if isnum == 5 || isnum == 6 {
            pass &= cblat3_cchk4(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g
            );
        } else if isnum == 7 || isnum == 8 {
            pass &= cblat3_cchk5(
                srnamt, eps, thresh, true,
                nidim, idim, nalf, alf.data, nbet, bet.data, nmax,
                *a, *aa, *as, *b, *bb, *bs, *c, *cc, *cs, *ct, *g, *w
            );
        } else {
            warn("cblat3 | % not yet implemented, skipping for now", srnamt);
        }
    }

    if !pass then err("cblat3 did not pass.");
    return pass;
}


//  ######  ##     ## ##     ##  ######  ##     ##
// ##    ## ###   ### ###   ### ##    ## ##     ##
// ##       #### #### #### #### ##       ##     ##
// ##       ## ### ## ## ### ## ##       #########
// ##       ##     ## ##     ## ##       ##     ##
// ##    ## ##     ## ##     ## ##    ## ##     ##
//  ######  ##     ## ##     ##  ######  ##     ##

cblat3_cmmch :: (
    transa: u8,
    transb: u8,
    m: int,
    n: int,
    kk: int,
    alpha: $T,
    a: $A/ArrayType2D,
    lda: int,
    b: $B/ArrayType2D,
    ldb: int,
    beta: T,
    c: $C/ArrayType2D,
    ldc: int,
    ct: *$CT/ArrayType1D,
    g: *$G/ArrayType1D,
    cc: $CC/ArrayType2D,
    ldcc: int,
    eps: $U,
    errv: *U,
    $fatal: bool,
    $mv: bool
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == B.data_type);
        assert(T == C.data_type);
        assert(T == CT.data_type);
        assert(U == G.data_type);
        assert(T == CC.data_type);
    }

    czero := zero(T);
    rzero := zero(U);
    rone  := one(U);

    trana  := char_eq(transa, #char "t") || char_eq(transa, #char "c");
    tranb  := char_eq(transb, #char "t") || char_eq(transb, #char "c");
    ctrana := char_eq(transa, #char "c");
    ctranb := char_eq(transb, #char "c");

    errv.* = rzero;

    for j:1..n {
        for i:1..m {
            fset(ct, i, czero);
            fset(g,  i, rzero);
        }

        if !trana && !tranb {
            for k:1..kk {
                for i:1..m {
                    fset(ct, i, fget(ct.*, i) + fget(a, i, k)*fget(b, k, j));
                    fset(g,  i, fget(g.*, i) + abs1(fget(a, i, k))*abs1(fget(b, k, j)));
                }
            }
        } else if trana && !tranb {
            if ctrana {
                for k:1..kk {
                    for i:1..m {
                        fset(ct, i, fget(ct.*, i) + conjugate(fget(a, k, i))*fget(b, k, j));
                        fset(g,  i, fget(g.*, i) + abs1(fget(a, k, i))*abs1(fget(b, k, j)));
                    }
                }
            } else {
                for k:1..kk {
                    for i:1..m {
                        fset(ct, i, fget(ct.*, i) + fget(a, k, i)*fget(b, k, j));
                        fset(g,  i, fget(g.*, i) + abs1(fget(a, k, i))*abs1(fget(b, k, j)));
                    }
                }
            }
        } else if !trana && tranb {
            if ctranb {
                for k:1..kk {
                    for i:1..m {
                        fset(ct, i, fget(ct.*, i) + fget(a, i, k)*conjugate(fget(b, j, k)));
                        fset(g,  i, fget(g.*, i) + abs1(fget(a, i, k))*abs1(fget(b, j, k)));
                    }
                }
            } else {
                for k:1..kk {
                    for i:1..m {
                        fset(ct, i, fget(ct.*, i) + fget(a, i, k)*fget(b, j, k));
                        fset(g,  i, fget(g.*, i) + abs1(fget(a, i, k))*abs1(fget(b, j, k)));
                    }
                }
            }
        } else {
            if ctrana {
                if ctranb {
                    for k:1..kk {
                        for i:1..m {
                            fset(ct, i, fget(ct.*, i) + conjugate(fget(a, k, i))*conjugate(fget(b, j, k)));
                            fset(g,  i, fget(g.*, i) + abs1(fget(a, k, i))*abs1(fget(b, j, k)));
                        }
                    }
                } else {
                    for k:1..kk {
                        for i:1..m {
                            fset(ct, i, fget(ct.*, i) + conjugate(fget(a, k, i))*fget(b, j, k));
                            fset(g,  i, fget(g.*, i) + abs1(fget(a, k, i))*abs1(fget(b, j, k)));
                        }
                    }
                }
            } else {
                if ctranb {
                    for k:1..kk {
                        for i:1..m {
                            fset(ct, i, fget(ct.*, i) + fget(a, k, i)*conjugate(fget(b, j, k)));
                            fset(g,  i, fget(g.*, i) + abs1(fget(a, k, i))*abs1(fget(b, j, k)));
                        }
                    }
                } else {
                    for k:1..kk {
                        for i:1..m {
                            fset(ct, i, fget(ct.*, i) + fget(a, k, i)*fget(b, j, k));
                            fset(g,  i, fget(g.*, i) + abs1(fget(a, k, i))*abs1(fget(b, j, k)));
                        }
                    }
                }
            }
        }

        for i:1..m {
            fset(ct, i, alpha*fget(ct.*, i) + beta*fget(c, i, j));
            fset(g,  i, abs1(alpha)*fget(g.*, i) + abs1(beta)*abs1(fget(c, i, j)));
        }

        for i:1..m {
            erri := abs1(fget(ct.*, i) - fget(cc, i, j))/eps;
            if fget(g.*, i) != rzero {
                erri /= fget(g.*, i);
            }
            errv.* = max(errv.*, erri);
            if errv.* * sqrt(eps) >= rone {
                if mv {
                    err("cmmch | i = %, expected = %, computed = %", i, fget(ct.*, i), fget(cc, i, j));
                } else {
                    err("cmmch | i = %, expected = %, computed = %", i, fget(cc, i, j), fget(ct.*, i));
                }
                if n > 1 {
                    err("cmmch | column = %", j);
                }
                err("cmmch | transa = %, transb = %, m = %, n = %, kk = %", char_to_string(transa), char_to_string(transb), m, n, kk);
                err("cmmch | alpha = %, beta = %, err = %", alpha, beta, errv.*);
                err("", exit_program=fatal);
                return false;
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

cblat3_cchk1 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nalf: int,
    alf: []$T,
    nbet: int,
    bet: []T,
    nmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    b: *$B/ArrayType2D,
    bb: *$BB/ArrayType1D,
    bs: *$BS/ArrayType1D,
    c: *$C/ArrayType2D,
    cc: *$CC/ArrayType1D,
    cs: *$CS/ArrayType1D,
    ct: *$CT/ArrayType1D,
    g: *$G/ArrayType1D
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == B.data_type);
        assert(T == BB.data_type);
        assert(T == BS.data_type);
        assert(T == C.data_type);
        assert(T == CC.data_type);
        assert(T == CS.data_type);
        assert(T == CT.data_type);
        assert(U == G.data_type);
    }

    println_push("cchk1");

    ich :: "NTC";

    pass := true;
    reset := true;
    errmax := zero(U);
    czero := zero(T);

    nc : int = 0;

    for im:1..nidim {
        m := idim[im-1];

        for in:1..nidim {
            n := idim[in-1];
            ldc := m;
            if ldc < nmax then ldc += 1;
            if ldc > nmax then continue;
            lcc := ldc*n;
            is_null := n <= 0 || m <= 0;

            for ik:1..nidim {
                k := idim[ik-1];

                for ica:1..3 {
                    transa := ich[ica-1];
                    trana := char_eq(transa, #char "t") || char_eq(transa, #char "c");

                    ma, na : int = ---;
                    if trana {
                        ma = k;
                        na = m;
                    } else {
                        ma = m;
                        na = k;
                    }

                    lda := ma;
                    if lda < nmax then lda += 1;
                    if lda > nmax then continue;
                    laa := lda*na;

                    cmake("ge", 0, 0, ma, na, a, nmax, aa, lda, ma-1, na-1, *reset, czero);

                    for icb:1..3 {
                        transb := ich[icb-1];
                        tranb := char_eq(transb, #char "t") || char_eq(transb, #char "c");

                        mb, nb : int = ---;
                        if tranb {
                            mb = n;
                            nb = k;
                        } else {
                            mb = k;
                            nb = n;
                        }

                        ldb := mb;
                        if ldb < nmax then ldb += 1;
                        if ldb > nmax then continue;
                        lbb := ldb*nb;

                        cmake("ge", 0, 0, mb, nb, b, nmax, bb, ldb, mb-1, nb-1, *reset, czero);

                        for ia:1..nalf {
                            alpha := alf[ia-1];

                            for ib:1..nbet {
                                beta := bet[ib-1];

                                cmake("ge", 0, 0, m, n, c, nmax, cc, ldc, m-1, n-1, *reset, czero);

                                nc += 1;

                                transas := transa;
                                transbs := transb;
                                ms := m;
                                ns := n;
                                ks := k;
                                als := alpha;
                                for i:1..laa {
                                    fset(as, i, fget(aa.*, i));
                                }
                                ldas := lda;
                                for i:1..lbb {
                                    fset(bs, i, fget(bb.*, i));
                                }
                                ldbs := ldb;
                                bls := beta;
                                for i:1..lcc {
                                    fset(cs, i, fget(cc.*, i));
                                }
                                ldcs := ldc;

                                aa_view := array2d_view(aa, lda, row_major=false);
                                bb_view := array2d_view(bb, ldb, row_major=false);
                                cc_view := array2d_view(cc, ldc, row_major=false);
                                cs_view := array2d_view(cs, ldc, row_major=false);

                                cgemm(transa, transb, m, n, k, alpha, aa_view, lda, bb_view, ldb, beta, *cc_view, ldc);

                                isame1  := transa == transas;
                                isame2  := transb == transbs;
                                isame3  := ms == m;
                                isame4  := ns == n;
                                isame5  := ks == k;
                                isame6  := als == alpha;
                                isame7  := lce(as.*, aa.*, laa);
                                isame8  := ldas == lda;
                                isame9  := lce(bs.*, bb.*, lbb);
                                isame10 := ldbs == ldb;
                                isame11 := bls == beta;
                                isame12 := ifx is_null then lce(cs.*, cc.*, lcc) else lceres("ge", #char " ", m, n, cs_view, cc_view, ldc);
                                isame13 := ldcs == ldc;

                                same := isame1 && isame2 && isame3 && isame4 && isame5 && isame6 &&
                                        isame7 && isame8 && isame9 && isame10 && isame11 && isame12 && isame13;
                                if !same {
                                    if !isame1  then err("cblat3 cchk1 | parameter 1 changed incorrectly");
                                    if !isame2  then err("cblat3 cchk1 | parameter 2 changed incorrectly");
                                    if !isame3  then err("cblat3 cchk1 | parameter 3 changed incorrectly");
                                    if !isame4  then err("cblat3 cchk1 | parameter 4 changed incorrectly");
                                    if !isame5  then err("cblat3 cchk1 | parameter 5 changed incorrectly");
                                    if !isame6  then err("cblat3 cchk1 | parameter 6 changed incorrectly");
                                    if !isame7  then err("cblat3 cchk1 | parameter 7 changed incorrectly");
                                    if !isame8  then err("cblat3 cchk1 | parameter 8 changed incorrectly");
                                    if !isame9  then err("cblat3 cchk1 | parameter 9 changed incorrectly");
                                    if !isame10 then err("cblat3 cchk1 | parameter 10 changed incorrectly");
                                    if !isame11 then err("cblat3 cchk1 | parameter 11 changed incorrectly");
                                    if !isame12 then err("cblat3 cchk1 | parameter 12 changed incorrectly");
                                    if !isame13 then err("cblat3 cchk1 | parameter 13 changed incorrectly");
                                    err("cblat3 cchk1 | sname = %", sname);
                                    err("cblat3 cchk1 | transa = %, transb = %, m = %, n = %, k = %, alpha = %, lda = %, ldb = %, beta = %, ldc = %",
                                        char_to_string(transa), char_to_string(transb), m, n, k, alpha, lda, ldb, beta, ldc,
                                        exit_program=fatal);
                                    return false;
                                }

                                if !is_null {
                                    errv := zero(U);
                                    ok := cblat3_cmmch(transa, transb, m, n, k, alpha, a.*, nmax, b.*, nmax, beta, c.*, nmax, ct, g, cc_view, ldc, eps, *errv, fatal, true);
                                    pass &= ok;
                                    errmax = max(errmax, errv);
                                    if !ok then return false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if errmax < thresh {
        println("% passed the computational tests (% calls)", sname, nc, color=.FG_GREEN);
    } else {
        err("% completed the computational tests (% calls) but with maximum test ratio % - suspect", sname, nc, errmax);
        pass = false;
    }

    return pass;
}


//  ######   ######  ##     ## ##    ##    #######  
// ##    ## ##    ## ##     ## ##   ##    ##     ## 
// ##       ##       ##     ## ##  ##            ## 
// ##       ##       ######### #####       #######  
// ##       ##       ##     ## ##  ##     ##        
// ##    ## ##    ## ##     ## ##   ##    ##        
//  ######   ######  ##     ## ##    ##   ######### 

cblat3_cchk2 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nalf: int,
    alf: []$T,
    nbet: int,
    bet: []T,
    nmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    b: *$B/ArrayType2D,
    bb: *$BB/ArrayType1D,
    bs: *$BS/ArrayType1D,
    c: *$C/ArrayType2D,
    cc: *$CC/ArrayType1D,
    cs: *$CS/ArrayType1D,
    ct: *$CT/ArrayType1D,
    g: *$G/ArrayType1D
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == B.data_type);
        assert(T == BB.data_type);
        assert(T == BS.data_type);
        assert(T == C.data_type);
        assert(T == CC.data_type);
        assert(T == CS.data_type);
        assert(T == CT.data_type);
        assert(U == G.data_type);
    }

    println_push("cchk2");

    ichs :: "LR";
    ichu :: "UL";

    pass := true;
    reset := true;
    errmax := zero(U);
    czero := zero(T);
    mat_type := String.slice(sname, 1, 2);
    conj := mat_type == "HE";

    nc : int = 0;

    for im:1..nidim {
        m := idim[im-1];

        for in:1..nidim {
            n := idim[in-1];

            ldc := m;
            if ldc < nmax then ldc += 1;
            if ldc > nmax then continue;
            lcc := ldc*n;
            is_null := n <= 0 || m <= 0;

            ldb := m;
            if ldb < nmax then ldb += 1;
            if ldb > nmax then continue;
            lbb := ldb*n;

            cmake("ge", 0, 0, m, n, b, nmax, bb, ldb, m-1, n-1, *reset, czero);

            for ics:1..2 {
                side := ichs[ics-1];
                left := char_eq(side, #char "l");

                na := ifx left then m else n;

                lda := na;
                if lda < nmax then lda += 1;
                if lda > nmax then continue;
                laa := lda*na;

                for icu:1..2 {
                    uplo := ichu[icu-1];

                    cmake(mat_type, uplo, 0, na, na, a, nmax, aa, lda, na-1, na-1, *reset, czero);

                    for ia:1..nalf {
                        alpha := alf[ia-1];

                        for ib:1..nbet {
                            beta := bet[ib-1];

                            cmake("ge", 0, 0, m, n, c, nmax, cc, ldc, m-1, n-1, *reset, czero);

                            nc += 1;

                            sides := side;
                            uplos := uplo;
                            ms := m;
                            ns := n;
                            als := alpha;
                            for i:1..laa {
                                fset(as, i, fget(aa.*, i));
                            }
                            ldas := lda;
                            for i:1..lbb {
                                fset(bs, i, fget(bb.*, i));
                            }
                            ldbs := ldb;
                            bls := beta;
                            for i:1..lcc {
                                fset(cs, i, fget(cc.*, i));
                            }
                            ldcs := ldc;

                            aa_view := array2d_view(aa, lda, row_major=false);
                            bb_view := array2d_view(bb, ldb, row_major=false);
                            cc_view := array2d_view(cc, ldc, row_major=false);
                            cs_view := array2d_view(cs, ldc, row_major=false);

                            if conj {
                                chemm(side, uplo, m, n, alpha, aa_view, lda, *bb_view, ldb, beta, *cc_view, ldc);
                            } else {
                                csymm(side, uplo, m, n, alpha, aa_view, lda, *bb_view, ldb, beta, *cc_view, ldc);
                            }

                            isame1  := side == sides;
                            isame2  := uplo == uplos;
                            isame3  := ms == m;
                            isame4  := ns == n;
                            isame5  := als == alpha;
                            isame6  := lce(as.*, aa.*, laa);
                            isame7  := ldas == lda;
                            isame8  := lce(bs.*, bb.*, lbb);
                            isame9  := ldbs == ldb;
                            isame10 := bls == beta;
                            isame11 := ifx is_null then lce(cs.*, cc.*, lcc) else lceres("ge", #char " ", m, n, cs_view, cc_view, ldc);
                            isame12 := ldcs == ldc;

                            same := isame1 && isame2 && isame3 && isame4 && isame5 && isame6 &&
                                    isame7 && isame8 && isame9 && isame10 && isame11 && isame12;
                            if !same {
                                if !isame1  then err("cblat3 cchk2 | parameter 1 changed incorrectly");
                                if !isame2  then err("cblat3 cchk2 | parameter 2 changed incorrectly");
                                if !isame3  then err("cblat3 cchk2 | parameter 3 changed incorrectly");
                                if !isame4  then err("cblat3 cchk2 | parameter 4 changed incorrectly");
                                if !isame5  then err("cblat3 cchk2 | parameter 5 changed incorrectly");
                                if !isame6  then err("cblat3 cchk2 | parameter 6 changed incorrectly");
                                if !isame7  then err("cblat3 cchk2 | parameter 7 changed incorrectly");
                                if !isame8  then err("cblat3 cchk2 | parameter 8 changed incorrectly");
                                if !isame9  then err("cblat3 cchk2 | parameter 9 changed incorrectly");
                                if !isame10 then err("cblat3 cchk2 | parameter 10 changed incorrectly");
                                if !isame11 then err("cblat3 cchk2 | parameter 11 changed incorrectly");
                                if !isame12 then err("cblat3 cchk2 | parameter 12 changed incorrectly");
                                err("cblat3 cchk2 | sname = %", sname);
                                err("cblat3 cchk2 | side = %, uplo = %, m = %, n = %, alpha = %, lda = %, ldb = %, beta = %, ldc = %",
                                    char_to_string(side), char_to_string(uplo), m, n, alpha, lda, ldb, beta, ldc,
                                    exit_program=fatal);
                                return false;
                            }

                            if !is_null {
                                errv := zero(U);
                                ok := false;
                                if left {
                                    ok = cblat3_cmmch(#char "n", #char "n", m, n, m, alpha, a.*, nmax, b.*, nmax, beta, c.*, nmax, ct, g, cc_view, ldc, eps, *errv, fatal, true);
                                } else {
                                    ok = cblat3_cmmch(#char "n", #char "n", m, n, n, alpha, b.*, nmax, a.*, nmax, beta, c.*, nmax, ct, g, cc_view, ldc, eps, *errv, fatal, true);
                                }
                                pass &= ok;
                                errmax = max(errmax, errv);
                                if !ok then return false;
                            }
                        }
                    }
                }
            }
        }
    }

    if errmax < thresh {
        println("% passed the computational tests (% calls)", sname, nc, color=.FG_GREEN);
    } else {
        err("% completed the computational tests (% calls) but with maximum test ratio % - suspect", sname, nc, errmax);
        pass = false;
    }

    return pass;
}


//  ######   ######  ##     ## ##    ##    #######  
// ##    ## ##    ## ##     ## ##   ##    ##     ## 
// ##       ##       ##     ## ##  ##            ## 
// ##       ##       ######### #####       #######  
// ##       ##       ##     ## ##  ##            ## 
// ##    ## ##    ## ##     ## ##   ##    ##     ## 
//  ######   ######  ##     ## ##    ##    #######  

cblat3_cchk3 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nalf: int,
    alf: []$T,
    nmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    b: *$B/ArrayType2D,
    bb: *$BB/ArrayType1D,
    bs: *$BS/ArrayType1D,
    ct: *$CT/ArrayType1D,
    g: *$G/ArrayType1D,
    c: *$C/ArrayType2D
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == B.data_type);
        assert(T == BB.data_type);
        assert(T == BS.data_type);
        assert(T == CT.data_type);
        assert(U == G.data_type);
        assert(T == C.data_type);
    }

    println_push("cchk3");

    ichu :: "UL";
    icht :: "NTC";
    ichd :: "UN";
    ichs :: "LR";

    pass := true;
    reset := true;
    errmax := zero(U);
    czero := zero(T);
    cone  := one(T);
    is_mm := String.slice(sname, 3, 5) == "MM";
    is_sm := String.slice(sname, 3, 5) == "SM";

    for j:1..nmax {
        for i:1..nmax {
            fset(c, i, j, czero);
        }
    }

    nc : int = 0;

    for im:1..nidim {
        m := idim[im-1];

        for in:1..nidim {
            n := idim[in-1];

            ldb := m;
            if ldb < nmax then ldb += 1;
            if ldb > nmax then continue;
            lbb := ldb*n;
            is_null := m <= 0 || n <= 0;

            for ics:1..2 {
                side := ichs[ics-1];
                left := char_eq(side, #char "l");
                na := ifx left then m else n;

                lda := na;
                if lda < nmax then lda += 1;
                if lda > nmax then continue;
                laa := lda*na;

                for icu:1..2 {
                    uplo := ichu[icu-1];

                    for ict:1..3 {
                        transa := icht[ict-1];

                        for icd:1..2 {
                            diag := ichd[icd-1];

                            for ia:1..nalf {
                                alpha := alf[ia-1];

                                cmake("tr", uplo, diag, na, na, a, nmax, aa, lda, na-1, na-1, *reset, czero);
                                cmake("ge", 0, 0, m, n, b, nmax, bb, ldb, m-1, n-1, *reset, czero);

                                nc += 1;

                                sides := side;
                                uplos := uplo;
                                transas := transa;
                                diags := diag;
                                ms := m;
                                ns := n;
                                als := alpha;
                                for i:1..laa {
                                    fset(as, i, fget(aa.*, i));
                                }
                                ldas := lda;
                                for i:1..lbb {
                                    fset(bs, i, fget(bb.*, i));
                                }
                                ldbs := ldb;

                                aa_view := array2d_view(aa, lda, row_major=false);
                                bb_view := array2d_view(bb, ldb, row_major=false);
                                bs_view := array2d_view(bs, ldb, row_major=false);

                                if is_mm {
                                    ctrmm(side, uplo, transa, diag, m, n, alpha, aa_view, lda, *bb_view, ldb);
                                } else if is_sm {
                                    ctrsm(side, uplo, transa, diag, m, n, alpha, aa_view, lda, *bb_view, ldb);
                                }

                                isame1  := side == sides;
                                isame2  := uplo == uplos;
                                isame3  := transa == transas;
                                isame4  := diag == diags;
                                isame5  := ms == m;
                                isame6  := ns == n;
                                isame7  := als == alpha;
                                isame8  := lce(as.*, aa.*, laa);
                                isame9  := ldas == lda;
                                isame10 := ifx is_null then lce(bs.*, bb.*, lbb) else lceres("ge", #char " ", m, n, bs_view, bb_view, ldb);
                                isame11 := ldbs == ldb;

                                same := isame1 && isame2 && isame3 && isame4 && isame5 && isame6 &&
                                        isame7 && isame8 && isame9 && isame10 && isame11;
                                if !same {
                                    if !isame1  then err("cblat3 cchk3 | parameter 1 changed incorrectly");
                                    if !isame2  then err("cblat3 cchk3 | parameter 2 changed incorrectly");
                                    if !isame3  then err("cblat3 cchk3 | parameter 3 changed incorrectly");
                                    if !isame4  then err("cblat3 cchk3 | parameter 4 changed incorrectly");
                                    if !isame5  then err("cblat3 cchk3 | parameter 5 changed incorrectly");
                                    if !isame6  then err("cblat3 cchk3 | parameter 6 changed incorrectly");
                                    if !isame7  then err("cblat3 cchk3 | parameter 7 changed incorrectly");
                                    if !isame8  then err("cblat3 cchk3 | parameter 8 changed incorrectly");
                                    if !isame9  then err("cblat3 cchk3 | parameter 9 changed incorrectly");
                                    if !isame10 then err("cblat3 cchk3 | parameter 10 changed incorrectly");
                                    if !isame11 then err("cblat3 cchk3 | parameter 11 changed incorrectly");
                                    err("cblat3 cchk3 | sname = %", sname);
                                    err("cblat3 cchk3 | side = %, uplo = %, transa = %, diag = %, m = %, n = %, alpha = %, lda = %, ldb = %",
                                        char_to_string(side), char_to_string(uplo), char_to_string(transa), char_to_string(diag), m, n, alpha, lda, ldb,
                                        exit_program=fatal);
                                    return false;
                                }

                                if !is_null {
                                    errv := zero(U);
                                    ok := false;
                                    if is_mm {
                                        if left {
                                            ok = cblat3_cmmch(transa, #char "n", m, n, m, alpha, a.*, nmax, b.*, nmax, czero, c.*, nmax, ct, g, bb_view, ldb, eps, *errv, fatal, true);
                                        } else {
                                            ok = cblat3_cmmch(#char "n", transa, m, n, n, alpha, b.*, nmax, a.*, nmax, czero, c.*, nmax, ct, g, bb_view, ldb, eps, *errv, fatal, true);
                                        }
                                    } else if is_sm {
                                        for j:1..n {
                                            for i:1..m {
                                                fset(c, i, j, fget(bb_view, i, j));
                                                fset(bb, i + (j-1)*ldb, alpha*fget(b.*, i, j));
                                            }
                                        }

                                        if left {
                                            ok = cblat3_cmmch(transa, #char "n", m, n, m, cone, a.*, nmax, c.*, nmax, czero, b.*, nmax, ct, g, bb_view, ldb, eps, *errv, fatal, false);
                                        } else {
                                            ok = cblat3_cmmch(#char "n", transa, m, n, n, cone, c.*, nmax, a.*, nmax, czero, b.*, nmax, ct, g, bb_view, ldb, eps, *errv, fatal, false);
                                        }
                                    }

                                    pass &= ok;
                                    errmax = max(errmax, errv);
                                    if !ok then return false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if errmax < thresh {
        println("% passed the computational tests (% calls)", sname, nc, color=.FG_GREEN);
    } else {
        err("% completed the computational tests (% calls) but with maximum test ratio % - suspect", sname, nc, errmax);
        pass = false;
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

cblat3_cchk4 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nalf: int,
    alf: []$T,
    nbet: int,
    bet: []T,
    nmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    b: *$B/ArrayType2D,
    bb: *$BB/ArrayType1D,
    bs: *$BS/ArrayType1D,
    c: *$C/ArrayType2D,
    cc: *$CC/ArrayType1D,
    cs: *$CS/ArrayType1D,
    ct: *$CT/ArrayType1D,
    g: *$G/ArrayType1D
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == B.data_type);
        assert(T == BB.data_type);
        assert(T == BS.data_type);
        assert(T == C.data_type);
        assert(T == CC.data_type);
        assert(T == CS.data_type);
        assert(T == CT.data_type);
        assert(U == G.data_type);
    }

    println_push("cchk4");

    icht :: "NC";
    ichu :: "UL";

    pass := true;
    reset := true;
    errmax := zero(U);
    czero := zero(T);
    rzero := zero(U);
    rone  := one(U);
    conj := String.slice(sname, 1, 2) == "HE";

    nc : int = 0;

    for in:1..nidim {
        n := idim[in-1];

        ldc := n;
        if ldc < nmax then ldc += 1;
        if ldc > nmax then continue;
        lcc := ldc*n;

        for ik:1..nidim {
            k := idim[ik-1];

            for ict:1..2 {
                trans := icht[ict-1];
                tran := char_eq(trans, #char "c");
                if tran && !conj then trans = #char "t";

                ma, na : int = ---;
                if tran {
                    ma = k;
                    na = n;
                } else {
                    ma = n;
                    na = k;
                }

                lda := ma;
                if lda < nmax then lda += 1;
                if lda > nmax then continue;
                laa := lda*na;

                cmake("ge", 0, 0, ma, na, a, nmax, aa, lda, ma-1, na-1, *reset, czero);

                for icu:1..2 {
                    uplo := ichu[icu-1];
                    upper := char_eq(uplo, #char "u");

                    for ia:1..nalf {
                        alpha := alf[ia-1];
                        ralpha : U = ---;
                        if conj {
                            ralpha = real(alpha);
                            alpha = cpx(ralpha, rzero);
                        }

                        for ib:1..nbet {
                            beta := bet[ib-1];
                            rbeta : U = ---;
                            if conj {
                                rbeta = real(beta);
                                beta = cpx(rbeta, rzero);
                            }

                            is_null := n <= 0;
                            if conj {
                                is_null = is_null || (((k <= 0 || ralpha == rzero) && rbeta == rone));
                            }

                            cmake(String.slice(sname, 1, 2), uplo, 0, n, n, c, nmax, cc, ldc, n-1, n-1, *reset, czero);

                            nc += 1;

                            uplos := uplo;
                            transs := trans;
                            ns := n;
                            ks := k;
                            rals : U = ---;
                            als  : T = ---;
                            if conj {
                                rals = ralpha;
                            } else {
                                als = alpha;
                            }
                            for i:1..laa {
                                fset(as, i, fget(aa.*, i));
                            }
                            ldas := lda;
                            rbets : U = ---;
                            bets  : T = ---;
                            if conj {
                                rbets = rbeta;
                            } else {
                                bets = beta;
                            }
                            for i:1..lcc {
                                fset(cs, i, fget(cc.*, i));
                            }
                            ldcs := ldc;

                            aa_view := array2d_view(aa, lda, row_major=false);
                            cc_view := array2d_view(cc, ldc, row_major=false);
                            cs_view := array2d_view(cs, ldc, row_major=false);

                            if conj {
                                cherk(uplo, trans, n, k, ralpha, aa_view, lda, rbeta, *cc_view, ldc);
                            } else {
                                csyrk(uplo, trans, n, k, alpha, aa_view, lda, beta, *cc_view, ldc);
                            }

                            isame1  := uplo == uplos;
                            isame2  := trans == transs;
                            isame3  := ns == n;
                            isame4  := ks == k;
                            isame5  := ifx conj then rals == ralpha else als == alpha;
                            isame6  := lce(as.*, aa.*, laa);
                            isame7  := ldas == lda;
                            isame8  := ifx conj then rbets == rbeta else bets == beta;
                            isame9  := ifx is_null then lce(cs.*, cc.*, lcc) else lceres(String.slice(sname, 1, 2), uplo, n, n, cs_view, cc_view, ldc);
                            isame10 := ldcs == ldc;

                            same := isame1 && isame2 && isame3 && isame4 && isame5 &&
                                    isame6 && isame7 && isame8 && isame9 && isame10;
                            if !same {
                                if !isame1  then err("cblat3 cchk4 | parameter 1 changed incorrectly");
                                if !isame2  then err("cblat3 cchk4 | parameter 2 changed incorrectly");
                                if !isame3  then err("cblat3 cchk4 | parameter 3 changed incorrectly");
                                if !isame4  then err("cblat3 cchk4 | parameter 4 changed incorrectly");
                                if !isame5  then err("cblat3 cchk4 | parameter 5 changed incorrectly");
                                if !isame6  then err("cblat3 cchk4 | parameter 6 changed incorrectly");
                                if !isame7  then err("cblat3 cchk4 | parameter 7 changed incorrectly");
                                if !isame8  then err("cblat3 cchk4 | parameter 8 changed incorrectly");
                                if !isame9  then err("cblat3 cchk4 | parameter 9 changed incorrectly");
                                if !isame10 then err("cblat3 cchk4 | parameter 10 changed incorrectly");
                                err("cblat3 cchk4 | sname = %", sname);
                                if conj {
                                    err("cblat3 cchk4 | uplo = %, trans = %, n = %, k = %, alpha = %, lda = %, beta = %, ldc = %",
                                        char_to_string(uplo), char_to_string(trans), n, k, ralpha, lda, rbeta, ldc,
                                        exit_program=fatal);
                                } else {
                                    err("cblat3 cchk4 | uplo = %, trans = %, n = %, k = %, alpha = %, lda = %, beta = %, ldc = %",
                                        char_to_string(uplo), char_to_string(trans), n, k, alpha, lda, beta, ldc,
                                        exit_program=fatal);
                                }
                                return false;
                            }

                            if !is_null {
                                transt : u8 = ---;
                                if conj {
                                    transt = #char "c";
                                } else {
                                    transt = #char "t";
                                }
                                jc := 1;
                                for j:1..n {
                                    jj, lj : int = ---;
                                    if upper {
                                        jj = 1;
                                        lj = j;
                                    } else {
                                        jj = j;
                                        lj = n - j + 1;
                                    }

                                    aw1 := array2d_view(as, nmax, row_major=false);
                                    aw2 := array2d_view(bs, nmax, row_major=false);
                                    cw1 := array2d_view(cs, nmax, row_major=false);
                                    cw2 := array2d_view(bb, nmax, row_major=false);

                                    if tran {
                                        for cj:1..lj {
                                            for ri:1..k {
                                                fset(*aw1, ri, cj, fget(a.*, ri, jj + cj - 1));
                                            }
                                        }
                                        for ri:1..k {
                                            fset(*aw2, ri, 1, fget(a.*, ri, j));
                                        }
                                    } else {
                                        for cj:1..k {
                                            for ri:1..lj {
                                                fset(*aw1, ri, cj, fget(a.*, jj + ri - 1, cj));
                                            }
                                        }
                                        for cj:1..k {
                                            fset(*aw2, 1, cj, fget(a.*, j, cj));
                                        }
                                    }

                                    for ri:1..lj {
                                        fset(*cw1, ri, 1, fget(c.*, jj + ri - 1, j));
                                        fset(*cw2, ri, 1, fget(cc_view, jj + ri - 1, j));
                                    }

                                    errv := zero(U);
                                    ok := false;
                                    if tran {
                                        ok = cblat3_cmmch(transt, #char "n", lj, 1, k, alpha, aw1, nmax, aw2, nmax, beta, cw1, nmax, ct, g, cw2, nmax, eps, *errv, fatal, true);
                                    } else {
                                        ok = cblat3_cmmch(#char "n", transt, lj, 1, k, alpha, aw1, nmax, aw2, nmax, beta, cw1, nmax, ct, g, cw2, nmax, eps, *errv, fatal, true);
                                    }

                                    pass &= ok;
                                    errmax = max(errmax, errv);
                                    if !ok {
                                        if n > 1 then err("cblat3 cchk4 | these are the results for column %", j);
                                        return false;
                                    }

                                    if upper {
                                        jc += ldc;
                                    } else {
                                        jc += ldc + 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if errmax < thresh {
        println("% passed the computational tests (% calls)", sname, nc, color=.FG_GREEN);
    } else {
        err("% completed the computational tests (% calls) but with maximum test ratio % - suspect", sname, nc, errmax);
        pass = false;
    }

    return pass;
}



//  ######   ######  ##     ## ##    ## ########
// ##    ## ##    ## ##     ## ##   ##  ##
// ##       ##       ##     ## ##  ##   ##
// ##       ##       ######### #####    #######
// ##       ##       ##     ## ##  ##         ##
// ##    ## ##    ## ##     ## ##   ##  ##    ##
//  ######   ######  ##     ## ##    ##  ######

cblat3_cchk5 :: (
    sname: string,
    eps: $U,
    thresh: U,
    $fatal: bool,
    nidim: int,
    idim: []int,
    nalf: int,
    alf: []$T,
    nbet: int,
    bet: []T,
    nmax: int,
    a: *$A/ArrayType2D,
    aa: *$AA/ArrayType1D,
    as: *$AS/ArrayType1D,
    b: *$B/ArrayType2D,
    bb: *$BB/ArrayType1D,
    bs: *$BS/ArrayType1D,
    c: *$C/ArrayType2D,
    cc: *$CC/ArrayType1D,
    cs: *$CS/ArrayType1D,
    ct: *$CT/ArrayType1D,
    g: *$G/ArrayType1D,
    w: *$W/ArrayType1D
) -> bool {
    #run {
        assert(U == get_number_inner_type(T));
        assert(T == A.data_type);
        assert(T == AA.data_type);
        assert(T == AS.data_type);
        assert(T == B.data_type);
        assert(T == BB.data_type);
        assert(T == BS.data_type);
        assert(T == C.data_type);
        assert(T == CC.data_type);
        assert(T == CS.data_type);
        assert(T == CT.data_type);
        assert(U == G.data_type);
        assert(T == W.data_type);
    }

    println_push("cchk5");

    icht :: "NC";
    ichu :: "UL";

    pass := true;
    reset := true;
    errmax := zero(U);
    czero := zero(T);
    cone  := one(T);
    rzero := zero(U);
    rone  := one(U);
    conj := String.slice(sname, 1, 2) == "HE";

    nc : int = 0;

    for in:1..nidim {
        n := idim[in-1];

        ldc := n;
        if ldc < nmax then ldc += 1;
        if ldc > nmax then continue;
        lcc := ldc*n;

        for ik:1..nidim {
            k := idim[ik-1];

            for ict:1..2 {
                trans := icht[ict-1];
                tran := char_eq(trans, #char "c");
                if tran && !conj then trans = #char "t";

                ma, na : int = ---;
                if tran {
                    ma = k;
                    na = n;
                } else {
                    ma = n;
                    na = k;
                }

                lda := ma;
                if lda < nmax then lda += 1;
                if lda > nmax then continue;
                laa := lda*na;

                cmake("ge", 0, 0, ma, na, a, nmax, aa, lda, ma-1, na-1, *reset, czero);

                ldb := lda;
                lbb := laa;
                cmake("ge", 0, 0, ma, na, b, nmax, bb, ldb, ma-1, na-1, *reset, czero);

                for icu:1..2 {
                    uplo := ichu[icu-1];
                    upper := char_eq(uplo, #char "u");

                    for ia:1..nalf {
                        alpha := alf[ia-1];

                        for ib:1..nbet {
                            beta := bet[ib-1];
                            rbeta : U = ---;
                            if conj {
                                rbeta = real(beta);
                                beta = cpx(rbeta, rzero);
                            }

                            is_null := n <= 0;
                            if conj {
                                is_null = is_null || (((k <= 0 || alpha == czero) && rbeta == rone));
                            }

                            cmake(String.slice(sname, 1, 2), uplo, 0, n, n, c, nmax, cc, ldc, n-1, n-1, *reset, czero);

                            nc += 1;

                            uplos := uplo;
                            transs := trans;
                            ns := n;
                            ks := k;
                            als := alpha;
                            for i:1..laa {
                                fset(as, i, fget(aa.*, i));
                            }
                            ldas := lda;
                            for i:1..lbb {
                                fset(bs, i, fget(bb.*, i));
                            }
                            ldbs := ldb;
                            rbets : U = ---;
                            bets  : T = ---;
                            if conj {
                                rbets = rbeta;
                            } else {
                                bets = beta;
                            }
                            for i:1..lcc {
                                fset(cs, i, fget(cc.*, i));
                            }
                            ldcs := ldc;

                            aa_view := array2d_view(aa, lda, row_major=false);
                            bb_view := array2d_view(bb, ldb, row_major=false);
                            cc_view := array2d_view(cc, ldc, row_major=false);
                            cs_view := array2d_view(cs, ldc, row_major=false);

                            if conj {
                                cher2k(uplo, trans, n, k, alpha, aa_view, lda, bb_view, ldb, rbeta, *cc_view, ldc);
                            } else {
                                csyr2k(uplo, trans, n, k, alpha, aa_view, lda, bb_view, ldb, beta, *cc_view, ldc);
                            }

                            isame1  := uplo == uplos;
                            isame2  := trans == transs;
                            isame3  := ns == n;
                            isame4  := ks == k;
                            isame5  := als == alpha;
                            isame6  := lce(as.*, aa.*, laa);
                            isame7  := ldas == lda;
                            isame8  := lce(bs.*, bb.*, lbb);
                            isame9  := ldbs == ldb;
                            isame10 := ifx conj then rbets == rbeta else bets == beta;
                            isame11 := ifx is_null then lce(cs.*, cc.*, lcc) else lceres("he", uplo, n, n, cs_view, cc_view, ldc);
                            isame12 := ldcs == ldc;

                            same := isame1 && isame2 && isame3 && isame4 && isame5 && isame6 &&
                                    isame7 && isame8 && isame9 && isame10 && isame11 && isame12;
                            if !same {
                                if !isame1  then err("cblat3 cchk5 | parameter 1 changed incorrectly");
                                if !isame2  then err("cblat3 cchk5 | parameter 2 changed incorrectly");
                                if !isame3  then err("cblat3 cchk5 | parameter 3 changed incorrectly");
                                if !isame4  then err("cblat3 cchk5 | parameter 4 changed incorrectly");
                                if !isame5  then err("cblat3 cchk5 | parameter 5 changed incorrectly");
                                if !isame6  then err("cblat3 cchk5 | parameter 6 changed incorrectly");
                                if !isame7  then err("cblat3 cchk5 | parameter 7 changed incorrectly");
                                if !isame8  then err("cblat3 cchk5 | parameter 8 changed incorrectly");
                                if !isame9  then err("cblat3 cchk5 | parameter 9 changed incorrectly");
                                if !isame10 then err("cblat3 cchk5 | parameter 10 changed incorrectly");
                                if !isame11 then err("cblat3 cchk5 | parameter 11 changed incorrectly");
                                if !isame12 then err("cblat3 cchk5 | parameter 12 changed incorrectly");
                                err("cblat3 cchk5 | sname = %", sname);
                                if conj {
                                    err("cblat3 cchk5 | uplo = %, trans = %, n = %, k = %, alpha = %, lda = %, ldb = %, beta = %, ldc = %",
                                        char_to_string(uplo), char_to_string(trans), n, k, alpha, lda, ldb, rbeta, ldc,
                                        exit_program=fatal);
                                } else {
                                    err("cblat3 cchk5 | uplo = %, trans = %, n = %, k = %, alpha = %, lda = %, ldb = %, beta = %, ldc = %",
                                        char_to_string(uplo), char_to_string(trans), n, k, alpha, lda, ldb, beta, ldc,
                                        exit_program=fatal);
                                }
                                return false;
                            }

                            if !is_null {
                                transt : u8 = ---;
                                if conj {
                                    transt = #char "c";
                                } else {
                                    transt = #char "t";
                                }

                                work_a := array2d_view(as, ifx tran then 2*nmax else nmax, row_major=false);
                                work_c := array2d_view(cs, nmax, row_major=false);
                                work_r := array2d_view(bs, nmax, row_major=false);
                                work_w := array2d_view(w, 2*nmax, row_major=false);

                                jc := 1;
                                for j:1..n {
                                    jj, lj : int = ---;
                                    if upper {
                                        jj = 1;
                                        lj = j;
                                    } else {
                                        jj = j;
                                        lj = n - j + 1;
                                    }

                                    for ri:1..lj {
                                        fset(*work_c, ri, 1, fget(c.*, jj + ri - 1, j));
                                        fset(*work_r, ri, 1, fget(cc_view, jj + ri - 1, j));
                                    }

                                    errv := zero(U);
                                    ok := false;
                                    if tran {
                                        for cj:1..lj {
                                            for ri:1..k {
                                                fset(*work_a, ri, cj, fget(a.*, ri, jj + cj - 1));
                                                fset(*work_a, k + ri, cj, fget(b.*, ri, jj + cj - 1));
                                            }
                                        }
                                        for ri:1..k {
                                            fset(*work_w, ri, 1, alpha*fget(b.*, ri, j));
                                            if conj {
                                                fset(*work_w, k + ri, 1, conjugate(alpha)*fget(a.*, ri, j));
                                            } else {
                                                fset(*work_w, k + ri, 1, alpha*fget(a.*, ri, j));
                                            }
                                        }

                                        ok = cblat3_cmmch(transt, #char "n", lj, 1, 2*k, cone, work_a, 2*nmax, work_w, 2*nmax, beta, work_c, nmax, ct, g, work_r, nmax, eps, *errv, fatal, true);
                                    } else {
                                        for cj:1..k {
                                            for ri:1..lj {
                                                fset(*work_a, ri, cj, fget(a.*, jj + ri - 1, cj));
                                                fset(*work_a, ri, k + cj, fget(b.*, jj + ri - 1, cj));
                                            }
                                        }
                                        for ri:1..k {
                                            if conj {
                                                fset(*work_w, ri, 1, alpha*conjugate(fget(b.*, j, ri)));
                                                fset(*work_w, k + ri, 1, conjugate(alpha*fget(a.*, j, ri)));
                                            } else {
                                                fset(*work_w, ri, 1, alpha*fget(b.*, j, ri));
                                                fset(*work_w, k + ri, 1, alpha*fget(a.*, j, ri));
                                            }
                                        }

                                        ok = cblat3_cmmch(#char "n", #char "n", lj, 1, 2*k, cone, work_a, nmax, work_w, 2*nmax, beta, work_c, nmax, ct, g, work_r, nmax, eps, *errv, fatal, true);
                                    }

                                    pass &= ok;
                                    errmax = max(errmax, errv);
                                    if !ok {
                                        if n > 1 then err("cblat3 cchk5 | these are the results for column %", j);
                                        return false;
                                    }

                                    if upper {
                                        jc += ldc;
                                    } else {
                                        jc += ldc + 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if errmax < thresh {
        println("% passed the computational tests (% calls)", sname, nc, color=.FG_GREEN);
    } else {
        err("% completed the computational tests (% calls) but with maximum test ratio % - suspect", sname, nc, errmax);
        pass = false;
    }

    return pass;
}

