// ? LAPACK cblat2 helper routines
// ? Loaded and driven by _blat2.jai


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

    #if DEBUG println_push("cmvch");

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
    DEBUG :: false;
    #if DEBUG println_push("cchk1");

    pass := true;

    ich :: "ntc";

    czero := zero(T);
    cone  := one(T);
    chalf := ncast(T, 0.5);
    rzero := zero(U);
    rone  := one(U);

    full   := char_eq(sname[2], #char "e", #char "y");
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
                                    // println("in %, im %, iku %, ic %, ix %, iy %, ia %, ib %", in, im, iku, ic, ix, iy, ia, ib);
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
                                        #if T == float32 {
                                            // println("sgemv");
                                            sgemv(trans, m, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == float64 {
                                            // println("dgemv");
                                            dgemv(trans, m, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C32 {
                                            // println("cgemv");
                                            cgemv(trans, m, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C64 {
                                            // println("zgemv");
                                            zgemv(trans, m, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                    } else if banded {
                                        #if T == float32 {
                                            // println("sgbmv");
                                            sgbmv(trans, m, n, kl, ku, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == float64 {
                                            // println("dgbmv");
                                            dgbmv(trans, m, n, kl, ku, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C32 {
                                            // println("cgbmv");
                                            cgbmv(trans, m, n, kl, ku, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C64 {
                                            // println("zgbmv");
                                            zgbmv(trans, m, n, kl, ku, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
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
    DEBUG :: false;
    #if DEBUG println_push("cchk2");

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
                                // println("in %, iku %, iu %, ix %, iy %, ia %, ib %", in, iku, iu, ix, iy, ia, ib);
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
                                    #if T == float32 {
                                        // println("sspmv");
                                        sspmv(uplo, n, alpha, aa, _xx, incx, beta, yy, incy);
                                    }
                                    #if T == float64 {
                                        // println("dspmv");
                                        dspmv(uplo, n, alpha, aa, _xx, incx, beta, yy, incy);
                                    }
                                    #if T == C32 {
                                        // println("chpmv");
                                        chpmv(uplo, n, alpha, aa, _xx, incx, beta, yy, incy);
                                    }
                                    #if T == C64 {
                                        // println("zhpmv");
                                        zhpmv(uplo, n, alpha, aa, _xx, incx, beta, yy, incy);
                                    }
                                } else {
                                    aa_view := array2d_view(aa, lda, row_major=false);
                                    if banded {
                                        #if T == float32 {
                                            // println("ssbmv");
                                            ssbmv(uplo, n, k, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == float64 {
                                            // println("dsbmv");
                                            dsbmv(uplo, n, k, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C32 {
                                            // println("chbmv");
                                            chbmv(uplo, n, k, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C64 {
                                            // println("zhbmv");
                                            zhbmv(uplo, n, k, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                    } else {
                                        #if T == float32 {
                                            // println("ssymv");
                                            ssymv(uplo, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == float64 {
                                            // println("dsymv");
                                            dsymv(uplo, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C32 {
                                            // println("chemv");
                                            chemv(uplo, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
                                        #if T == C64 {
                                            // println("zhemv");
                                            zhemv(uplo, n, alpha, aa_view, lda, _xx, incx, beta, yy, incy);
                                        }
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
    DEBUG :: false;
    #if DEBUG println_push("cchk3");

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
    DEBUG :: false;
    #if DEBUG println_push("cchk4");

    pass := true;

    czero := zero(T);
    cone  := one(T);
    chalf := ncast(T, 0.5);
    rzero := zero(U);

    conjug := sname.count > 4 && char_eq(sname[4], #char "c");
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
                        #if T == float32 {
                            sger(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                        }
                        #if T == float64 {
                            dger(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                        }
                        #if T == C32 {
                            if conjug {
                                cgerc(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            } else {
                                cgeru(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            }
                        }
                        #if T == C64 {
                            if conjug {
                                zgerc(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            } else {
                                zgeru(m, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            }
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
    DEBUG :: false;
    #if DEBUG println_push("cchk5");

    pass := true;

    ichu :: "ul";

    czero := zero(T);
    cone  := one(T);
    chalf := ncast(T, 0.5);
    rzero := zero(U);
    rone  := one(U);

    full   := char_eq(sname[2], #char "e", #char "y");
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
                    ralpha : U = ---;
                    alpha  : T = ---;
                    #if T == C32 || T == C64 {
                        ralpha = real(alf[ia]);
                        if imag(alf[ia]) != rzero then continue;
                        alpha = cpx(ralpha, rzero);
                    } else {
                        ralpha = alf[ia];
                        alpha = alf[ia];
                    }
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
                        #if T == float32 {
                            ssyr(uplo, n, ralpha, _xx, incx, *aa_view, lda);
                        }
                        #if T == float64 {
                            dsyr(uplo, n, ralpha, _xx, incx, *aa_view, lda);
                        }
                        #if T == C32 {
                            cher(uplo, n, ralpha, _xx, incx, *aa_view, lda);
                        }
                        #if T == C64 {
                            zher(uplo, n, ralpha, _xx, incx, *aa_view, lda);
                        }
                    } else if packed {
                        #if T == float32 {
                            sspr(uplo, n, ralpha, _xx, incx, aa);
                        }
                        #if T == float64 {
                            dspr(uplo, n, ralpha, _xx, incx, aa);
                        }
                        #if T == C32 {
                            chpr(uplo, n, ralpha, _xx, incx, aa);
                        }
                        #if T == C64 {
                            zhpr(uplo, n, ralpha, _xx, incx, aa);
                        }
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
                        ok = lceres(String.slice(sname, 1, 2), uplo, n, n, aa_view, as_view, lda);
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


//  ######   ######  ##     ## ##    ##  #######
// ##    ## ##    ## ##     ## ##   ##  ##     ##
// ##       ##       ##     ## ##  ##   ##
// ##       ##       ######### #####    ########
// ##       ##       ##     ## ##  ##   ##     ##
// ##    ## ##    ## ##     ## ##   ##  ##     ##
//  ######   ######  ##     ## ##    ##  #######

cchk6 :: (
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
    DEBUG :: false;
    #if DEBUG println_push("cchk6");

    pass := true;

    ichu :: "ul";

    czero := zero(T);
    cone  := one(T);
    chalf := ncast(T, 0.5);
    rzero := zero(U);

    full   := char_eq(sname[2], #char "e", #char "y");
    packed := char_eq(sname[2], #char "p");

    nargs := ifx full then 9 else 8;

    nc := 0;
    reset := true;
    errmax := rzero;

    w : Array1D(2, T);

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
                        _null := n <= 0 || alpha == czero;

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
                        for i:0..ly-1 {
                            set(ys, i, get(yy, i));
                        }

                        if full {
                            aa_view := array2d_view(aa, lda, row_major=false);
                            #if T == float32 {
                                ssyr2(uplo, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            }
                            #if T == float64 {
                                dsyr2(uplo, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            }
                            #if T == C32 {
                                cher2(uplo, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            }
                            #if T == C64 {
                                zher2(uplo, n, alpha, _xx, incx, yy, incy, *aa_view, lda);
                            }
                        } else if packed {
                            #if T == float32 {
                                sspr2(uplo, n, alpha, _xx, incx, yy, incy, aa);
                            }
                            #if T == float64 {
                                dspr2(uplo, n, alpha, _xx, incx, yy, incy, aa);
                            }
                            #if T == C32 {
                                chpr2(uplo, n, alpha, _xx, incx, yy, incy, aa);
                            }
                            #if T == C64 {
                                zhpr2(uplo, n, alpha, _xx, incx, yy, incy, aa);
                            }
                        }

                        ok := lce(xs, _xx, lx);
                        pass &= ok;
                        if !ok {
                            err(" lce | xs = %", xs);
                            err(" lce | xx = %", _xx);
                            err(" lce | lx = %", lx, exit_program=fatal);
                        }

                        ok = lce(ys, yy, ly);
                        pass &= ok;
                        if !ok {
                            err(" lce | ys = %", ys);
                            err(" lce | yy = %", yy);
                            err(" lce | ly = %", ly, exit_program=fatal);
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
                            ok = lceres(String.slice(sname, 1, 2), uplo, n, n, aa_view, as_view, lda);
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
                                    fset(x, i, fget(_xx, 1 + (i-1)*abs(incx)));
                                }
                            } else {
                                for i:1..n {
                                    fset(x, i, fget(_xx, 1 + (n-i)*abs(incx)));
                                }
                            }

                            if incy > 0 {
                                for i:1..n {
                                    fset(y, i, fget(yy, 1 + (i-1)*abs(incy)));
                                }
                            } else {
                                for i:1..n {
                                    fset(y, i, fget(yy, 1 + (n-i)*abs(incy)));
                                }
                            }

                            z_view := array2d_view(z, nmax, row_major=false);
                            ja := 1;
                            for j:1..n {
                                fset(*w, 1, alpha*conjugate(fget(y, j)));
                                fset(*w, 2, conjugate(alpha*fget(x, j)));

                                jj, lj : int = ---;
                                if upper {
                                    jj = 1;
                                    lj = j;
                                } else {
                                    jj = j;
                                    lj = n-j+1;
                                }

                                for i:1..lj {
                                    fset(z, i, fget(x, jj+i-1));
                                    fset(z, i + nmax, fget(y, jj+i-1));
                                    fset(_xx, i, fget(a.*, jj+i-1, j));
                                    fset(ys, i, fget(aa, ja+i-1));
                                }

                                ok = cmvch(#char "n", lj, 2, cone, z_view, nmax, w, 1, cone, _xx, 1, yt, g, ys, eps, fatal);
                                pass &= ok;

                                erri := rzero;
                                for i:1..lj {
                                    gi := fget(g.*, i);
                                    diff := abs(fget(yt.*, i) - fget(ys, i))/eps;
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

    rogue  : T = ---;
    #if T == C32 || T == C64 {
        rogue = cpx(ncast(U, -1.0e10), ncast(U, 1.0e10));
    } else {
        rogue = ncast(T, ncast(U, -1.0e10));
    }
    rrogue := ncast(U, -1.0e10);

    gen := char_eq(type[0], #char "g");
    her := char_eq(type[0], #char "h");
    sym := her || char_eq(type[0], #char "s");
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

            if her {
                fset(a, j, i, conjugate(fget(a.*,i,j)));
            } else if sym {
                fset(a, j, i, fget(a.*,i,j));
            } else {
                fset(a, j, i, czero);
            }
        }
        if her  then fset(a, j, j, real_nc(fget(a.*,j,j)));
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
    } else if String.equal_nocase(type, "he") || String.equal_nocase(type, "sy") || String.equal_nocase(type, "tr") {
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
            if her {
                #if T == C32 || T == C64 {
                    jj := j + (j-1)*lda;
                    fset(aa,jj, cpx(real(fget(aa.*,jj)), rrogue));
                }
            }
        }
    } else if String.equal_nocase(type, "hb") || String.equal_nocase(type, "sb") || String.equal_nocase(type, "tb") {
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
            if her {
                #if T == C32 || T == C64 {
                    jj := kk + (j-1)*lda;
                    fset(aa,jj, cpx(real(fget(aa.*,jj)), rrogue));
                }
            }
        }
    } else if String.equal_nocase(type, "hp") || String.equal_nocase(type, "sp") || String.equal_nocase(type, "tp") {
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
                    if her {
                        #if T == C32 || T == C64 {
                            fset(aa,ioff, cpx(real(fget(aa.*,ioff)), rrogue));
                        } else {
                            fset(aa,ioff, fget(a.*,i,j));
                        }
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


// Generates test values uniformly distributed between -0.5 and 0.5.
// Complex types use independent real and imaginary parts.

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
    #if T == C32 || T == C64 {
        return cpx(ncast(X, cast(float64)(i-500)/1001.0), ncast(X, cast(float64)(j-500)/1001.0));
    } else {
        return ncast(T, cast(float64)(i-500)/1001.0);
    }
}
