// ? LAPACK *gsvj1()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgsvj1/cgsvj1/dgsvj1/zgsvj1
// *GSVJ1 is a Jacobi-SVD preprocessing kernel targeting selected pivots.

sgsvj1 :: #bake_arguments _gsvj1_real(RT=float32);
dgsvj1 :: #bake_arguments _gsvj1_real(RT=float64);
cgsvj1 :: #bake_arguments _gsvj1_complex(RT=float32, MT=C32);
zgsvj1 :: #bake_arguments _gsvj1_complex(RT=float64, MT=C64);

#scope_file
_gsvj1_real :: (
    $$jobv: u8,
    $$m: int,
    $$n: int,
    $$n1: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    sva: *$SVA/ArrayType1D,
    $$mv: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    $$eps: $RT,
    $$sfmin: RT,
    $$tol: RT,
    $$nsweep: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           D.data_type == RT &&
           SVA.data_type == RT &&
           V.data_type == RT &&
           WORK.data_type == RT;
} {
    zero_t :: #run zero(RT);
    one_t  :: #run one(RT);
    routine_name :: #run if RT == float32 { "SGSVJ1" } else { "DGSVJ1" };

    applv, rsvec : bool = ---;
    big, bigtheta, large, mxaapq, mxsinj, rootbig, rooteps, rootsfmin, roottol, small, temp1 : RT = ---;
    blskip, emptsw, i, ierr, iswrot, kbl, mvl, nblc, nblr, notrot, pskippped, rowskip, swband : int = ---;

    info.* = 0;
    applv = lsame(jobv, #char "A");
    rsvec = lsame(jobv, #char "V");

    if !(rsvec || applv || lsame(jobv, #char "N")) {
        info.* = -1;
    } else if m < 0 {
        info.* = -2;
    } else if n < 0 || n > m {
        info.* = -3;
    } else if n1 < 0 {
        info.* = -4;
    } else if lda < m {
        info.* = -6;
    } else if (rsvec || applv) && mv < 0 {
        info.* = -9;
    } else if (rsvec && ldv < n) || (applv && ldv < mv) {
        info.* = -11;
    } else if tol <= eps {
        info.* = -14;
    } else if nsweep < 0 {
        info.* = -15;
    } else if lwork < m {
        info.* = -17;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if rsvec {
        mvl = n;
    } else if applv {
        mvl = mv;
    }
    rsvec = rsvec || applv;

    rooteps = sqrt(eps);
    rootsfmin = sqrt(sfmin);
    small = sfmin / eps;
    big = one_t / sfmin;
    rootbig = one_t / rootsfmin;
    large = big / sqrt(ncast(RT, m*n));
    bigtheta = one_t / rooteps;
    roottol = sqrt(tol);

    emptsw = n1*(n - n1);
    notrot = 0;
    swband = 0;
    kbl = min(8, n);
    nblr = n1 / kbl;
    if (nblr*kbl) != n1 then nblr += 1;
    nblc = (n - n1) / kbl;
    if (nblc*kbl) != (n - n1) then nblc += 1;
    blskip = (kbl*kbl) + 1;
    rowskip = min(5, kbl);
    pskippped = 0;
    i = 1;
    iswrot = 0;
    mxaapq = zero_t;
    mxsinj = zero_t;
    temp1 = zero_t;
    ierr = 0;

    // Staged port note:
    // The Jacobi sweep body from *gsvj1.f is still pending.
    _ = a; _ = d; _ = sva; _ = v; _ = work;
    _ = rootbig; _ = large; _ = bigtheta; _ = roottol;
    _ = emptsw; _ = notrot; _ = blskip; _ = mvl; _ = nblc; _ = nblr;
    _ = rowskip; _ = swband; _ = pskippped; _ = i; _ = iswrot; _ = mxaapq; _ = mxsinj; _ = temp1; _ = ierr;
}

_gsvj1_complex :: (
    $$jobv: u8,
    $$m: int,
    $$n: int,
    $$n1: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    sva: *$SVA/ArrayType1D,
    $$mv: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    $$eps: $RT,
    $$sfmin: RT,
    $$tol: RT,
    $$nsweep: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           D.data_type == MT &&
           SVA.data_type == RT &&
           V.data_type == MT &&
           WORK.data_type == MT;
} {
    zero_t :: #run zero(RT);
    one_t  :: #run one(RT);
    routine_name :: #run if MT == C32 { "CGSVJ1" } else { "ZGSVJ1" };

    applv, rsvec : bool = ---;
    big, bigtheta, large, mxaapq, mxsinj, rootbig, rooteps, rootsfmin, roottol, small, temp1 : RT = ---;
    blskip, emptsw, i, ierr, iswrot, kbl, mvl, nblc, nblr, notrot, pskippped, rowskip, swband : int = ---;

    info.* = 0;
    applv = lsame(jobv, #char "A");
    rsvec = lsame(jobv, #char "V");

    if !(rsvec || applv || lsame(jobv, #char "N")) {
        info.* = -1;
    } else if m < 0 {
        info.* = -2;
    } else if n < 0 || n > m {
        info.* = -3;
    } else if n1 < 0 {
        info.* = -4;
    } else if lda < m {
        info.* = -6;
    } else if (rsvec || applv) && mv < 0 {
        info.* = -9;
    } else if (rsvec && ldv < n) || (applv && ldv < mv) {
        info.* = -11;
    } else if tol <= eps {
        info.* = -14;
    } else if nsweep < 0 {
        info.* = -15;
    } else if lwork < m {
        info.* = -17;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if rsvec {
        mvl = n;
    } else if applv {
        mvl = mv;
    }
    rsvec = rsvec || applv;

    rooteps = sqrt(eps);
    rootsfmin = sqrt(sfmin);
    small = sfmin / eps;
    big = one_t / sfmin;
    rootbig = one_t / rootsfmin;
    large = big / sqrt(ncast(RT, m*n));
    bigtheta = one_t / rooteps;
    roottol = sqrt(tol);

    emptsw = n1*(n - n1);
    notrot = 0;
    swband = 0;
    kbl = min(8, n);
    nblr = n1 / kbl;
    if (nblr*kbl) != n1 then nblr += 1;
    nblc = (n - n1) / kbl;
    if (nblc*kbl) != (n - n1) then nblc += 1;
    blskip = (kbl*kbl) + 1;
    rowskip = min(5, kbl);
    pskippped = 0;
    i = 1;
    iswrot = 0;
    mxaapq = zero_t;
    mxsinj = zero_t;
    temp1 = zero_t;
    ierr = 0;

    // Staged port note:
    // The complex Jacobi sweep body from *gsvj1.f is still pending.
    _ = a; _ = d; _ = sva; _ = v; _ = work;
    _ = rootbig; _ = large; _ = bigtheta; _ = roottol;
    _ = emptsw; _ = notrot; _ = blskip; _ = mvl; _ = nblc; _ = nblr;
    _ = rowskip; _ = swband; _ = pskippped; _ = i; _ = iswrot; _ = mxaapq; _ = mxsinj; _ = temp1; _ = ierr;
}
