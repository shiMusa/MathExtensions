// ? LAPACK *gsvj0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgsvj0/cgsvj0/dgsvj0/zgsvj0
// *GSVJ0 is a Jacobi-SVD preprocessing kernel used by *GESVJ.

sgsvj0 :: #bake_arguments _gsvj0_real(RT=float32);
dgsvj0 :: #bake_arguments _gsvj0_real(RT=float64);
cgsvj0 :: #bake_arguments _gsvj0_complex(RT=float32, MT=C32);
zgsvj0 :: #bake_arguments _gsvj0_complex(RT=float64, MT=C64);

#scope_file
_gsvj0_real :: (
    $$jobv: u8,
    $$m: int,
    $$n: int,
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
    half_t :: #run one(RT) / ncast(RT, 2);
    routine_name :: #run if RT == float32 { "SGSVJ0" } else { "DGSVJ0" };

    applv, rsvec : bool = ---;
    aapp, big, bigtheta, mxaapq, mxsinj, rootbig, rooteps, rootsfmin, roottol, small, temp1 : RT = ---;
    blskip, emptsw, i, ierr, iswrot, kbl, lkahead, mvl, nbl, notrot, pskippped, rowskip, swband : int = ---;
    fastr : [5] RT;

    info.* = 0;
    applv = lsame(jobv, #char "A");
    rsvec = lsame(jobv, #char "V");

    if !(rsvec || applv || lsame(jobv, #char "N")) {
        info.* = -1;
    } else if m < 0 {
        info.* = -2;
    } else if n < 0 || n > m {
        info.* = -3;
    } else if lda < m {
        info.* = -5;
    } else if (rsvec || applv) && mv < 0 {
        info.* = -8;
    } else if (rsvec && ldv < n) || (applv && ldv < mv) {
        info.* = -10;
    } else if tol <= eps {
        info.* = -13;
    } else if nsweep < 0 {
        info.* = -14;
    } else if lwork < m {
        info.* = -16;
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
    bigtheta = one_t / rooteps;
    roottol = sqrt(tol);

    emptsw = (n*(n - 1)) / 2;
    notrot = 0;
    fastr[0] = zero_t;
    swband = 0;
    kbl = min(8, n);
    nbl = n / kbl;
    if (nbl*kbl) != n then nbl += 1;
    blskip = (kbl*kbl) + 1;
    rowskip = min(5, kbl);
    lkahead = 1;
    pskippped = 0;
    i = 1;
    iswrot = 0;
    mxaapq = zero_t;
    mxsinj = zero_t;
    aapp = zero_t;
    temp1 = zero_t;
    ierr = 0;

    // Staged port note:
    // The Jacobi sweep body from *gsvj0.f is still pending.
    // This file lands the family interface, argument validation,
    // and the source-shaped setup needed for the next pass.
    _ = d; _ = sva; _ = v; _ = work;
    _ = half_t; _ = rootbig; _ = bigtheta; _ = roottol;
    _ = emptsw; _ = notrot; _ = blskip; _ = lkahead; _ = mvl; _ = nbl;
    _ = rowskip; _ = swband; _ = pskippped; _ = i; _ = iswrot; _ = mxaapq; _ = mxsinj; _ = aapp; _ = temp1; _ = ierr;
}

_gsvj0_complex :: (
    $$jobv: u8,
    $$m: int,
    $$n: int,
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
    routine_name :: #run if MT == C32 { "CGSVJ0" } else { "ZGSVJ0" };

    applv, rsvec : bool = ---;
    big, bigtheta, mxaapq, mxsinj, rootbig, rooteps, rootsfmin, roottol, small, temp1 : RT = ---;
    blskip, emptsw, i, ierr, iswrot, kbl, lkahead, mvl, nbl, notrot, pskippped, rowskip, swband : int = ---;

    info.* = 0;
    applv = lsame(jobv, #char "A");
    rsvec = lsame(jobv, #char "V");

    if !(rsvec || applv || lsame(jobv, #char "N")) {
        info.* = -1;
    } else if m < 0 {
        info.* = -2;
    } else if n < 0 || n > m {
        info.* = -3;
    } else if lda < m {
        info.* = -5;
    } else if (rsvec || applv) && mv < 0 {
        info.* = -8;
    } else if (rsvec && ldv < n) || (applv && ldv < mv) {
        info.* = -10;
    } else if tol <= eps {
        info.* = -13;
    } else if nsweep < 0 {
        info.* = -14;
    } else if lwork < m {
        info.* = -16;
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
    bigtheta = one_t / rooteps;
    roottol = sqrt(tol);

    emptsw = (n*(n - 1)) / 2;
    notrot = 0;
    swband = 0;
    kbl = min(8, n);
    nbl = n / kbl;
    if (nbl*kbl) != n then nbl += 1;
    blskip = (kbl*kbl) + 1;
    rowskip = min(5, kbl);
    lkahead = 1;
    pskippped = 0;
    i = 1;
    iswrot = 0;
    mxaapq = zero_t;
    mxsinj = zero_t;
    temp1 = zero_t;
    ierr = 0;

    // Staged port note:
    // The complex Jacobi sweep body from *gsvj0.f is still pending.
    _ = a; _ = d; _ = sva; _ = v; _ = work;
    _ = rootbig; _ = bigtheta; _ = roottol;
    _ = emptsw; _ = notrot; _ = blskip; _ = lkahead; _ = mvl; _ = nbl;
    _ = rowskip; _ = swband; _ = pskippped; _ = i; _ = iswrot; _ = mxaapq; _ = mxsinj; _ = temp1; _ = ierr;
}
