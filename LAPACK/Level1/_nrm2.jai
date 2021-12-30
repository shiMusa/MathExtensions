
// ? LAPACK *nrm2
// ? http://www.netlib.org/lapack/explore-html/df/d28/group__single__blas__level1_gaee5779d5d216a7cd8cf83488fb6bb175.html#gaee5779d5d216a7cd8cf83488fb6bb175


// ? original comment scnrm2
// SCNRM2 returns the euclidean norm of a vector via the function
// name, so that
// 
// SCNRM2 := sqrt( x**H*x )

//  Anderson E. (2017)
// Algorithm 978: Safe Scaling in the Level 1 BLAS
// ACM Trans Math Softw 44:1--28
// https://doi.org/10.1145/3061665
// 
// Blue, James L. (1978)
// A Portable Fortran Program to Find the Euclidean Norm of a Vector
// ACM Trans Math Softw 4:15--23
// https://doi.org/10.1145/355769.355771

dznrm2 :: #bake_arguments _nrm2(R=float64);
dnrm2  :: dznrm2;
scnrm2 :: #bake_arguments _nrm2(R=float32);
snrm   :: scnrm2;
_nrm2 :: ($R: Type, n: int, x: $X/ArrayType1D, incx: int) -> R {

    XT   :: #run get_number_inner_type(X.data_type);
    T    :: UpCast(XT, float).T;
    one  :: one(R);
    zero :: zero(R);
    #if R == float {
        maxn :: Math.FLOAT32_MAX;
    } else {
        maxn :: Math.FLOAT64_MAX;
    }

    // local scalars
    ix : int = ---;
    notbig : bool = ---;
    abig, amed, asml, ax, scl, sumsq, ymax, ymin : float64 = ---;

    // quick return if possible
    if n <= 0 return 0.0;

    scl = one;
    sumsq = zero;

    //  Compute the sum of squares in 3 accumulators:
    //     abig -- sums of squares scaled down to avoid overflow
    //     asml -- sums of squares scaled up to avoid underflow
    //     amed -- sums of squares that do not require scaling
    //  The thresholds and multipliers are
    //     tbig -- values bigger than this are scaled down by sbig
    //     tsml -- values smaller than this are scaled up by ssml

    notbig = true;
    asml = zero;
    amed = zero;
    abig = zero;
    ix = 0;

    if incx < 0 then ix = (n-1)*(-incx);
    for i: 0..n-1 {
        ax = abs(real(get(x,ix)));
        
        if ax > tbig {
            abig += (ax*sbig)*(ax*sbig);
            nobig = false;
        } else if ax < tsml {
            if notbig then asml += (ax*ssml)*(ax*ssml);
        } else {
            amed += ax*ax;
        }
        
        ax = abs(imag(get(x,ix)));

        if ax > tbig {
            abig += (ax*sbig)*(ax*sbig);
            notbig = false;
        } else if ax < tsml {
            if notbig then asml += (ax*ssml)*(ax*ssml);
        } else {
            amed += ax*ax;
        }

        ix += incx;
    }

    // Combine abig and amed or amed and asml if more than one
    // accumulator was used.

    if abig > zero {
        // Combine abig and amed if abig > 0.
        if amed > zero || amed > maxn /* || amed != amed WTF? */ {
            abig += (amed*sbig)*sbig;
        }
        scl = one/sbig;
        sumsq = abig;
    } else if asml > zero {
        // Combine amed and asml if asml > 0.
        if amed > zero || amed > maxn /* || amed != amed WTF? */ {
            amed = sqrt(amed);
            asml = sqrt(asml)/ssml;
            if asml > amed {
                ymin = amed;
                ymax = asml;
            } else {
                ymin = asml;
                ymax = amed;
            }
            scl = one;
            sumsq = ymax*ymax*(one + (ymin/ymax)*(ymin/ymax));
        } else {
            scl = one / ssml;
            sumsq = asml;
        }
    } else {
        // Otherwise all values are mid-range
        scl = one;
        sumsq = amed;
    }

    return scl * sqrt(sumsq);
}