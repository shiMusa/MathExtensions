
// ? LAPACK drotmg
// ? http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga13e351a3dfafa2cd8dc5302dcf53f69a.html#ga13e351a3dfafa2cd8dc5302dcf53f69a

// ? Original comment
// CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
// THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*>    DY2)**T.
// WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..

// DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0

//     (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
// H=(          )    (          )    (          )    (          )
//     (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
// LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
// RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
// VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)

// THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
// INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
// OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.

#scope_file
gam    : float64 : 4096.0;
gamsq  : float64 : 16777216.0;
rgamsq : float64 : 5.9604645e-8;
#scope_export

drotmg :: (dd1: *float64, dd2: *float64, dx1: *float64, dy1: float64, dparam: [5]float64) {
    dflag : float64 = ---;
    dh11  : float64 = ---;
    dh12  : float64 = ---;
    dh21  : float64 = ---;
    dh22  : float64 = ---;

    if <<dd1 < 0.0 {
        // go zero-h-d-and-dx1..
        dflag = -1.0;
        dh11 = 0.0;
        dh12 = 0.0;
        dh21 = 0.0;
        dh22 = 0.0;

        <<dd1 = 0.0;
        <<dd2 = 0.0;
        <<dx1 = 0.0;
    } else {
        // case-dd1-nonnegative
        dp2 := <<dd2 * dy1;
        if dp2 == 0.0 {
            dparam[0] = -2.0;
            return;
        }
        // regular-case..
        dp1 := <<dd1 * <<dx1;
        dq2 := dp2 * dy1;
        dq1 := dp1 * <<dx1;

        if dabs(dq1) > dabs(dq2) {
            dh21 = -dy1/<<dx1;
            dh12 = dp2/dp1;

            du := 1.0 - dh12*dh21;

            if du > 0.0 {
                dflag = 0.0;
                <<dd1 /= du;
                <<dd2 /= du;
                <<dx1 /= du;
            } else {
                // This code path is here for safety. We do not expect this
                // condition to ever hold except in edge cases with rounding
                // errors. See DOI: 10.1145/35841.355847
                dflag = -1.0;
                dh11 = 0.0;
                dh12 = 0.0;
                dh21 = 0.0;
                dh22 = 0.0;

                <<dd1 = 0.0;
                <<dd2 = 0.0;
                <<dx1 = 0.0;
            }
        } else {
            if dq2 < 0.0 {
                // go zero-h-d-and-dx1..
                dflag = -1.0;
                dh11 = 0.0;
                dh12 = 0.0;
                dh21 = 0.0;
                dh22 = 0.0;

                <<dd1 = 0.0;
                <<dd2 = 0.0;
                <<dx1 = 0.0;
            } else {
                dflag = 1.0;
                dh11 = dp1/dp2;
                dh22 = <<dx1/dy1;
                du := 1.0 + dh11*dh22;
                dtmp := <<dd1/du;
                <<dd2 = <<dd1/du;
                <<dd1 = dtmp;
                <<dx1 = dy1*du;
            }
        }

        // procedure..scale-check
        if <<dd1 != 0.0 {
            while <<dd1 <= rgamsq || <<dd1 >= gamsq {
                if dflag == 0.0 {
                    dh11 = 1.0;
                    dh22 = 1.0;
                    dflag = -1.0;
                } else {
                    dh21 = -1.0;
                    dh12 = 1.0;
                    dflag = -1.0;
                }
                if <<dd1 <= rgamsq {
                    <<dd1 *= gam*gam;
                    <<dx1 /= gam;
                    dh11 /= gam;
                    dh12 /= gam;
                } else {
                    <<dd1 /= gam*gam;
                    <<dx1 *= gam;
                    dh11 *= gam;
                    dh12 *= gam;
                }
            }
        }
        if <<dd2 != 0.0 {
            while dabs(<<dd2) <= rgamsq || dabs(<<dd2) >= gamsq {
                if dflag == 0.0 {
                    dh11 = 1.0;
                    dh22 = 1.0;
                    dflag = -1.0;
                } else {
                    dh21 = -1.0;
                    dh12 = 1.0;
                    dflag = -1.0;
                }
                if dabs(<<dd2) <= rgamsq {
                    <<dd2 *= gam*gam;
                    dh21 /= gam;
                    dh22 /= gam;
                } else {
                    <<dd2 /= gam*gam;
                    dh21 *= gam;
                    dh22 *= gam;
                }
            }
        }
    }

    if dflag < 0.0 {
        dparam[1] = dh11;
        dparam[2] = dh21;
        dparam[3] = dh12;
        dparam[4] = dh22;
    } else if dflag == 0.0 {
        dparam[2] = dh21;
        dparam[3] = dh12;
    } else {
        dparam[1] = dh11;
        dparam[4] = dh22;
    }
    dparam[0] = dflag;
}