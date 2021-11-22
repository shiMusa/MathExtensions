
// ? LAPACK drotm
// ? http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga0d579faa55a493933032c5dcddbe5f4f.html#ga0d579faa55a493933032c5dcddbe5f4f

// ? Original comment:
// APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
//
// (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
// (DY**T)
//
// DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
// LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
// WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
//
// DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
//
//   (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
// H=(          )    (          )    (          )    (          )
//   (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
// SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.

drotm :: ($$n: int, dx: *$X/ArrayType1D, $$incx: int, dy: *$Y/ArrayType1D, $$incy: int, dparam: [5]float64) {
    #assert(X.data_type == float64 && Y.data_type == float64);
    dflag := dparam[0];
    if n <= 0 || dflag+2.0 == 0.0 then return;
    if incx == incy && incx > 0 {
        if dflag < 0.0 {
            dh11 := dparam[1];
            dh12 := dparam[3];
            dh21 := dparam[2];
            dh22 := dparam[4];
            ix := 0;
            for i:0..n-1 {
                w := get(dx,ix);
                z := get(dy,ix);
                set(dx,ix, w*dh11 + z*dh12);
                set(dy,ix, w*dh21 + z*dh22);
                ix += incx;
            }
        } else if dflag == 0.0 {
            dh12 := dparam[3];
            dh21 := dparam[2];
            ix := 0;
            for i: 0..n-1 {
                w := get(dx,ix);
                z := get(dy,ix);
                set(dx,ix, w + z*dh12);
                set(dy,ix, w*dh21 + z);
                ix += incx;
            }
        } else {
            dh11 := dparam[1];
            dh22 := dparam[4];
            ix := 0;
            for i: 0..n-1 {
                w := get(dx,ix);
                z := get(dy,ix);
                set(dx,ix,  w*dh11 + z);
                set(dy,ix, -w + dh22*z);
            }
        }
    } else {
        kx := 0;
        ky := 0;
        if incx < 0 then kx = (1-n)*incx;
        if incy < 0 then ky = (1-n)*incy;

        if dflag < 0.0 {
            dh11 := dparam[1];
            dh12 := dparam[3];
            dh21 := dparam[2];
            dh22 := dparam[4];
            for i: 0..n-1 {
                w := get(dx,kx);
                z := get(dy,ky);
                set(dx,kx, w*dh11 + z*dh12);
                set(dy,ky, w*dh21 + z*dh22);
                kx += incx;
                ky += incy;
            }
        } else if dflag == 0.0 {
            dh12 := dparam[3];
            dh21 := dparam[2];
            for i: 0..n-1 {
                w := get(dx,kx);
                z := get(dy,ky);
                set(dx,kx, w + z*dh12);
                set(dy,ky, w*dh21 + z);
                kx += incx;
                ky += incy;
            }
        } else {
            dh11 := dparam[1];
            dh22 := dparam[4];
            for i: 0..n-1 {
                w := get(dx,kx);
                z := get(dy,ky);
                set(dx,kx,  w*dh11 + z);
                set(dy,ky, -w + dh22*z);
                kx += incx;
                ky += incy;
            }
        }
    }
}