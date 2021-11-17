

#load "Level1/_axpy.jai";
#load "Level1/_copy.jai";
#load "Level1/_dotc.jai";
#load "Level1/_dotu.jai";
#load "Level1/_scal.jai";
#load "Level1/_rot.jai";
#load "Level1/_swap.jai";
#load "Level1/_asum.jai";
#load "Level1/_abs.jai";
#load "Level1/drotm.jai";
#load "Level1/drotmg.jai";


test_lapack :: () {
    println_push("LAPACK", color = .FG_WHITE);
    test_lapack_axpy();
    test_lapack_copy();
    test_lapack_dotc();
    test_lapack_dotu();
    test_lapack_scal();
    test_lapack_rot();
    test_lapack_swap();
    test_lapack_asum();
}