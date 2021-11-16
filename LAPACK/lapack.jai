

#load "Level1/_axpy.jai";
#load "Level1/_copy.jai";
#load "Level1/_dotc.jai";
#load "Level1/_dotu.jai";
#load "Level1/_scal.jai";
#load "Level1/_srot.jai";
#load "Level1/_swap.jai";
#load "Level1/_asum.jai";


test_lapack :: () {
    println_push("LAPACK", color = .FG_WHITE);
    test_lapack_axpy();
    test_lapack_copy();
    test_lapack_dotc();
    test_lapack_dotu();
    test_lapack_scal();
    test_lapack_srot();
    test_lapack_swap();
    test_lapack_asum();
}