#load "Console.jai";

#load "Utils.jai";
#load "Approx.jai";

#load "Types.jai";
#load "Elementary.jai";
#load "Complex.jai";
#load "Quaternion.jai";
#load "Polynomials.jai";
#load "Approximations.jai";

#load "VectorTypes.jai";
#load "MatrixTypes.jai";
#load "Vector.jai";
#load "Matrix.jai";

#load "Checks.jai";
#load "LinearAlgebra.jai";

#load "LAPACK/lapack.jai";

#module_parameters(CHECKS := true);


run_all_tests :: (verbose := false) {
    println("\nRunning tests ///////////////////////////////////////////////////\n", color = .FG_RED);

    test_console();

    test_utils();
    test_complex();
    test_quaternion();

    test_elementary();

    test_polynomials();
    test_approximations();
    
    test_MatrixTypes();

    test_vector();
    test_matrix();

    test_checks();
    test_linear_algebra();

    test_lapack();
}


#scope_module

#import "Basic";
Math   :: #import "Math";
Pool   :: #import "Pool";
String :: #import "String";
File   :: #import "File";
