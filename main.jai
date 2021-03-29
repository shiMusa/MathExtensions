#import "Basic";
Math :: #import "Math";


#load "Utils.jai";
#load "Complex.jai";
#load "Polynomials.jai";
#load "Approximations.jai";
#load "Matrix.jai";
#load "ComplexMatrix.jai";
#load "Vector.jai";
#load "ComplexVector.jai";
#load "LinearAlgebra.jai";



main :: () {
    print("\n\ntesting some math extensions\n\n");

    test_utils();
    // test_complex();
    // test_polynomials();
    // test_approximations();
    test_matrix();
    test_complex_matrix();
    // test_vector();
    test_linear_algebra();
}