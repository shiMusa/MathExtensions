#import "Basic";
Math :: #import "Math";


#load "Utils.jai";
#load "Complex.jai";
#load "Polynomials.jai";



main :: () {
    print("\n\ntesting some math extensions\n\n");

    test_utils();
    test_complex();
    test_polynomials();
}