
#load "module.jai";



main :: () {
    print("\n\ntesting some math extensions\n\n");

    test_utils();
    // test_complex();
    test_complex_generic();

    test_elementary();

    test_polynomials();
    test_approximations();
    
    // test_scalar();
    // test_vector();
    // test_matrix();

    // test_checks();
    // test_linear_algebra();
}