
#load "module.jai";



main :: () {
    print("\n\ntesting some math extensions\n\n");

    test_utils();
    test_complex();

    test_elementary();

    test_polynomials();
    test_approximations();
    
    // test_real_matrix();
    // test_real_vector();

    // test_complex_vector();
    // test_complex_matrix();
    
    // test_checks();
    // test_linear_algebra();

    test_vector();
    test_matrix();
}