
#load "module.jai";



main :: () {
    print("\n\ntesting some math extensions\n\n");

    test_utils();
    test_complex();

    test_elementary();

    test_polynomials();
    test_approximations();
    
    test_checks();
    test_linear_algebra();

    test_scalar();
    test_vector();
    test_matrix();
}

#run main();