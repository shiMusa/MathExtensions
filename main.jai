
#load "module.jai";



main :: () {
    println("\nRunning tests ///////////////////////////////////////////////////\n", color = .FG_RED);

    test_utils();
    test_complex();

    test_elementary();

    test_polynomials();
    test_approximations();
    
    test_vector();
    test_matrix();

    test_checks();
    test_linear_algebra();
}