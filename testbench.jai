Complex64 :: struct {
    real: float64;
    imag: float64;
}

foo :: (numbers: ..float64) {}
foo :: (numbers: ..Complex64) {}

main :: () {
    // both lines lead to the same error!
    foo(cast(float64)1.0, 2.0, 3.0, 4.0);
    // foo(Complex64.{1.0, 0.0}, Complex64.{2.0, 0.0});
}