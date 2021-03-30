#load "module.jai";


Foo :: struct (T: Type) {
    type : Type : T;
    data : [..]T;
}

main :: () {
    x : u8 = zero(u8);
    y := zero(Complex64);
    print("% %\n", x, y);

    one_u8 := ncast(u8, 0);
    one_float64 := ncast(float64, 0);
    one_complex := ncast(Complex64, 0);
    print("% % %\n", one_u8, one_float64, one_complex);
}