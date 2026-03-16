// ? LAPACK BLAT1 driver

#load "cblat1.jai";

cblat1 :: () -> bool {
    return cblat1_impl();
}
