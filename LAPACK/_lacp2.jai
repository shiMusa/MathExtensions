// ? LAPACK *lacp2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clacp2/zlacp2
// *LACP2 copies all or part of a real matrix to a complex matrix.

clacp2 :: #bake_arguments _lacp2(RT=float32, MT=C32);
zlacp2 :: #bake_arguments _lacp2(RT=float64, MT=C64);

#scope_file
_lacp2 :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == RT &&
           B.data_type == MT;
} {
    i, j : int = ---;

    if lsame(uplo, #char "U") {
        for j:1..n {
            for i:1..min(j, m) {
                fset(b, i, j, ncast(MT, fget(a, i, j)));
            }
        }
    } else if lsame(uplo, #char "L") {
        for j:1..n {
            for i:j..m {
                fset(b, i, j, ncast(MT, fget(a, i, j)));
            }
        }
    } else {
        for j:1..n {
            for i:1..m {
                fset(b, i, j, ncast(MT, fget(a, i, j)));
            }
        }
    }
}
