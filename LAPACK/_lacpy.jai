// ? LAPACK *lacpy()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slacpy/dlacpy/clacpy/zlacpy
// *LACPY copies all or part of one matrix to another.

slacpy :: #bake_arguments _lacpy(MT=float32);
dlacpy :: #bake_arguments _lacpy(MT=float64);
clacpy :: #bake_arguments _lacpy(MT=C32);
zlacpy :: #bake_arguments _lacpy(MT=C64);

#scope_file
_lacpy :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    i, j : int = ---;

    if lsame(uplo, #char "U") {
        for j:1..n {
            for i:1..min(j, m) {
                fset(b, i, j, fget(a, i, j));
            }
        }
    } else if lsame(uplo, #char "L") {
        for j:1..n {
            for i:j..m {
                fset(b, i, j, fget(a, i, j));
            }
        }
    } else {
        for j:1..n {
            for i:1..m {
                fset(b, i, j, fget(a, i, j));
            }
        }
    }
}
