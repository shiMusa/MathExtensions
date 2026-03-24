// ? LAPACK *lapmr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slapmr/clapmr/dlapmr/zlapmr
// *LAPMR rearranges rows of a matrix as specified by a permutation
// vector.

slapmr :: #bake_arguments _lapmr(MT=float32);
dlapmr :: #bake_arguments _lapmr(MT=float64);
clapmr :: #bake_arguments _lapmr(MT=C32);
zlapmr :: #bake_arguments _lapmr(MT=C64);

#scope_file
_lapmr :: (
    $$forwrd: bool,
    $$m: int,
    $$n: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    k: *$K/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           K.data_type == int;
} {
    i, in, j, jj : int = ---;
    temp : MT = ---;

    if m <= 1 {
        return;
    }

    for i:1..m {
        fset(k, i, -fget(k.*, i));
    }

    if forwrd {
        for i:1..m {
            if fget(k.*, i) > 0 {
                continue;
            }

            j = i;
            fset(k, j, -fget(k.*, j));
            in = fget(k.*, j);

            while true {
                if fget(k.*, in) > 0 {
                    break;
                }

                for jj:1..n {
                    temp = fget(x.*, j, jj);
                    fset(x, j, jj, fget(x.*, in, jj));
                    fset(x, in, jj, temp);
                }

                fset(k, in, -fget(k.*, in));
                j = in;
                in = fget(k.*, in);
            }
        }
    } else {
        for i:1..m {
            if fget(k.*, i) > 0 {
                continue;
            }

            fset(k, i, -fget(k.*, i));
            j = fget(k.*, i);
            while j != i {
                for jj:1..n {
                    temp = fget(x.*, i, jj);
                    fset(x, i, jj, fget(x.*, j, jj));
                    fset(x, j, jj, temp);
                }

                fset(k, j, -fget(k.*, j));
                j = fget(k.*, j);
            }
        }
    }
}
