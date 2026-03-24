// ? LAPACK *lapmt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slapmt/clapmt/dlapmt/zlapmt
// *LAPMT performs a forward or backward permutation of the columns of a
// matrix.

slapmt :: #bake_arguments _lapmt(MT=float32);
dlapmt :: #bake_arguments _lapmt(MT=float64);
clapmt :: #bake_arguments _lapmt(MT=C32);
zlapmt :: #bake_arguments _lapmt(MT=C64);

#scope_file
_lapmt :: (
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
    i, ii, j, in : int = ---;
    temp : MT = ---;

    if n <= 1 {
        return;
    }

    for i:1..n {
        fset(k, i, -fget(k.*, i));
    }

    if forwrd {
        for i:1..n {
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

                for ii:1..m {
                    temp = fget(x.*, ii, j);
                    fset(x, ii, j, fget(x.*, ii, in));
                    fset(x, ii, in, temp);
                }

                fset(k, in, -fget(k.*, in));
                j = in;
                in = fget(k.*, in);
            }
        }
    } else {
        for i:1..n {
            if fget(k.*, i) > 0 {
                continue;
            }

            fset(k, i, -fget(k.*, i));
            j = fget(k.*, i);
            while j != i {
                for ii:1..m {
                    temp = fget(x.*, ii, i);
                    fset(x, ii, i, fget(x.*, ii, j));
                    fset(x, ii, j, temp);
                }

                fset(k, j, -fget(k.*, j));
                j = fget(k.*, j);
            }
        }
    }
}
