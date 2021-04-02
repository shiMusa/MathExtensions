




is_quadratic :: (m: Matrix64) -> bool {
    return m.cols == m.rows;
}
is_quadratic :: (m: CMatrix64) -> bool {
    return m.cols == m.rows;
}

is_diagonal_unit :: (m: Matrix64) -> bool {
    res := true;
    n := min(m.cols, m.rows);
    for 0..n-1 {
        if m[it][it] != 1.0 then return false;
    }
    return true;
}
is_diagonal_unit :: (m: CMatrix64) -> bool {
    n := min(m.cols, m.rows);
    for 0..n-1 {
        if m[it][it] != R then return false;
    }
    return true;
}


// checks if values != 0 are only in the bottom left corner
is_left_triangular :: (m: Matrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: i+1..m.cols-1 {
            if m[i][j] != 0.0 then return false;
        }
    }
    return true;
}
is_left_triangular :: (m: CMatrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: i+1..m.cols-1 {
            if m[i][j] != CZERO then return false;
        }
    }
    return true;
}

// checks if values != 0 are only in the top right corner
is_right_triangular :: (m: Matrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: 0..i-1 {
            if m[i][j] != 0.0 then return false;
        }
    }
    return true;
}
is_right_triangular :: (m: CMatrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: 0..i-1 {
            if m[i][j] != CZERO then return false;
        }
    }
    return true;
}

is_right_quasi_triangular :: (R: Matrix64) -> bool {
    for i: 0..R.rows-1 {
        for j: 0..i-2 {
            if i>j+1 { 
                if R[i][j] != 0.0 {
                    print("non-zero element found: R[%][%] = %\n", i, j, R[i][j]);
                    return false;
                }
            }
        }
        if i>=2 {
            // diagonal blocks
            if R[i][i-1] != 0.0 {
                if R[i-1][i-2] != 0.0 || R[i+1][i] != 0.0 {
                    print("diagonal block error at i = %\n", i);
                    return false;
                }
            }
            // 1x1 diagonal blocks invertible
            if R[i-1][i] == 0.0 && R[i+1][i] == 0.0 {
                if R[i][i] == 0.0 {
                    print("1x1 block non-invertible, i = %\n", i);
                    return false;
                }
            }
            // 2x2 diagonal blocks are invertible
            if R[i][i-1] != 0.0 {
                if R[i-1][i-1]*R[i][i] - R[i-1][i]*R[i][i-1] == 0.0 {
                    print("2x2 block non-invertible, i = %\n", i);
                    return false;
                }
            }
        }
    }
    return true;
}
is_right_quasi_triangular :: (R: CMatrix64) -> bool {
    for i: 0..R.rows-1 {
        for j: 0..i-2 {
            if i>j+1 {
                if R[i][j] != CZERO {
                    print("non-zero element found: R[%][%] = %\n", i, j, R[i][j]);
                    return false;
                }
            } 
        }
        if i>=2 {
            // diagonal blocks
            if R[i][i-1] != CZERO {
                if R[i-1][i-2] != CZERO || R[i+1][i] != CZERO {
                    print("diagonal block error at i = %\n", i);
                    return false;
                }
            }
            // 1x1 diagonal blocks invertible
            if R[i-1][i] == CZERO && R[i+1][i] == CZERO {
                if R[i][i] == CZERO {
                    print("1x1 block non-invertible, i = %\n", i);
                    return false;
                }
            }
            // 2x2 diagonal blocks are invertible
            if R[i][i-1] != CZERO {
                if R[i-1][i-1]*R[i][i] - R[i-1][i]*R[i][i-1] == CZERO {
                    print("2x2 block non-invertible, i = %\n", i);
                    return false;
                }
            }
        }
    }
    return true;
}

is_unitary :: (U: CMatrix64) -> bool {
    return is_quadratic(U) && dagger(U)*U == cmat64_id(U.rows);
}

is_left_trapezoidal :: inline (M: CMatrix64) -> bool {
    return M.rows >= M.cols && is_left_triangular(M);
}

is_right_trapezoidal :: inline (M: CMatrix64) -> bool {
    return M.cols >= M.rows && is_right_triangular(M);
}





test_checks :: () {
    print("\n\nChecks -------------------------------\n\n");
    {
        assert(is_quadratic(cmat64_ones(5, 5)));
        assert(!is_quadratic(cmat64_ones(2, 5)));
    }
    {
        m := cmat64_id(3,3);
        assert(is_diagonal_unit(m));
        m[1][1] = 2.0*R;
        assert(!is_diagonal_unit(m));
    }
    {
        m := cmat64_zeros(4, 4);
        m[3][0] = 42.0*R;
        assert(is_left_triangular(m));
        m[0][3] = 42.0*R;
        assert(!is_left_triangular(m));
    }
    {
        m := cmat64_zeros(4, 4);
        m[0][3] = 42.0*R;
        assert(is_right_triangular(m));
        m[3][0] = 42.0*R;
        assert(!is_right_triangular(m));
    }
    {
        m := mat64(4, 4,
            1, 3, 4, 5,
            0, 2, 3, 4,
            0, 1, 2, 3,
            0, 0, 0, 2
        );
        cm := cmat64(m);

        assert(is_right_quasi_triangular(m));
        m = transpose(m);
        assert(!is_right_quasi_triangular(m));

        assert(is_right_quasi_triangular(cm));
        cm = transpose(cm);
        assert(!is_right_quasi_triangular(cm));
    }

    {
        
    }
}