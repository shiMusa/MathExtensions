




is_quadratic :: (m: Matrix($T)) -> bool {
    return m.cols == m.rows;
}

is_diagonal_unit :: (m: Matrix($T)) -> bool {
    res := true;
    n := min(m.cols, m.rows);
    for 0..n-1 {
        if m[it][it] != one(T) then return false;
    }
    return true;
}


// checks if values != 0 are only in the bottom left corner
is_left_triangular :: (m: Matrix($T)) -> bool {
    for i: 0..m.rows-1 {
        for j: i+1..m.cols-1 {
            if m[i][j] != zero(T) then return false;
        }
    }
    return true;
}

// checks if values != 0 are only in the top right corner
is_right_triangular :: (m: Matrix($T)) -> bool {
    for i: 0..m.rows-1 {
        for j: 0..i-1 {
            if m[i][j] != zero(T) then return false;
        }
    }
    return true;
}

is_right_quasi_triangular :: (R: Matrix($T)) -> bool {
    for i: 0..R.rows-1 {
        for j: 0..i-2 {
            if i>j+1 { 
                if R[i][j] != zero(T) {
                    print("non-zero element found: R[%][%] = %\n", i, j, R[i][j]);
                    return false;
                }
            }
        }
        if i>=2 {
            // diagonal blocks
            if R[i][i-1] != zero(T) {
                if R[i-1][i-2] != zero(T) || R[i+1][i] != zero(T) {
                    print("diagonal block error at i = %\n", i);
                    return false;
                }
            }
            // 1x1 diagonal blocks invertible
            if R[i-1][i] == zero(T) && R[i+1][i] == zero(T) {
                if R[i][i] == zero(T) {
                    print("1x1 block non-invertible, i = %\n", i);
                    return false;
                }
            }
            // 2x2 diagonal blocks are invertible
            if R[i][i-1] != zero(T) {
                if R[i-1][i-1]*R[i][i] - R[i-1][i]*R[i][i-1] == zero(T) {
                    print("2x2 block non-invertible, i = %\n", i);
                    return false;
                }
            }
        }
    }
    return true;
}

is_unitary :: (U: Matrix($T)) -> bool {
    return is_quadratic(U) && is_diagonal_unit(dagger(U)*U);
}

is_left_trapezoidal :: inline (M: Matrix($T)) -> bool {
    return M.rows >= M.cols && is_left_triangular(M);
}

is_right_trapezoidal :: inline (M: Matrix($T)) -> bool {
    return M.cols >= M.rows && is_right_triangular(M);
}





test_checks :: () {
    print("\n\nChecks -------------------------------\n\n");
    {
        assert(is_quadratic(mat_ones(Complex64, 5, 5)));
        assert(!is_quadratic(mat_ones(Complex64, 2, 5)));
    }
    {
        m := mat_id(Complex64, 3,3);
        assert(is_diagonal_unit(m));
        m[1][1] = 2.0*R;
        assert(!is_diagonal_unit(m));
    }
    {
        m := mat_zeros(Complex64, 4, 4);
        m[3][0] = 42.0*R;
        assert(is_left_triangular(m));
        m[0][3] = 42.0*R;
        assert(!is_left_triangular(m));
    }
    {
        m := mat_zeros(Complex64, 4, 4);
        m[0][3] = 42.0*R;
        assert(is_right_triangular(m));
        m[3][0] = 42.0*R;
        assert(!is_right_triangular(m));
    }
    {
        m := mat(4, 4,
            1, 3, 4, 5,
            0, 2, 3, 4,
            0, 1, 2, 3,
            0, 0, 0, 2
        );
        cm := mat(Complex64, m);

        assert(is_right_quasi_triangular(m));
        m = transpose(m);
        assert(!is_right_quasi_triangular(m));

        assert(is_right_quasi_triangular(cm));
        cm = transpose(cm);
        assert(!is_right_quasi_triangular(cm));
    }
}