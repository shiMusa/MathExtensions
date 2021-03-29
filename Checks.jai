




mat64_is_quadratic :: (m: Matrix64) -> bool {
    return m.cols == m.rows;
}
cmat64_is_quadratic :: (m: CMatrix64) -> bool {
    return m.cols == m.rows;
}

mat64_is_diagonal_unit :: (m: Matrix64) -> bool {
    res := true;
    n := min(m.cols, m.rows);
    for 0..n-1 {
        if m[it][it] != 1.0 then return false;
    }
    return true;
}
cmat64_is_diagonal_unit :: (m: CMatrix64) -> bool {
    n := min(m.cols, m.rows);
    for 0..n-1 {
        if m[it][it] != R then return false;
    }
    return true;
}

mat64_is_left_triangular :: (m: Matrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: i+1..m.cols-1 {
            if m[i][j] != 0.0 then return false;
        }
    }
    return true;
}
cmat64_is_left_triangular :: (m: CMatrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: i+1..m.cols-1 {
            if m[i][j] != CZERO then return false;
        }
    }
    return true;
}

mat64_is_right_triangular :: (m: Matrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: 0..i-1 {
            if m[i][j] != 0.0 then return false;
        }
    }
    return true;
}
cmat64_is_right_triangular :: (m: CMatrix64) -> bool {
    for i: 0..m.rows-1 {
        for j: 0..i-1 {
            if m[i][j] != CZERO then return false;
        }
    }
    return true;
}

mat64_is_right_quasi_triangular :: (R: Matrix64) -> bool {
    for i: 0..R.rows-1 {
        for j: 0..R.cols-1 {
            if i>j+1 && R[i][j] != 0.0 then return false;
        }
        if i>=2 {
            // diagonal blocks
            if R[i-1][i] != 0.0 {
                if R[i-2][i-1] != 0.0 || R[i][i+1] != 0.0 then return false;
            }
            // 1x1 diagonal blocks invertible
            if R[i-1][i] == 0.0 && R[i+1][i] == 0.0 {
                if R[i][i] == 0.0 then return false;
            }
            // 2x2 diagonal blocks are invertible
            if R[i][i-1] != 0.0 {
                if R[i-1][i-1]*R[i][i] - R[i-1][i]*R[i][i-1] == 0.0 then return false;
            }
        }
    }
    return true;
}
cmat64_is_right_quasi_triangular :: (R: CMatrix64) -> bool {
    for i: 0..R.rows-1 {
        for j: 0..R.cols-1 {
            if i>j+1 && R[i][j] != CZERO then return false;
        }
        if i>=2 {
            // diagonal blocks
            if R[i-1][i] != CZERO {
                if R[i-2][i-1] != CZERO || R[i][i+1] != CZERO then return false;
            }
            // 1x1 diagonal blocks invertible
            if R[i-1][i] == CZERO && R[i+1][i] == CZERO {
                if R[i][i] == CZERO then return false;
            }
            // 2x2 diagonal blocks are invertible
            if R[i][i-1] != CZERO {
                if R[i-1][i-1]*R[i][i] - R[i-1][i]*R[i][i-1] == CZERO then return false;
            }
        }
    }
    return true;
}



test_checks :: () {
    {
        assert()
    }
}