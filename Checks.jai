




// is_quadratic :: (m: Matrix($T)) -> bool {
//     return m.cols == m.rows;
// }

is_diagonal_unit :: (m: $M/MatrixType) -> bool {
    res := true;
    n := min(rows(m), cols(m));
    eins := one(M.data_type);
    for 0..n-1 {
        if get(m,it,it) != eins then return false;
    }
    return true;
}


// checks if values != 0 are only in the bottom left corner
is_left_triangular :: (m: $M/MatrixType) -> bool {
    zro := zero(M.data_type);
    for i: 0..rows(m)-1 {
        for j: i+1..cols(m)-1 {
            #if #run is_float_number(M.data_type) {
                if !approx(get(m,i,j), zro) then return false;
            } else {
                if get(m,i,j) != zro then return false;
            }
        }
    }
    return true;
}

// checks if values != 0 are only in the top right corner
is_right_triangular :: (m: $M/MatrixType) -> bool {
    zro := zero(M.data_type);
    for i: 0..rows(m)-1 {
        for j: 0..i-1 {
            #if #run is_float_number(M.data_type) {
                if !approx(get(m,i,j), zro) then return false;
            } else {
                if get(m,i,j) != zro then return false;
            }
        }
    }
    return true;
}

is_right_quasi_triangular :: (R: $M/MatrixType) -> bool {
    zro := zero(M.data_type);
    for i: 0..rows(R)-1 {
        for j: 0..i-2 {
            if i>j+1 { 
                if get(R,i,j) != zro {
                    println("non-zero element found: R[%][%] = %", i, j, get(R,i,j));
                    return false;
                }
            }
        }
        if i>=2 {
            // diagonal blocks
            if get(R,i,i-1) != zro {
                if get(R,i-1,i-2) != zro || get(R,i+1,i) != zro {
                    println("diagonal block error at i = %", i);
                    return false;
                }
            }
            // 1x1 diagonal blocks invertible
            if get(R,i-1,i) == zro && get(R,i+1,i) == zro {
                if get(R,i,i) == zro {
                    println("1x1 block non-invertible, i = %", i);
                    return false;
                }
            }
            // 2x2 diagonal blocks are invertible
            if get(R,i,i-1) != zro {
                if get(R,i-1,i-1)*get(R,i,i) - get(R,i-1,i)*get(R,i,i-1) == zro {
                    println("2x2 block non-invertible, i = %", i);
                    return false;
                }
            }
        }
    }
    return true;
}

is_unitary :: (m: $M/MatrixType) -> bool {
    dm := dhmat(M.data_type, rows(m), cols(m));
    defer free(dm);
    copy(*dm, m);
    dagger(*d);
    dmm := dhmat(M.data_type, cols(m), cols(m));
    defer free(dmm);
    mat_mul_mat(dm, m, *dmm);
    return is_quadratic(m) && is_diagonal_unit(dmm);
}

// more rows than columns
is_left_trapezoidal :: inline (m: $M/MatrixType) -> bool {
    return rows(m) >= cols(m) && is_left_triangular(m);
}

// more columns than rows
is_right_trapezoidal :: inline (m: $M/MatrixType) -> bool {
    return cols(m) >= rows(m) && is_right_triangular(m);
}





test_checks :: () {
    println_push("Checks", color = .FG_WHITE);
    {
        m : DenseMatrix(C64, 3, 3);
        mat_init_unit(*m);
        assert(is_diagonal_unit(m));
        set(*m,1,1, 2.0*R);
        assert(!is_diagonal_unit(m));
    }
    {
        m : DenseMatrix(C64, 4, 4);
        set(*m,3,0, 42.0*R);
        assert(is_left_triangular(m));
        set(*m,0,3, 42.0*R);
        assert(!is_left_triangular(m));
    }
    {
        m : DenseMatrix(C64, 4, 4);
        set(*m,0,3, 42.0*R);
        assert(is_right_triangular(m));
        set(*m,3,0, 42.0*R);
        assert(!is_right_triangular(m));
    }
    {
        m : DenseMatrix(int, 4, 4);
        mat_init(*m,
            1, 3, 4, 5,
            0, 2, 3, 4,
            0, 1, 2, 3,
            0, 0, 0, 2
        );
        cm : DenseMatrix(C64, 4, 4);
        cast_to(*cm, m);

        assert(is_right_quasi_triangular(m));
        transpose(*m);
        assert(!is_right_quasi_triangular(m));

        assert(is_right_quasi_triangular(cm));
        transpose(*cm);
        assert(!is_right_quasi_triangular(cm));
    }
}