




Matrix64 :: struct {
    rows: int;
    cols: int;
    data: [..]float64;
    // transp: bool = false;
}

copy :: (m: Matrix64) -> Matrix64 {
    res := mat64(m.rows, m.cols);
    array_copy(*res.data, m.data);
    // res.transp = m.transp;
    return res;
}


operator [] :: inline (m: Matrix64, i: int) -> []float64 {
    #if CHECKS {
        assert(i >= 0, "Row-index negative! i = %", i);
        assert(i < m.rows, "Row-index out of bounds: % rows", m.rows);
    }

    return array_view_count(m.data, i * m.cols, m.cols);
}

row :: (m: Matrix64, i: int) -> Vector64 {
    res := vec64(m.cols);
    for 0..m.cols-1 {
        res.data[it] = m[i][it];
    }
    return res;
}

operator == :: inline (a: Matrix64, b: Matrix64) -> bool {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    r := true;
    for a.data {
        r = r & (it == b.data[it_index]);
    }
    return a.rows == b.rows && a.cols == b.cols && r;
}


operator + :: (a: Matrix64, b: Matrix64) -> Matrix64 {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] + b.data[it];
    }
    return res;
}

operator - :: (a: Matrix64, b: Matrix64) -> Matrix64 {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}

operator - :: (a: Matrix64) -> Matrix64 {
    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = -a.data[it];
    }
    return res;
}

operator * :: (a: Matrix64, val: float64) -> Matrix64 #symmetric {
    res := mat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = val * it;
    }
    return res;
}

operator / :: (a: Matrix64, val: float64) -> Matrix64 {
    res := mat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = it / val;
    }
    return res;
}


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 2.11 Case Study: Matrix-Matrix Multiplication, page 118
//
// According to Table 2.12, the best performing loop order is j/k/i with
// the least number of cache misses.
//
// This is still a very naive implementation and there are known algorithms
// that have much better perfomances! 
//
// WARNING: no dimension checks! Make sure you  know what you're doing.
operator * :: (a: Matrix64, b: Matrix64) -> Matrix64 {
    #if CHECKS {
        assert(a.cols == b.rows, "Matrix dimensions don't match: a.cols = %, b.rows = %", a.cols, b.rows);
    }

    res := mat64_zeros(a.rows, b.cols);
    for j: 0..b.cols-1 {
        for k: 0..a.cols-1 {
            for i: 0..a.rows-1 {
                res[i][j] = res[i][j] + a[i][k] * b[k][j];
                // set(res, i, j, get(res, i, j) + get(a, i, k) * get(b, k, j));
            }
        }
    }
    return res;
}






operator * :: (a: Matrix64, b: Vector64) -> Vector64 {
    #if CHECKS {
        assert(a.cols == b.dim, "Dimensions mismatch: a.cols = %, b.dim = %", a.cols, b.dim);
    }

    res := vec64_zeros(a.rows);
    for i: 0..a.rows-1 {
        for j: 0..a.cols-1 {
            // res.data[i] += get(a, i, j) * b[j];
            res.data[i] += a[i][j] * b[j];
        }
    }
    return res;
}

operator * :: (a: Vector64, b: Matrix64) -> Vector64 {
    #if CHECKS {
        assert(b.cols == a.dim, "Dimensions mismatch: b.cols = %, a.dim = %", b.cols, a.dim);
    }

    res := vec64_zeros(a.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.rows-1 {
            res.data[i] += b[i][j] * a[j];
            // res.data[i] += get(b, i, j) * a[j];
        }
    }
    return res;
}

operator * :: (a: Matrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(a.cols == b.dim, "Dimensions mismatch: a.cols = %, b.dim = %", a.cols, b.dim);
    }

    res := cvec64_zeros(a.rows);
    for i: 0..a.rows-1 {
        for j: 0..a.cols-1 {
            // res.data[i] += get(a, i, j) * b[j];
            res.data[i] += a[i][j] * b[j];
        }
    }
    return res;
}

operator * :: (a: CVector64, b: Matrix64) -> CVector64 {
    #if CHECKS {
        assert(b.cols == a.dim, "Dimensions mismatch: b.cols = %, a.dim = %", b.cols, a.dim);
    }

    res := cvec64_zeros(a.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.rows-1 {
            // res.data[i] += get(b, i, j) * a[j];
            res.data[i] += b[i][j] * a[j];
        }
    }
    return res;
}








mat64 :: (rows: int, cols: int) -> Matrix64 {
    res : Matrix64;
    res.rows = rows;
    res.cols = cols;
    array_resize(*res.data, rows * cols, false);
    return res;
}

mat64 :: inline (dim: int) -> Matrix64 {
    return mat64(dim, dim);
}

mat64 :: (rows: int, cols: int, init: float64) -> Matrix64 {
    m := mat64(rows, cols);
    for 0..rows*cols-1 {
        m.data[it] = init;
    }
    return m;
}

mat64 :: (rows: int, cols:int, values: ..float64) -> Matrix64 {
    #if CHECKS {
        assert(values.count <= rows * cols, "Too many values: rows*cols = %, values.count = %", rows*cols, values.count);
    }

    m := mat64_zeros(rows, cols);
    for values {
        m.data[it_index] = it;
    }
    return m;
}


mat64 :: inline (dim: int, init: float64) -> Matrix64 {
    return mat64(dim, dim, init);
}





mat64_zeros :: inline (rows: int, cols: int) -> Matrix64 {
    return mat64(rows, cols, 0.0);
}

mat64_zeros :: inline (dim: int) -> Matrix64 {
    return mat64(dim, dim, 0.0);
}





mat64_ones :: inline (rows: int, cols: int) -> Matrix64 {
    return mat64(rows, cols, 1.0);
}

mat64_ones :: inline (dim: int) -> Matrix64 {
    return mat64(dim, dim, 1.0);
}


get :: inline (m: Matrix64, row: int, col: int) -> float64 {
    #if CHECKS {
        assert(row <= m.rows && col <= m.cols,
            "Out of bounds: m: %x%, row = %, col = %", m.rows, m.cols, row, col);
    }

    // if m.transp {
        // return m.data[col * m.cols + row];
    // }
    return m.data[row * m.cols + col];
}

set :: inline (m: Matrix64, row: int, col: int, val: float64) {
    #if CHECKS {
        assert(row <= m.rows && col <= m.cols,
            "Out of bounds: m: %x%, row = %, col = %", m.rows, m.cols, row, col);
    }

    // if m.transp {
    //     m.data[col * m.cols + row] = val;
    // } else {
    m.data[row * m.cols + col] = val;
    // }
}


mat64_id :: (rows: int, cols: int) -> Matrix64 {
    m := mat64_zeros(rows, cols);
    dim := ifx rows <= cols then rows else cols;
    for 0..dim-1 {
        set(m, it, it, 1.0);
    }
    return m;
}

mat64_id :: inline (dim: int) -> Matrix64 {
    return mat64_id(dim, dim);
}





mat64_hadamard_2 :: inline () -> Matrix64 {
    return mat64(2, 2, 1.0, 1.0, 1.0, -1.0);
}

mat64_hadamard :: inline (nestings: int) -> Matrix64 {
    #if CHECKS {
        assert(nestings >= 0, "Nestings must be >=0 but = %", nestings);
    }
    if nestings == 0 {
        return mat64_hadamard_2();
    }
    return tensor(mat64_hadamard_2(), mat64_hadamard(nestings-1));
}




// https://www.wikiwand.com/en/Tensor_product
tensor :: (a: Matrix64, b: Matrix64) -> Matrix64 {
    res := mat64(a.rows * b.rows, a.cols * b.cols);
    for r: 0..res.rows-1 {
        ar := r/b.rows;
        br := r % b.rows;
        
        for c: 0..res.cols-1 {
            ac := c/b.cols;
            bc := c % b.cols;
            res[r][c] = a[ar][ac] * b[br][bc];
            // set(res, r, c, get(a, ar, ac) * get(b, br, bc));
        }
    }
    return res;
}









transpose :: (m: Matrix64) -> Matrix64 #must {
    // res := copy(m);
    // res.transp = !m.transp;
    // return res;

    res := mat64(m.cols, m.rows);
    for r: 0..m.rows-1 {
        for c: 0..m.cols-1 {
            res[c][r] = m[r][c];
            // set(res, c, r, get(m, r, c));
        }
    }
    return res;
}

// transpose :: (m: *Matrix64) {
//     m.transp

//     // for r: 0..m.rows-1 {
//     //     for c: 0..m.cols-1 {
//     //         i1 := r * m.cols + c;
//     //         i2 := 
//     //         set(<<m, c, r, get(<<m, r, c));
//     //     }
//     // }
// }






mat64_reflector :: (normal: Vector64) -> Matrix64 {
    tau := 2.0/(normal*normal);
    return mat64_id(normal.dim) - tau * normal / normal;
}








test_matrix :: () {
    print("\n\nMatrix -------------------------------\n\n");
    {
        m := mat64(3, 3);
        print("m = %\n", m);
    }
    {
        a := 3;
        m := mat64(a, a);
        print("m = %\n", m);
    }



    {
        m := mat64(3, 3, 42.0);
        print("m = %\n", m);
    }
    {
        a := 3;
        m := mat64(a, a, 42.0);
        print("m = %\n", m);
    }


    {
        m := mat64(3, 3, 42.0);
        set(m, 1, 2, 1.0);
        assert(get(m, 1, 2) == 1.0);
        print("m = %\n", m);
    }

    {
        m := mat64_id(4);
        print("m = %\n", m);
        assert(get(m, 0, 1) == 0.0);
        assert(get(m, 1, 1) == 1.0);
        assert(get(m, 2, 1) == 0.0);
        assert(get(m, 3, 3) == 1.0);
    }

    {
        m := mat64_zeros(3);
        set(m, 0, 2, 1.0);
        m = transpose(m);
        assert(get(m, 2, 0) == 1.0);
    }

    {
        m := mat64(2, 2, 1.0, 2.0, 3.0, 4.0);
        print("m = %\n", m);
        m = m * m;
        print("m = %\n", m);
    }
    {
        m := mat64(1, 3, 1.0, 2.0, 3.0);
        print("m = %\n", m);
        mt := transpose(m);
        m = m * mt;
        print("m = %\n", m);
    }

    {
        m := mat64(2, 2, 1.0, 2.0, 3.0, 4.0);
        t := tensor(m, m);
        assert(t == mat64(4, 4,
            1.0,  2.0,  2.0,  4.0,
            3.0,  4.0,  6.0,  8.0,
            3.0,  6.0,  4.0,  8.0,
            9.0, 12.0, 12.0, 16.0
        ));
        print("t = %\n", t);
    }
    {
        assert(mat64_hadamard_2() == mat64(2,2, 1.0, 1.0, 1.0, -1.0));
    }
    {
        assert(mat64_hadamard(1) == mat64(4,4,
            1, 1, 1, 1,
            1,-1, 1,-1,
            1, 1,-1,-1,
            1,-1,-1, 1
        ));
    }
}