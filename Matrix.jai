




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


operator + :: (a: Matrix64, b: Matrix64) -> Matrix64 {
    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] + b.data[it];
    }
    return res;
}

operator - :: inline (a: Matrix64, b: Matrix64) -> Matrix64 {
    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}

operator - :: inline (a: Matrix64) -> Matrix64 {
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
// According to Table 2.12, the best performing loop order is j/k/i with
// the least number of cache misses.
// WARNING: no dimension checks! Make sure you  know what you're doing.
operator * :: (a: Matrix64, b: Matrix64) -> Matrix64 {
    res := mat64_zeros(a.rows, b.cols);
    for j: 0..b.cols-1 {
        for k: 0..a.cols-1 {
            for i: 0..a.rows-1 {
                set(res, i, j, get(res, i, j) + get(a, i, k) * get(b, k, j));
            }
        }
    }
    return res;
}






operator * :: (a: Matrix64, b: Vector64) -> Vector64 {
    res := vec64_zeros(a.rows);
    for i: a.rows {
        for j: a.cols {
            res[i] += get(a, i, j) * b[j];
        }
    }
    return res;
}

operator * :: (a: Vector64, b: Matrix64) -> Vector64 {
    res := vec64_zeros(a.cols);
    for i: a.cols {
        for j: a.rows {
            res[i] += get(a, i, j) * b[j];
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

mat64 :: inline (rows: int, cols: int, init: float64) -> Matrix64 {
    m := mat64(rows, cols);
    for 0..rows*cols-1 {
        m.data[it] = init;
    }
    return m;
}

mat64 :: inline(rows: int, cols:int, values: ..float64) -> Matrix64 {
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




operator [] :: inline (m: Matrix64, i: int) -> float64 {
    return m.data[i];
}

get :: inline (m: Matrix64, row: int, col: int) -> float64 {
    // if m.transp {
        // return m.data[col * m.cols + row];
    // }
    return m.data[row * m.cols + col];
}

set :: inline (m: Matrix64, row: int, col: int, val: float64) {
    // if m.transp {
    //     m.data[col * m.cols + row] = val;
    // } else {
    m.data[row * m.cols + col] = val;
    // }
}


mat64_id :: inline (rows: int, cols: int) -> Matrix64 {
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






transpose :: (m: Matrix64) -> Matrix64 #must {
    // res := copy(m);
    // res.transp = !m.transp;
    // return res;

    res := mat64(m.cols, m.rows);
    for r: 0..m.rows-1 {
        for c: 0..m.cols-1 {
            set(res, c, r, get(m, r, c));
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
}