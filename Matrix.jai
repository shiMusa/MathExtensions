




Matrix64 :: struct {
    rows: int;
    cols: int;
    data: [..]float64;
}




operator + :: (a: Matrix64, b: Matrix64) -> Matrix64 {
    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] + b.data[it];
    }
}

operator - :: inline (a: Matrix64, b: Matrix64) -> Matrix64 {
    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] - b.data[it];
    }
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
    return m.data[row * m.cols + col];
}

set :: inline (m: Matrix64, row: int, col: int, val: float64) {
    m.data[row * m.cols + col] = val;
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
    res := mat64(m.cols, m.rows);
    for r: 0..m.rows-1 {
        for c: 0..m.cols-1 {
            set(res, c, r, get(m, r, c));
        }
    }
    return res;
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

}