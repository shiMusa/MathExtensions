




CMatrix64 :: struct {
    rows: int;
    cols: int;
    data: [..]Complex64;
    // transp: bool = false;
}

copy :: (m: CMatrix64) -> CMatrix64 {
    res := cmat64(m.rows, m.cols);
    array_copy(*res.data, m.data);
    // res.transp = m.transp;
    return res;
}


operator == :: inline (a: CMatrix64, b: CMatrix64) -> bool {
    r := true;
    for a.data {
        r = r & (it == b.data[it_index]);
    }
    return a.rows == b.rows && a.cols == b.cols && r;
}


operator + :: (a: CMatrix64, b: CMatrix64) -> CMatrix64 {
    res := cmat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] + b.data[it];
    }
    return res;
}

operator - :: (a: CMatrix64, b: CMatrix64) -> CMatrix64 {
    res := cmat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}

operator - :: (a: CMatrix64) -> CMatrix64 {
    res := cmat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = -a.data[it];
    }
    return res;
}



operator * :: (a: CMatrix64, val: float64) -> CMatrix64 #symmetric {
    res := cmat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = val * it;
    }
    return res;
}

operator * :: (a: CMatrix64, val: Complex64) -> CMatrix64 #symmetric {
    res := cmat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = val * it;
    }
    return res;
}



operator / :: (a: CMatrix64, val: float64) -> CMatrix64 {
    res := cmat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = it / val;
    }
    return res;
}

operator / :: (a: CMatrix64, val: Complex64) -> CMatrix64 {
    res := cmat64(a.rows, a.cols);
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
operator * :: (a: CMatrix64, b: CMatrix64) -> CMatrix64 {
    res := cmat64_zeros(a.rows, b.cols);
    for j: 0..b.cols-1 {
        for k: 0..a.cols-1 {
            for i: 0..a.rows-1 {
                set(res, i, j, get(res, i, j) + get(a, i, k) * get(b, k, j));
            }
        }
    }
    return res;
}






operator * :: (a: CMatrix64, b: CVector64) -> CVector64 {
    res := cvec64_zeros(a.rows);
    for i: 0..a.rows-1 {
        for j: 0..a.cols-1 {
            res.data[i] += get(a, i, j) * b[j];
        }
    }
    return res;
}

operator * :: (a: CVector64, b: CMatrix64) -> CVector64 {
    res := cvec64_zeros(a.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.rows-1 {
            res.data[i] += get(b, i, j) * b[j];
        }
    }
    return res;
}








cmat64 :: (rows: int, cols: int) -> CMatrix64 {
    res : CMatrix64;
    res.rows = rows;
    res.cols = cols;
    array_resize(*res.data, rows * cols, false);
    return res;
}

cmat64 :: inline (dim: int) -> CMatrix64 {
    return cmat64(dim, dim);
}

cmat64 :: inline (rows: int, cols: int, init: float64) -> CMatrix64 {
    return cmat64(rows, cols, c64(init, 0.0));
}
cmat64 ::  (rows: int, cols: int, init: Complex64) -> CMatrix64 {
    m := cmat64(rows, cols);
    for 0..rows*cols-1 {
        m.data[it] = init;
    }
    return m;
}



// TODO: wait for compiler bug fix
// cmat64 :: inline(rows: int, cols: int, values: ..float64) -> CMatrix64 {
//     m := cmat64_zeros(rows, cols);
//     for values {
//         m.data[it_index] = c64(it, 0.0);
//     }
//     return m;
// }

cmat64 :: (rows: int, cols: int, values: ..Complex64) -> CMatrix64 {
    m := cmat64_zeros(rows, cols);
    for values {
        m.data[it_index] = it;
    }
    return m;
}


cmat64 :: inline (dim: int, init: float64) -> CMatrix64 {
    return cmat64(dim, dim, c64(init, 0.0));
}
cmat64 :: inline (dim: int, init: Complex64) -> CMatrix64 {
    return cmat64(dim, dim, init);
}





cmat64_zeros :: inline (rows: int, cols: int) -> CMatrix64 {
    return cmat64(rows, cols, 0.0);
}

cmat64_zeros :: inline (dim: int) -> CMatrix64 {
    return cmat64(dim, dim, 0.0);
}





cmat64_ones :: inline (rows: int, cols: int) -> CMatrix64 {
    return cmat64(rows, cols, R);
}

cmat64_ones :: inline (dim: int) -> CMatrix64 {
    return cmat64(dim, dim, R);
}




operator [] :: inline (m: CMatrix64, i: int) -> Complex64 {
    return m.data[i];
}

get :: inline (m: CMatrix64, row: int, col: int) -> Complex64 {
    return m.data[row * m.cols + col];
}

set :: inline (m: CMatrix64, row: int, col: int, val: float64) {
    m.data[row * m.cols + col] = c64(val, 0.0);
}
set :: inline (m: CMatrix64, row: int, col: int, val: Complex64) {
    m.data[row * m.cols + col] = val;
}


cmat64_id :: (rows: int, cols: int) -> CMatrix64 {
    m := cmat64_zeros(rows, cols);
    dim := ifx rows <= cols then rows else cols;
    for 0..dim-1 {
        set(m, it, it, R);
    }
    return m;
}

cmat64_id :: inline (dim: int) -> CMatrix64 {
    return cmat64_id(dim, dim);
}






transpose :: (m: CMatrix64) -> CMatrix64 #must {
    // res := copy(m);
    // res.transp = !m.transp;
    // return res;

    res := cmat64(m.cols, m.rows);
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



conjugate :: (m: CMatrix64) -> CMatrix64 {
    res := cmat64(m.cols, m.rows);
    for m.data {
        res.data[it_index] = conjugate(it);
    }
    return res;
}


conjugate_transpose :: inline (m: CMatrix64) -> CMatrix64 {
    return conjugate(transpose(m));
}
dagger :: conjugate_transpose;






test_complex_matrix :: () {
    print("\n\nComplex Matrix -------------------------------\n\n");
    {
        m := cmat64_zeros(3, 3);
        print("m = %\n", m);
    }
    {
        a := 3;
        m := cmat64_zeros(a, a);
        print("m = %\n", m);
    }



    {
        m := cmat64(3, 3, 42.0);
        print("m = %\n", m);
    }
    {
        a := 3;
        m := cmat64(a, a, 42.0);
        print("m = %\n", m);
    }


    {
        m := cmat64(3, 3, 42.0);
        set(m, 1, 2, 1.0);
        assert(get(m, 1, 2) == c64(1.0, 0.0));
        print("m = %\n", m);
    }

    {
        m := cmat64_id(4);
        print("m = %\n", m);
        assert(get(m, 0, 1) == c64(0.0, 0.0));
        assert(get(m, 1, 1) == c64(1.0, 0.0));
        assert(get(m, 2, 1) == c64(0.0, 0.0));
        assert(get(m, 3, 3) == c64(1.0, 0.0));
    }

    {
        m := cmat64_zeros(3);
        set(m, 0, 2, 1.0);
        m = transpose(m);
        assert(get(m, 2, 0) == c64(1.0, 0.0));
    }

    {
        m := cmat64(2, 2, R, c64(2.0), c64(3.0), c64(4.0));
        print("m = %\n", m);
        m = m * m;
        print("m = %\n", m);
    }
    {
        m := cmat64(2, 2, R, 2.0*R, 3.0*R, 4.0*R);
        print("m = %\n", m);
        m = m * m;
        print("m = %\n", m);
    }
    // {
    //     m := cmat64(1, 3, 1.0, 2.0, 3.0);
    //     print("m = %\n", m);
    //     mt := transpose(m);
    //     m = m * mt;
    //     print("m = %\n", m);
    // }
    {
        m := cmat64(2, 2, 
            CZERO, -I, 
            I, CZERO
        );
        assert( dagger(m) == m );
    }
}