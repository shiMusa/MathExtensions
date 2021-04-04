








str :: (m: CMatrix64, indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    append(*builder, "\n");
    for i: 0..m.rows-1 {
        append(*builder, indent);
        for j: 0..m.cols-1 {
            print_to_builder(*builder, "%, ", str(m[i][j]) );
        }
        append(*builder, "\n");
    }
    return builder_to_string(*builder);
}







CMatrix64 :: struct {
    rows: int;
    cols: int;
    data: []Complex64;
}

cmat64 :: (rows: int, cols: int) -> CMatrix64 #must {
    res : CMatrix64 = ---;
    res.rows = rows;
    res.cols = cols;
    res.data = NewArray(rows*cols, Complex64, initialized=false);
    return res;
}

copy :: (m: CMatrix64) -> CMatrix64 #must {
    res := cmat64(m.rows, m.cols);
    res.data = array_copy(m.data);
    return res;
}

copy :: (to: *CMatrix64, from: CMatrix64) {
    to.rows = from.rows;
    to.cols = from.cols;
    to.data = array_copy(from.data);
}

operator [] :: inline (m: CMatrix64, i: int) -> []Complex64 {
    #if CHECKS {
        assert(i >= 0, "Row-index negative! i = %", i);
        assert(i < m.rows, "Row-index out of bounds: % rows", m.rows);
    }

    return array_view_count(m.data, i * m.cols, m.cols);
}

row :: (m: CMatrix64, i: int) -> CVector64 #must {
    #if CHECKS {
        assert(i >= 0 && i < m.rows);
    }
    res := cvec64(m.cols);
    for 0..m.cols-1 {
        res.data[it] = m[i][it];
    }
    return res;
}

column :: (m: CMatrix64, i: int) -> CVector64 #must {
    #if CHECKS {
        assert(i >= 0 && i < m.cols);
    }
    res := cvec64(m.rows);
    for 0..m.rows-1 {
        res.data[it] = m[it][i];
    }
    return res;
}

operator == :: inline (a: CMatrix64, b: CMatrix64) -> bool {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    r := true;
    for a.data {
        r = r & (it == b.data[it_index]);
    }
    return a.rows == b.rows && a.cols == b.cols && r;
}


operator + :: (a: CMatrix64, b: CMatrix64) -> CMatrix64 #must {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := cmat64(a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] + b.data[it];
    }
    return res;
}

add :: (a: *CMatrix64, b: CMatrix64) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    for 0..a.data.count-1 {
        a.data[it] += b.data[it];
    }
}

operator - :: (a: CMatrix64, b: CMatrix64) -> CMatrix64 #must {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := cmat64(a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}

sub :: (a: *CMatrix64, b: CMatrix64) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    for 0..a.data.count-1 {
        a.data[it] -= b.data[it];
    }
}

operator - :: (a: CMatrix64) -> CMatrix64 #must {
    res := cmat64(a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = -a.data[it];
    }
    return res;
}

neg :: (a: *CMatrix64) {
    for a.data {
        a.data[it_index] = -it;
    }
}



operator * :: (a: CMatrix64, val: float64) -> CMatrix64 #must #symmetric {
    res := cmat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = val * it;
    }
    return res;
}

operator * :: (a: CMatrix64, val: Complex64) -> CMatrix64 #must #symmetric {
    res := cmat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = val * it;
    }
    return res;
}

mul :: (a: *CMatrix64, val: float64) {
    for 0..a.data.count-1 {
        a.data[it] *= val;
    }
}
mul :: (a: *CMatrix64, val: Complex64) {
    for 0..a.data.count-1 {
        a.data[it] *= val;
    }
}


operator / :: (a: CMatrix64, val: float64) -> CMatrix64 #must {
    res := cmat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = it / val;
    }
    return res;
}

operator / :: (a: CMatrix64, val: Complex64) -> CMatrix64 #must {
    res := cmat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = it / val;
    }
    return res;
}

div :: (a: *CMatrix64, val: float64) {
    for 0..a.data.count-1 {
        a.data[it] /= val;
    }
}
div :: (a: *CMatrix64, val: Complex64) {
    for 0..a.data.count-1 {
        a.data[it] /= val;
    }
}


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 2.11 Case Study: Matrix-Matrix Multiplication, page 118
// According to Table 2.12, the best performing loop order is j/k/i with
// the least number of cache misses.
// WARNING: no dimension checks! Make sure you  know what you're doing.
operator * :: (a: CMatrix64, b: CMatrix64) -> CMatrix64 #must {
    #if CHECKS {
        assert(a.cols == b.rows, "Matrix dimensions don't match: a.cols = %, b.rows = %", a.cols, b.rows);
    }

    res := cmat64_zeros(a.rows, b.cols);
    for j: 0..b.cols-1 {
        for k: 0..a.cols-1 {
            tmp := b[k][j];
            for i: 0..a.rows-1 {
                res[i][j] += a[i][k] * tmp;
            }
        }
    }
    return res;
}


// TODO: think of a smart way to do in-place multiplication
//
// mul :: (a: CMatrix64, b: CMatrix64) {
//     #if CHECKS {
//         assert(a.cols == b.rows, "Matrix dimensions don't match: a.cols = %, b.rows = %", a.cols, b.rows);
//     }

//     res := cmat64_zeros(a.rows, b.cols);
//     for j: 0..b.cols-1 {
//         for k: 0..a.cols-1 {
//             tmp := b[k][j];
//             for i: 0..a.rows-1 {
//                 res[i][j] += a[i][k] * tmp;
//             }
//         }
//     }
//     return res;
// }




operator * :: (a: CMatrix64, b: CVector64) -> CVector64 #must {
    #if CHECKS {
        assert(a.cols == b.dim, "Dimensions mismatch: a.cols = %, b.dim = %", a.cols, b.dim);
    }

    res := cvec64_zeros(a.rows);
    for i: 0..a.rows-1 {
        for j: 0..a.cols-1 {
            res.data[i] += a[i][j] * b[j];
        }
    }
    return res;
}

operator * :: (a: CVector64, b: CMatrix64) -> CVector64 #must {
    #if CHECKS {
        assert(b.cols == a.dim, "Dimensions mismatch: b.cols = %, a.dim = %", b.cols, a.dim);
    }

    res := cvec64_zeros(a.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.rows-1 {
            res.data[i] += b[i][j] * a[j];
        }
    }
    return res;
}

operator * :: (a: CMatrix64, b: Vector64) -> CVector64 #must {
    #if CHECKS {
        assert(a.cols == b.dim, "Dimensions mismatch: a.cols = %, b.dim = %", a.cols, b.dim);
    }

    res := cvec64_zeros(a.rows);
    for i: 0..a.rows-1 {
        for j: 0..a.cols-1 {
            res.data[i] += a[i][j] * b[j];
        }
    }
    return res;
}

operator * :: (a: Vector64, b: CMatrix64) -> CVector64 #must {
    #if CHECKS {
        assert(b.cols == a.dim, "Dimensions mismatch: b.cols = %, a.dim = %", b.cols, a.dim);
    }
    
    res := cvec64_zeros(a.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.rows-1 {
            res.data[i] += b[i][j] * a[j];
        }
    }
    return res;
}






cmat64 :: (m: Matrix64) -> CMatrix64 #must {
    res := cmat64(m.rows, m.cols);
    for 0..res.data.count-1 {
        res.data[it] = c64(m.data[it]);
    }
    return res;
}

cmat64 :: inline (dim: int) -> CMatrix64 #must {
    return cmat64(dim, dim);
}

cmat64 :: inline (rows: int, cols: int, init: float64) -> CMatrix64 #must {
    return cmat64(rows, cols, c64(init, 0.0));
}
cmat64 ::  (rows: int, cols: int, init: Complex64) -> CMatrix64 #must {
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

cmat64 :: (rows: int, cols: int, values: ..Complex64) -> CMatrix64 #must {
    #if CHECKS {
        assert(values.count <= rows * cols, "Too many values: rows*cols = %, values.count = %", rows*cols, values.count);
    }

    m := cmat64_zeros(rows, cols);
    for values {
        m.data[it_index] = it;
    }
    return m;
}


cmat64 :: inline (dim: int, init: float64) -> CMatrix64 #must {
    return cmat64(dim, dim, c64(init, 0.0));
}
cmat64 :: inline (dim: int, init: Complex64) -> CMatrix64 #must {
    return cmat64(dim, dim, init);
}





cmat64_zeros :: inline (rows: int, cols: int) -> CMatrix64 #must {
    return cmat64(rows, cols, 0.0);
}

cmat64_zeros :: inline (dim: int) -> CMatrix64 #must {
    return cmat64(dim, dim, 0.0);
}





cmat64_ones :: inline (rows: int, cols: int) -> CMatrix64 #must {
    return cmat64(rows, cols, R);
}

cmat64_ones :: inline (dim: int) -> CMatrix64 #must {
    return cmat64(dim, dim, R);
}



cmat64_id :: (rows: int, cols: int) -> CMatrix64 #must {
    m := cmat64_zeros(rows, cols);
    dim := ifx rows <= cols then rows else cols;
    for 0..dim-1 {
        m[it][it] = R;
    }
    return m;
}

cmat64_id :: inline (dim: int) -> CMatrix64 #must {
    return cmat64_id(dim, dim);
}





cmat64_reflector :: (normal: CVector64) -> CMatrix64 #must {

    tau := c64(2.0)/(conjugate(normal)*normal);

    res := cmat64_id(normal.dim);
    sub(res, tau * normal / normal);
    return res;
}


// inclusive first bound, exclusive last bound
cmat64_submatrix :: (M: CMatrix64, row_from: int, row_to: int, column_from: int, column_to: int) -> CMatrix64 {
    #if CHECKS {
        assert(row_to > row_from, "\ncannot create submatrix with rows %->%\n", row_to, row_from);
        assert(column_to > column_from, "\ncannot create submatrix with columns %->%\n", column_to, column_from);
    }
    res := cmat64(row_to-row_from, column_to-column_from);

    for i: 0..row_to-row_from-1 {
        for j: 0..column_to-column_from-1 {
            res[i][j] = M[i+row_from][j+column_from];
        }
    }

    return res;
}




// there are faster ways of doing this.
permute_rows :: (M: CMatrix64, order: ..int) -> CMatrix64 {
    #if CHECKS {
        assert(order.count == M.rows);
    }

    res := cmat64(M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[order[i]][j];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute_columns :: (M: CMatrix64, order: ..int) -> CMatrix64 {
    #if CHECKS {
        assert(order.count == M.cols);
    }

    res := cmat64(M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[i][order[j]];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute :: (M: CMatrix64, row_order: []int, col_order: []int) -> CMatrix64 {
    #if CHECKS {
        assert(row_order.count == M.rows);
        assert(col_order.count == M.cols);
    }

    res := cmat64(M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[row_order[i]][col_order[j]];
        }
    }

    return res;
}


swap_columns :: (M: CMatrix64, i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0 && i != j);
        assert(i < M.cols && j < M.cols);
    }

    for 0..M.rows-1 {
        M[it][i], M[it][j] = swap(M[it][i], M[it][j]);
    }
}

swap_rows :: (M: CMatrix64, i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0 && i != j);
        assert(i < M.rows && j < M.rows);
    }

    for 0..M.cols-1 {
        M[i][it], M[j][it] = swap(M[i][it], M[j][it]);
    }
}










transpose :: (m: CMatrix64) -> CMatrix64 #must {
    // res := copy(m);
    // res.transp = !m.transp;
    // return res;

    res := cmat64(m.cols, m.rows);
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
// hermitian :: conjugate_transpose;








// https://www.wikiwand.com/en/Tensor_product
tensor :: (a: CMatrix64, b: CMatrix64) -> CMatrix64 {
    res := cmat64(a.rows * b.rows, a.cols * b.cols);
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

// https://www.wikiwand.com/en/Tensor_product
tensor :: (a: CMatrix64, b: Matrix64) -> CMatrix64 {
    res := cmat64(a.rows * b.rows, a.cols * b.cols);
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

// https://www.wikiwand.com/en/Tensor_product
tensor :: (a: Matrix64, b: CMatrix64) -> CMatrix64 {
    res := cmat64(a.rows * b.rows, a.cols * b.cols);
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












// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Lemma 3.5.4 Matrix Norms, page 176

// this needs the largest eigenvalue
// norm_2 :: inline (v: CVector64) -> float64 {
//     return norm(v, 2.0);
// }

norm_1 :: (m: CMatrix64) -> float64 {
    res : float64 = 0.0;

    for j: 0..m.cols {
        sum : float64 = 0.0;
        for i: 0..m.rows {
            sum += abs(m[i][j]);
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}

norm_inf :: (m: CMatrix64) -> float64 {
    res : float64 = 0.0;

    for i: 0..m.rows {
        sum : float64 = 0.0;
        for j: 0..m.cols {
            sum += abs(m[i][j]);
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}

norm_frobenius :: (m: CMatrix64) -> float64 {
    sum : float64 = 0.0;
    for m.data {
        sum += abs(it);
    }
    return sum;
}


// condition_frobenius :: (m: CMatrix64) -> float64 {
// }










test_complex_matrix :: () {
    print("\n\nComplex Matrix -------------------------------\n\n");
    {
        m := cmat64_zeros(3, 3);
        print("%\n", str(m));
    }
    {
        a := 3;
        m := cmat64_zeros(a, a);
        print("%\n", str(m));
    }



    {
        m := cmat64(3, 3, 42.0);
        print("%\n", str(m));
    }
    {
        a := 3;
        m := cmat64(a, a, 42.0);
        print("%\n", str(m));
    }


    {
        m := cmat64(3, 3, 42.0);
        m[1][2] = R;
        assert(m[1][2] == c64(1.0, 0.0));
        print("m = \n%\n", str(m));
    }

    {
        m := cmat64_id(4);
        print("m = \n%\n", str(m));
        assert(m[0][1] == c64(0.0, 0.0));
        assert(m[1][1] == c64(1.0, 0.0));
        assert(m[2][1] == c64(0.0, 0.0));
        assert(m[3][3] == c64(1.0, 0.0));
    }

    {
        m := cmat64_zeros(3);
        m[0][2] = R;
        m = transpose(m);
        assert(m[2][0] == c64(1.0, 0.0));
    }

    {
        a := cmat64_id(2,2);
        b := cmat64_id(2,2);
        add(a, b);
        assert(a == 2.0 * b);

        a = a + b;
        assert(a == 3.0 * b);
    }

    {
        m := cmat64(2, 2, R, c64(2.0), c64(3.0), c64(4.0));
        print("m = \n%\n", str(m));
        m = m * m;
        print("m = \n%\n", str(m));
    }
    {
        m := cmat64(2, 2, R, 2.0*R, 3.0*R, 4.0*R);
        print("m = \n%\n", str(m));
        m = m * m;
        print("m = \n%\n", str(m));
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


    {
        m := cmat64(2, 2, R, I, c64(2.0, -1.0), 2.0 + 3.0*I);
        print("m = \n%\n", str(m));
        t1 := tensor(m, m);
        print("t1 = \n%\n", str(t1));
        t2 := cmat64(4, 4,
            R            ,             I,             I,   -R           ,
            2.0 -       I, 2.0 + 3.0 * I, 1.0 + 2.0 * I, -3.0 +  2.0 * I,
            2.0 -       I, 1.0 + 2.0 * I, 2.0 + 3.0 * I, -3.0 +  2.0 * I,
            3.0 - 4.0 * I, 7.0 + 4.0 * I, 7.0 + 4.0 * I, -5.0 + 12.0 * I
        );
        print("t2 = \n%\n", str(t2));
        assert(t1 == t2);
    }

    {
        M := cmat64(3,3, 
            CZERO, R, 2*R,
            3*R, 4*R, 5*R,
            6*R, 7*R, 8*R
        );

        S := cmat64_submatrix(M, 1, 3, 1, 3);
        assert(S == cmat64(2,2, 4*R, 5*R, 7*R, 8*R));
    }
}