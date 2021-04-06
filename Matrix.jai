







str :: (m: Matrix64, indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    append(*builder, "\n");
    for i: 0..m.rows-1 {
        append(*builder, indent);
        for j: 0..m.cols-1 {
            print_to_builder(*builder, "%, ", m[i][j] );
        }
        append(*builder, "\n");
    }
    return builder_to_string(*builder);
}







Matrix64 :: struct {
    rows: int;
    cols: int;
    data: []float64;
}

mat64 :: (rows: int, cols: int) -> Matrix64 #must {
    res : Matrix64 = ---;
    res.rows = rows;
    res.cols = cols;
    res.data = NewArray(rows*cols, float64, initialized=false);
    return res;
}

copy :: (m: Matrix64) -> Matrix64 #must {
    res : Matrix64 = ---;
    res.cols = m.cols;
    res.rows = m.rows;
    res.data = array_copy(m.data);
    return res;
}

copy :: (to: *Matrix64, from: Matrix64) {
    to.rows = from.rows;
    to.cols = from.cols;
    to.data = array_copy(from.data);
}


operator [] :: inline (m: Matrix64, i: int) -> []float64 {
    #if CHECKS {
        assert(i >= 0, "Row-index negative! i = %", i);
        assert(i < m.rows, "Row-index out of bounds: % rows", m.rows);
    }

    return array_view_count(m.data, i * m.cols, m.cols);
}
row :: (m: Matrix64, i: int) -> Vector64 #must {
    #if CHECKS {
        assert(i >= 0 && i < m.rows);
    }
    res := vec64(m.cols);
    for 0..m.cols-1 {
        res.data[it] = m[i][it];
    }
    return res;
}

column :: (m: Matrix64, i: int) -> Vector64 #must {
    #if CHECKS {
        assert(i >= 0 && i < m.cols);
    }
    res := vec64(m.rows);
    for 0..m.rows-1 {
        res.data[it] = m[it][i];
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


operator + :: (a: Matrix64, b: Matrix64) -> Matrix64 #must {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] + b.data[it];
    }
    return res;
}
add :: (a: Matrix64, b: Matrix64) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    for 0..a.data.count-1 {
        a.data[it] += b.data[it];
    }
}

operator - :: (a: Matrix64, b: Matrix64) -> Matrix64 #must {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}
sub :: (a: Matrix64, b: Matrix64) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    for 0..a.data.count-1 {
        a.data[it] -= b.data[it];
    }
}

operator - :: (a: Matrix64) -> Matrix64 #must {
    res := mat64(a.rows, a.cols);
    for 0..a.data.count {
        res.data[it] = -a.data[it];
    }
    return res;
}
neg :: (a: Matrix64) {
    for a.data {
        a.data[it_index] = -it;
    }
}

operator * :: (a: Matrix64, val: float64) -> Matrix64 #must #symmetric {
    res := mat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = val * it;
    }
    return res;
}
mul :: (a: Matrix64, val: float64) {
    for 0..a.data.count-1 {
        a.data[it] *= val;
    }
}

operator / :: (a: Matrix64, val: float64) -> Matrix64 #must {
    res := mat64(a.rows, a.cols);
    for a.data {
        res.data[it_index] = it / val;
    }
    return res;
}
div :: (a: Matrix64, val: float64) {
    for 0..a.data.count-1 {
        a.data[it] /= val;
    }
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
operator * :: (a: Matrix64, b: Matrix64) -> Matrix64 #must {
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






operator * :: (a: Matrix64, b: Vector64) -> Vector64 #must {
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

operator * :: (a: Vector64, b: Matrix64) -> Vector64 #must {
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

operator * :: (a: Matrix64, b: CVector64) -> CVector64 #must {
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

operator * :: (a: CVector64, b: Matrix64) -> CVector64 #must {
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







operator / :: inline (a: Matrix64, b: Matrix64) -> Matrix64 #must {
    return inverse(b) * a;
}







mat64 :: inline (dim: int) -> Matrix64 {
    return mat64(dim, dim);
}

mat64 :: (rows: int, cols: int, init: float64) -> Matrix64 #must {
    m := mat64(rows, cols);
    for 0..rows*cols-1 {
        m.data[it] = init;
    }
    return m;
}

mat64 :: (rows: int, cols:int, values: ..float64) -> Matrix64 #must {
    #if CHECKS {
        assert(values.count <= rows * cols, "Too many values: rows*cols = %, values.count = %", rows*cols, values.count);
    }

    m := mat64_zeros(rows, cols);
    for values {
        m.data[it_index] = it;
    }
    return m;
}


mat64 :: inline (dim: int, init: float64) -> Matrix64 #must {
    return mat64(dim, dim, init);
}





mat64_zeros :: inline (rows: int, cols: int) -> Matrix64 #must {
    return mat64(rows, cols, 0.0);
}

mat64_zeros :: inline (dim: int) -> Matrix64 #must {
    return mat64(dim, dim, 0.0);
}





mat64_ones :: inline (rows: int, cols: int) -> Matrix64 #must {
    return mat64(rows, cols, 1.0);
}

mat64_ones :: inline (dim: int) -> Matrix64 #must {
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


mat64_id :: (rows: int, cols: int) -> Matrix64 #must {
    m := mat64_zeros(rows, cols);
    dim := ifx rows <= cols then rows else cols;
    for 0..dim-1 {
        set(m, it, it, 1.0);
    }
    return m;
}

mat64_id :: inline (dim: int) -> Matrix64 #must {
    return mat64_id(dim, dim);
}





mat64_hadamard_2 :: inline () -> Matrix64 #must {
    return mat64(2, 2, 1.0, 1.0, 1.0, -1.0);
}

mat64_hadamard :: inline (nestings: int) -> Matrix64 #must {
    #if CHECKS {
        assert(nestings >= 0, "Nestings must be >=0 but = %", nestings);
    }
    if nestings == 0 {
        return mat64_hadamard_2();
    }
    return tensor(mat64_hadamard_2(), mat64_hadamard(nestings-1));
}





mat64_reflector :: (normal: Vector64) -> Matrix64 #must {

    tau := 2.0/(normal*normal);

    res := mat64_id(normal.dim);
    sub(res, tau * normal / normal);
    return res;
}


// inclusive first bound, exclusive last bound
mat64_submatrix :: (M: Matrix64, row_from: int, row_to: int, column_from: int, column_to: int) -> Matrix64 #must {
    #if CHECKS {
        assert(row_to > row_from, "\ncannot create submatrix with rows %->%\n", row_to, row_from);
        assert(column_to > column_from, "\ncannot create submatrix with columns %->%\n", column_to, column_from);
    }
    res := mat64(row_to-row_from, column_to-column_from);

    for i: 0..row_to-row_from-1 {
        for j: 0..column_to-column_from-1 {
            res[i][j] = M[i+row_from][j+column_from];
        }
    }

    return res;
}










// there are faster ways of doing this.
permute_rows :: (M: Matrix64, order: ..int) -> Matrix64 #must {
    #if CHECKS {
        assert(order.count == M.rows);
    }

    res := mat64(M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[order[i]][j];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute_columns :: (M: Matrix64, order: ..int) -> Matrix64 #must {
    #if CHECKS {
        assert(order.count == M.cols);
    }

    res := mat64(M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[i][order[j]];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute :: (M: Matrix64, row_order: []int, col_order: []int) -> Matrix64 #must {
    #if CHECKS {
        assert(row_order.count == M.rows);
        assert(col_order.count == M.cols);
    }

    res := mat64(M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[row_order[i]][col_order[j]];
        }
    }

    return res;
}


swap_columns :: (M: Matrix64, i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0 && i != j);
        assert(i < M.cols && j < M.cols);
    }

    for 0..M.rows-1 {
        M[it][i], M[it][j] = swap(M[it][i], M[it][j]);
    }
}

swap_rows :: (M: Matrix64, i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0 && i != j);
        assert(i < M.rows && j < M.rows);
    }

    for 0..M.cols-1 {
        M[i][it], M[j][it] = swap(M[i][it], M[j][it]);
    }
}








// https://www.wikiwand.com/en/Tensor_product
tensor :: (a: Matrix64, b: Matrix64) -> Matrix64 #must {
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

transpose :: (m: *Matrix64) {
    <<m = transpose(<<m);
}



conjugate :: (m: Matrix64) -> Matrix64 #must {
    return copy(m);
}

conjugate :: (m: *Matrix64) {}


conjugate_transpose :: inline (m: Matrix64) -> Matrix64 #must {
    return transpose(m);
}

conjugate_transpose :: inline (m: *Matrix64) {
    transpose(m);
}






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Lemma 3.5.4 Matrix Norms, page 176

// this needs the largest eigenvalue
// norm_2 :: inline (v: CVector64) -> float64 {
//     return norm(v, 2.0);
// }

norm_1 :: (m: Matrix64) -> float64 {
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

norm_inf :: (m: Matrix64) -> float64 {
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

norm_frobenius :: (m: Matrix64) -> float64 {
    sum : float64 = 0.0;
    for m.data {
        sum += abs(it);
    }
    return sum;
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