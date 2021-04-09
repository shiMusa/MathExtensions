








str :: (m: Matrix($T), indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    append(*builder, "\n");
    for i: 0..m.rows-1 {
        append(*builder, indent);
        for j: 0..m.cols-1 {
            #if T == Complex64 {
                print_to_builder(*builder, "%, ", str(m[i][j]) );
            } else {
                print_to_builder(*builder, "%, ", m[i][j] );
            }
        }
        append(*builder, "\n");
    }
    return builder_to_string(*builder);
}







Matrix :: struct(T: Type) {
    rows: int;
    cols: int;
    data: []T;
}

operator [] :: inline (m: Matrix($T), i: int) -> []T {
    #if CHECKS {
        assert(i >= 0, "Row-index negative! i = %", i);
        assert(i < m.rows, "Row-index out of bounds: % rows", m.rows);
    }

    return array_view_count(m.data, i * m.cols, m.cols);
}





mat :: ($type: Type, rows: int, cols: int) -> Matrix(type) #must {
    res : Matrix(type) = ---;
    res.rows = rows;
    res.cols = cols;
    res.data = NewArray(rows*cols, type, initialized=false);
    return res;
}

copy :: (m: Matrix($T)) -> Matrix(T) #must {
    res := mat(T, m.rows, m.cols);
    res.data = array_copy(m.data);
    return res;
}

copy :: (to: *Matrix($T), from: Matrix(T)) {
    to.rows = from.rows;
    to.cols = from.cols;
    to.data = array_copy(from.data);
}


row :: (m: Matrix($T), i: int) -> Vector(T) #must {
    #if CHECKS {
        assert(i >= 0 && i < m.rows);
    }
    res := vec(T, m.cols);
    for 0..m.cols-1 {
        res.data[it] = m[i][it];
    }
    return res;
}

column :: (m: Matrix($T), i: int) -> Vector(T) #must {
    #if CHECKS {
        assert(i >= 0 && i < m.cols);
    }
    res := vec(T, m.rows);
    for 0..m.rows-1 {
        res.data[it] = m[it][i];
    }
    return res;
}

operator == :: inline (a: Matrix($T), b: Matrix(T)) -> bool {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    r := true;
    for a.data {
        r = r & (it == b.data[it_index]);
    }
    return a.rows == b.rows && a.cols == b.cols && r;
}








#insert #run specialize("_add", "operator +", COMBI_M_M_M);
_add :: (a: Matrix($A), b: Matrix($B), res: *Matrix($C)) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    <<res = mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] + b.data[it];
    }
}
add :: (a: *Matrix($A), b:Matrix($B)) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
        assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] += ncast(A, b.data[it]);
    }
}


#insert #run specialize("_sub", "operator -", COMBI_M_M_M);
_sub :: (a: Matrix($A), b: Matrix($B), res: *Matrix($C)) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    <<res = mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] - b.data[it];
    }
}
sub :: (a: *Matrix($A), b:Matrix($B)) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
        assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] -= ncast(A, b.data[it]);
    }
}








operator - :: (a: Matrix($T)) -> Matrix(T) #must {
    res := mat(T, a.rows, a.cols);

    for 0..a.data.count-1 {
        res.data[it] = -a.data[it];
    }

    return res;
}
neg :: (a: Matrix($T)) {
    for a.data {
        a.data[it_index] = -it;
    }
}








#insert #run specialize("_mulx", "operator *", COMBI_M_X_M, .["#symmetric"]);
_mulx :: (a: Matrix($A), b: $B, res: *Matrix($C)) {
    <<res = mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = ncast(C, a.data[it]) * ncast(C, b);
    }
}
mul :: (a: *Matrix($A), b: Complex64) {
    #if CHECKS {
        assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] *= ncast(A, b);
    }
}
mul :: (a: *Matrix($A), b: float64) {
    #if CHECKS {
        assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] *= ncast(A, b);
    }
}


#insert #run specialize("_divx", "operator /", COMBI_M_X_M);
_divx :: (a: Matrix($A), b: $B, res: *Matrix($C)) {
    <<res = mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = ncast(C, a.data[it]) / ncast(C, b);
    }
}
div :: (a: *Matrix($A), b: Complex64) {
    #if CHECKS {
        assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] /= ncast(A, b);
    }
}
div :: (a: *Matrix($A), b: float64) {
    #if CHECKS {
        assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] /= ncast(A, b);
    }
}









// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 2.11 Case Study: Matrix-Matrix Multiplication, page 118
// According to Table 2.12, the best performing loop order is j/k/i with
// the least number of cache misses.




#insert #run specialize("_mul", "operator *", COMBI_M_M_M);
_mul :: (a: Matrix($A), b: Matrix($B), res: *Matrix($C)) {
    #if CHECKS {
        assert(a.cols == b.rows, "Matrix dimensions don't match: a.cols = %, b.rows = %", a.cols, b.rows);
    }

    <<res = mat_zeros(C, a.rows, b.cols);
    for j: 0..b.cols-1 {
        for k: 0..a.cols-1 {
            tmp := ncast(C, b[k][j]);
            for i: 0..a.rows-1 {
                (<<res)[i][j] += ncast(C, a[i][k]) * tmp;
            }
        }
    }
}
// TODO: think of a smart way to do in-place multiplication



#insert #run specialize("_mul", "operator *", COMBI_M_V_V);
_mul :: (a: Matrix($A), b: Vector($B), res: *Vector($C)) {
    #if CHECKS {
        assert(a.cols == b.dim, "Dimensions mismatch: a.cols = %, b.dim = %", a.cols, b.dim);
    }

    <<res = vec_zeros(C, a.rows);
    for i: 0..a.rows-1 {
        for j: 0..a.cols-1 {
            res.data[i] += ncast(C, a[i][j]) * ncast(C, b[j]);
        }
    }
}
#insert #run specialize("_mul", "operator *", COMBI_V_M_V);
_mul :: (a: Vector($A), b: Matrix($B), res: *Vector($C)) {
    #if CHECKS {
        assert(b.cols == a.dim, "Dimensions mismatch: b.cols = %, a.dim = %", b.cols, a.dim);
    }

    <<res = vec_zeros(C, a.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.rows-1 {
            res.data[i] += ncast(C, b[i][j]) * ncast(C, a[j]);
        }
    }
}





#insert #run specialize("_div", "operator /", COMBI_M_M_M);
_div :: inline (a: Matrix($A), b: Matrix($B), res: *Matrix($C)) {
    <<res = inverse(b) * a;
}








mat :: inline ($type: Type, dim: int) -> Matrix(type) #must {
    return mat(type, dim, dim);
}

mat ::  (rows: int, cols: int, init: $T) -> Matrix(T) #must {
    m := mat(T, rows, cols);
    for 0..rows*cols-1 {
        m.data[it] = init;
    }
    return m;
}



// TODO: wait for compiler bug fix
// mat :: inline(rows: int, cols: int, values: ..float64) -> Matrix {
//     m := mat_zeros(rows, cols);
//     for values {
//         m.data[it_index] = c64(it, 0.0);
//     }
//     return m;
// }

mat :: (rows: int, cols: int, values: ..$T) -> Matrix(T) #must {
    #if CHECKS {
        assert(values.count <= rows * cols, "Too many values: rows*cols = %, values.count = %", rows*cols, values.count);
    }

    m := mat_zeros(T, rows, cols);
    for values {
        m.data[it_index] = it;
    }
    return m;
}


mat :: inline (dim: int, init: $T) -> Matrix(T) #must {
    return mat(dim, dim, init);
}




mat_zeros :: inline ($type: Type, rows: int, cols: int) -> Matrix(type) #must {
    return mat(rows, cols, zero(type));
}

mat_zeros :: inline ($type: Type, dim: int) -> Matrix(type) #must {
    return mat(dim, dim, zero(type));
}





mat_ones :: inline ($type: Type, rows: int, cols: int) -> Matrix(type) #must {
    return mat(rows, cols, one(type));
}

mat_ones :: inline ($type: Type, dim: int) -> Matrix(type) #must {
    return mat(dim, dim, one(type));
}



mat_id :: ($type: Type, rows: int, cols: int) -> Matrix(type) #must {
    m := mat_zeros(type, rows, cols);
    dim := ifx rows <= cols then rows else cols;
    for 0..dim-1 {
        m[it][it] = one(type);
    }
    return m;
}

mat_id :: inline ($type: Type, dim: int) -> Matrix(type) #must {
    return mat_id(type, dim, dim);
}



mat_hadamard_2 :: inline ($type: Type) -> Matrix(type) {
    return mat(2, 2, one(type), one(type), one(type), -one(type));
}

mat_hadamard :: inline ($type: Type, nestings: int) -> Matrix(type) {
    #if CHECKS {
        assert(nestings >= 0, "Nestings must be >=0 but = %", nestings);
    }
    if nestings == 0 {
        return mat_hadamard_2(type);
    }
    return tensor(mat_hadamard_2(type), mat_hadamard(type, nestings-1));
}




mat_reflector :: (normal: Vector($T)) -> Matrix(T) #must {
    tau := ncast(T, 2.0)/(conjugate(normal)*normal);

    res := mat_id(T, normal.dim);
    sub(res, tau * normal / normal);
    return res;
}


// inclusive first bound, exclusive last bound
mat_submatrix :: (M: Matrix($T), row_from: int, row_to: int, column_from: int, column_to: int) -> Matrix(T) #must {
    #if CHECKS {
        assert(row_to > row_from, "\ncannot create submatrix with rows %->%\n", row_to, row_from);
        assert(column_to > column_from, "\ncannot create submatrix with columns %->%\n", column_to, column_from);
    }
    res := mat(T, row_to-row_from, column_to-column_from);

    for i: 0..row_to-row_from-1 {
        for j: 0..column_to-column_from-1 {
            res[i][j] = M[i+row_from][j+column_from];
        }
    }

    return res;
}




// there are faster ways of doing this.
permute_rows :: (M: Matrix($T), order: ..int) -> Matrix(T) #must {
    #if CHECKS {
        assert(order.count == M.rows);
    }

    res := mat(T, M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[order[i]][j];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute_columns :: (M: Matrix($T), order: ..int) -> Matrix(T) #must {
    #if CHECKS {
        assert(order.count == M.cols);
    }

    res := mat(T, M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[i][order[j]];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute :: (M: Matrix($T), row_order: []int, col_order: []int) -> Matrix(T) #must {
    #if CHECKS {
        assert(row_order.count == M.rows);
        assert(col_order.count == M.cols);
    }

    res := mat(T, M.rows, M.cols);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            res[i][j] = M[row_order[i]][col_order[j]];
        }
    }

    return res;
}


swap_columns :: (M: *Matrix($T), i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0 && i != j);
        assert(i < M.cols && j < M.cols);
    }

    for 0..M.rows-1 {
        M[it][i], M[it][j] = swap(M[it][i], M[it][j]);
    }
}

swap_rows :: (M: *Matrix($T), i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0 && i != j);
        assert(i < M.rows && j < M.rows);
    }

    for 0..M.cols-1 {
        M[i][it], M[j][it] = swap(M[i][it], M[j][it]);
    }
}










transpose :: (m: Matrix($T)) -> Matrix(T) #must {
    res := mat(T, m.cols, m.rows);
    for r: 0..m.rows-1 {
        for c: 0..m.cols-1 {
            res[c][r] = m[r][c];
        }
    }
    return res;
}

// There are better ways to do this!
transpose :: (m: *Matrix($T)) {
    <<m = transpose(<<m);
}




conjugate :: (m: Matrix($T)) -> Matrix(T) #must {
    #if T != Complex64 {
        return m;
    } else {
        res := mat(T, m.cols, m.rows);
        for m.data {
            res.data[it_index] = conjugate(it);
        }
        return res;
    }
}

conjugate :: (m: *Matrix($T)) {
    #if T != Complex64 {
        return;
    } else {
        for m.data {
            m.data[it_index] = conjugate(it);
        }
    }
}


conjugate_transpose :: inline (m: Matrix($T)) -> Matrix(T) #must {
    return conjugate(transpose(m));
}
conjugate_transpose :: inline (m: *Matrix($T)) {
    transpose(m);
    conjugate(m);
}
dagger :: conjugate_transpose;
// hermitian :: conjugate_transpose;








// https://www.wikiwand.com/en/Tensor_product
#insert #run specialize("_tensor", "tensor", COMBI_M_M_M);
_tensor :: (a: Matrix($A), b: Matrix($B), res: *Matrix($C)) {
    <<res = mat(C, a.rows * b.rows, a.cols * b.cols);
    for r: 0..res.rows-1 {
        ar := r/b.rows;
        br := r % b.rows;
        
        for c: 0..res.cols-1 {
            ac := c/b.cols;
            bc := c % b.cols;

            (<<res)[r][c] = ncast(C, a[ar][ac]) * ncast(C, b[br][bc]);
            // set(res, r, c, get(a, ar, ac) * get(b, br, bc));
        }
    }
}









// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Lemma 3.5.4 Matrix Norms, page 176

// this needs the largest eigenvalue
// norm_2 :: inline (v: Vector) -> float64 {
//     return norm(v, 2.0);
// }

norm_1 :: (m: Matrix($T)) -> float64 {
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

norm_inf :: (m: Matrix($T)) -> float64 {
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

norm_frobenius :: (m: Matrix($T)) -> float64 {
    sum : float64 = 0.0;
    for m.data {
        sum += abs(it);
    }
    return sum;
}







#scope_file
// compiler bug: 
// the compiler won't find the operators generated with this method
// unless a stupid custom-named function is called first
// DUMMY_IN : [1][0]string : .[.[]];
// DUMMY_OUT : [1][1]string : .[.["int"]];
#insert #run specialize("_init", "init", .["=> int"]);
_init :: (ret: *int) {}



#scope_module

test_matrix :: () {
    print("\n\nMatrix(T) -------------------------------\n\n");
    init();
    {
        m := mat_zeros(float64, 3, 3);
        print("%\n", str(m));
    }
    {
        a := 3;
        m := mat_zeros(float64, a, a);
        print("%\n", str(m));
    }



    {
        m := mat(3, 3, 42.0);
        print("%\n", str(m));
    }
    {
        a := 3;
        m := mat(a, a, 42.0);
        print("%\n", str(m));
    }


    {
        m := mat(3, 3, c64(42.0));
        m[1][2] = R;
        assert(m[1][2] == c64(1.0, 0.0));
        print("m = \n%\n", str(m));
    }

    {
        m := mat_id(Complex64, 4);
        print("m = \n%\n", str(m));
        assert(m[0][1] == CZERO);
        assert(m[1][1] == R);
        assert(m[2][1] == CZERO);
        assert(m[3][3] == R);
    }

    {
        m := mat_zeros(Complex64, 3);
        m[0][2] = R;
        m = transpose(m);
        assert(m[2][0] == R);
    }

    {
        a := mat_id(float64, 2,2);
        b := mat_id(float64, 2,2);
        add(*a, b);
        assert(a == 2.0 * b);

        a = a + b;
        assert(a == 3.0 * b);
    }

    {
        m := mat(2, 2, R, c64(2.0), c64(3.0), c64(4.0));
        print("m = \n%\n", str(m));
        m = m * m;
        print("m = \n%\n", str(m));
    }
    {
        m := mat(2, 2, R, 2.0*R, 3.0*R, 4.0*R);
        print("m = \n%\n", str(m));
        m = m * m;
        print("m = \n%\n", str(m));
    }
    // {
    //     m := mat(1, 3, 1.0, 2.0, 3.0);
    //     print("m = %\n", m);
    //     mt := transpose(m);
    //     m = m * mt;
    //     print("m = %\n", m);
    // }
    {
        m := mat(2, 2, 
            CZERO, -I, 
            I, CZERO
        );
        assert( dagger(m) == m );
    }


    {
        m := mat(2, 2, R, I, c64(2.0, -1.0), 2.0 + 3.0*I);
        print("m = \n%\n", str(m));
        t1 := tensor(m, m);
        print("t1 = \n%\n", str(t1));
        t2 := mat(4, 4,
            R            ,             I,             I,   -R           ,
            2.0 -       I, 2.0 + 3.0 * I, 1.0 + 2.0 * I, -3.0 +  2.0 * I,
            2.0 -       I, 1.0 + 2.0 * I, 2.0 + 3.0 * I, -3.0 +  2.0 * I,
            3.0 - 4.0 * I, 7.0 + 4.0 * I, 7.0 + 4.0 * I, -5.0 + 12.0 * I
        );
        print("t2 = \n%\n", str(t2));
        assert(t1 == t2);
    }

    {
        M := mat(3,3, 
            CZERO, R, 2*R,
            3*R, 4*R, 5*R,
            6*R, 7*R, 8*R
        );

        S := mat_submatrix(M, 1, 3, 1, 3);
        assert(S == mat(2,2, 4*R, 5*R, 7*R, 8*R));
    }
}