






/*

 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##

*/

str :: (m: TriangularMatrix($T), indent: string = "    ") -> string {
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







/*
 ######  ######## ########  ##     ##  ######  ########
##    ##    ##    ##     ## ##     ## ##    ##    ##
##          ##    ##     ## ##     ## ##          ##
 ######     ##    ########  ##     ## ##          ##
      ##    ##    ##   ##   ##     ## ##          ##
##    ##    ##    ##    ##  ##     ## ##    ##    ##
 ######     ##    ##     ##  #######   ######     ##
*/



TriangularMatrix :: struct(T: Type) {
    using M : Matrix;
    isUpper : bool;
}

get :: inline (m: TriangularMatrix($T), i: int, j: int) -> T {
    #if CHECKS {
        assert(
            (m.isUpper && j >= i) || (!m.isUpper && i >= j),
            "index out of bounds!"
        );
    }
    r : int = ---;
    if m.isUpper {
        // 1 + 2 + ... + n   = n*(n+1)/2
        // 0 + 1 + ... + n-1 = (n-1)*n/2
        r = i * m.cols - (i-1) * i/2
    } else {
        r = i * (i+1) / 2;
    }
    return m.data[r + j];  
}

get :: inline (m: TriangularMatrix($T), i: int) -> T {
    return m.data[i];
}

set :: inline (m: TriangularMatrix($T), i: int, j: int, val: $U) {
    #if CHECKS {
        assert(
            (m.isUpper && j >= i) || (!m.isUpper && i >= j),
            "index out of bounds!"
        );
    }
    r : int = ---; 
    if m.isUpper {
        // 1 + 2 + ... + n   = n*(n+1)/2
        // 0 + 1 + ... + n-1 = (n-1)*n/2
        r = i * m.cols - (i-1) * i/2
    } else {
        r = i * (i+1) / 2;
    }
    m.data[r + j] = ncast(T, val);
} 

set :: inline (m: TriangularMatrix($T), i: int, val: $U) {
    m.data[i] = ncast(T, val);
}









/*
#### ##    ## #### ######## ####    ###    ##       #### ######## ########
 ##  ###   ##  ##     ##     ##    ## ##   ##        ##       ##  ##
 ##  ####  ##  ##     ##     ##   ##   ##  ##        ##      ##   ##
 ##  ## ## ##  ##     ##     ##  ##     ## ##        ##     ##    ######
 ##  ##  ####  ##     ##     ##  ######### ##        ##    ##     ##
 ##  ##   ###  ##     ##     ##  ##     ## ##        ##   ##      ##
#### ##    ## ####    ##    #### ##     ## ######## #### ######## ########
*/


#scope_file
data_size :: (rows: int, cols: int) -> int {
    minDim := min(rows, cols);
    maxDim := max(rows, cols);
    trisize := minDim * (minDim+1) / 2;
    return trisize + minDim * (maxDim - minDim);
    // e.g.
    // M = | 1 2 : 3 4 5|
    //     |   6 : 7 8 9|
}
#scope_export


tri_mat :: ($type: Type, rows: int, cols: int, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(type) #must {
    res : TriangularMatrix(type) = ---;
    res.rows = rows;
    res.cols = cols;
    res.flags = flags;
    res.isUpper = isUpper;
    res.data = NewArray(data_size(rows, cols), type, initialized=false);
    return res;
}

tri_mat :: ($type: Type, m: $M/TriangularMatrix) -> TriangularMatrix(type) #must {
    res : TriangularMatrix(type) = ---;
    res.rows = m.rows;
    res.cols = m.cols;
    res.isUpper = m.isUpper;
    res.flags = m.flags;
    res.data = NewArray(m.data.count, type, initialized=false);
    for 0..m.data.count-1 {
        set(res, it, ncast(type, get(m, it)));
    }
    return res;
}

copy :: (m: TriangularMatrix($T)) -> TriangularMatrix(T) #must {
    res := tri_mat(T, m.rows, m.cols, m.flags, m.isUpper);
    res.data = array_copy(m.data);
    return res;
}

copy :: (to: *TriangularMatrix($T), from: TriangularMatrix(T)) {
    to.rows = from.rows;
    to.cols = from.cols;
    to.isUpper = from.isUpper;
    to.data = array_copy(from.data);
}


tri_mat :: inline ($type: Type, dim: int, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(type) #must {
    return tri_mat(type, dim, dim, isUpper, flags);
}


tri_mat :: (rows: int, cols: int, init: $T, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(T) #must {
    m := tri_mat(T, rows, cols, isUpper, flags);
    for 0..m.data.count-1 {
        m.data[it] = init;
    }
    return m;
}



tri_mat :: (rows: int, cols: int, values: ..$T, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(T) #must {
    #if CHECKS {
        assert(values.count <= data_size(rows, cols), "Too many values: rows*cols = %, values.count = %", rows*cols, values.count);
    }

    m := tri_mat_zeros(T, rows, cols, isUpper, flags);
    for values {
        m.data[it_index] = it;
    }
    return m;
}


tri_mat :: inline (dim: int, init: $T, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(T) #must {
    return tri_mat(dim, dim, init, isUpper, flags);
}




tri_mat_zeros :: inline ($type: Type, rows: int, cols: int, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(type) #must {
    return tri_mat(rows, cols, zero(type), isUpper, flags);
}

tri_mat_zeros :: inline ($type: Type, dim: int, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(type) #must {
    return tri_mat(dim, dim, zero(type), isUpper, flags);
}




tri_mat_ones :: inline ($type: Type, rows: int, cols: int, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(type) #must {
    return tri_mat(rows, cols, one(type), isUpper, flags);
}

tri_mat_ones :: inline ($type: Type, dim: int, isUpper: bool, flags: MATRIX_FLAGS = .NONE) -> TriangularMatrix(type) #must {
    return tri_mat(dim, dim, one(type), isUpper, flags);
}



tri_mat_id :: ($type: Type, rows: int, cols: int, isUpper: bool) -> TriangularMatrix(type) #must {
    m := tri_mat_zeros(type, rows, cols, isUpper, .DIAGONAL);
    dim := ifx rows <= cols then rows else cols;
    for 0..dim-1 {
        set(m, it, it, one(type));
    }
    return m;
}

tri_mat_id :: inline ($type: Type, dim: int, isUpper: bool) -> TriangularMatrix(type) #must {
    return tri_mat_id(type, dim, dim, isUpper);
}







/*
 #######  ########  ######## ########     ###    ########  #######  ########   ######
##     ## ##     ## ##       ##     ##   ## ##      ##    ##     ## ##     ## ##    ##
##     ## ##     ## ##       ##     ##  ##   ##     ##    ##     ## ##     ## ##
##     ## ########  ######   ########  ##     ##    ##    ##     ## ########   ######
##     ## ##        ##       ##   ##   #########    ##    ##     ## ##   ##         ##
##     ## ##        ##       ##    ##  ##     ##    ##    ##     ## ##    ##  ##    ##
 #######  ##        ######## ##     ## ##     ##    ##     #######  ##     ##  ######
*/




operator == :: (a: Matrix($A), b: Matrix($B)) -> bool {
    if a.rows != b.rows || a.cols != b.cols {
            return false;
    }
    for a.data {
        if it != b.data[it_index] {
            return false;
        }
    }
    return true;
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 2.11 Case Study: Matrix-Matrix Multiplication, page 118
// According to Table 2.12, the best performing loop order is j/k/i with
// the least number of cache misses.
//
// Similar discussion in:
// Matrix Computations, 4th edition
// The Johns Hopkins University Press
// Algorithm 1.1.5 ijk Matrix Multiplication, page 8
operator * :: (a: Matrix($A), b: Matrix($B)) -> Matrix($C) #modify {
    C = upcast(A, B);
}{
    #if CHECKS {
        assert(a.cols == b.rows, "Matrix dimensions don't match: a.cols = %, b.rows = %", a.cols, b.rows);
    }

    res := mat_zeros(C, a.rows, b.cols);
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

operator * :: (a: Matrix($A), b: Vector($B)) -> Vector($C) #modify {
    C = upcast(A, B);
}{
    #if CHECKS {
        assert(a.cols == b.dim, "Dimensions mismatch: a.cols = %, b.dim = %", a.cols, b.dim);
    }

    res := vec_zeros(C, a.rows);
    for i: 0..a.rows-1 {
        for j: 0..a.cols-1 {
            res.data[i] += a[i][j] * b[j];
        }
    }

    return res;
}

operator * :: (a: Vector($A), b: Matrix($B)) -> Vector($C) #modify {
    C = upcast(A, B);
}{
    #if CHECKS {
        assert(a.dim == b.rows, "Dimensions mismatch: a.dim = %, b.rows = %", a.dim, b.rows);
    }

    res := vec_zeros(C, b.cols);
    for j: 0..b.rows-1 {
        for i: 0..a.dim-1 {
            res.data[j] += a[i] * b[i][j];
        }
    }

    return res;
}

operator * :: (a: Matrix($A), b: $B) -> Matrix($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A, B);
} #symmetric {

    res := mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] * b;
    }

    return res;
}


mul :: (a: *Matrix($A), b: $B) {
    #if CHECKS {
        assert(is_minor_number(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] *= b;
    }
}






operator / :: inline (a: Matrix($A), b: Matrix($B)) -> Matrix($C) #modify {
    C = upcast(A, B);
} { 
    return inverse(b) * a;
}

operator / :: (a: Vector($A), b: Matrix($B)) -> Matrix($C) #modify {
    C = upcast(A, B);
} { 
    return inverse(b) * a;
}

operator / :: inline (a: Matrix($A), b: $B) -> Matrix($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A, B);
} { 
    res := mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] / b;
    }
    return res;
}



div :: (a: *Matrix($A), b: Complex64) {
    #if CHECKS {
        assert(is_minor_number(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] /= b;
    }
}
div :: (a: *Matrix($A), b: float64) {
    #if CHECKS {
        assert(is_minor_number(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] /= b;
    }
}







operator + :: (a: Matrix($A), b: Matrix($B)) -> Matrix($C) #modify {
    C = upcast(A, B);
}{
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] + b.data[it];
    }
    return res;
}





add :: (a: *Matrix($A), b:Matrix($B)) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
        assert(is_minor_number(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] += b.data[it];
    }
}








operator - :: (a: Matrix($A), b: Matrix($B)) -> Matrix($C) #modify {
    C = upcast(A, B);
}{
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
    }

    res := mat(C, a.rows, a.cols);
    for 0..a.data.count-1 {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}



sub :: (a: *Matrix($A), b:Matrix($B)) {
    #if CHECKS {
        assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
        assert(is_minor_number(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] -= b.data[it];
    }
}





operator - :: (a: Matrix($A)) -> Matrix(A) {
    res := mat(A, a.rows, a.cols);
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















// Matrix Computations, 4th edition
// The Johns Hopkins University Press
// Algorithm 1.1.6 Dot Product Matrix Multiplication, page 10
// Algorithm 1.1.7 Saxpy Matrix Multiplication, page 10
// C = C + A*B
saxpy :: (c: *Matrix($C), a: Matrix($A), b: Matrix($B)) {
    #if CHECKS {
        assert(is_minor_number(C, A) && is_minor_number(C, B));
        assert(c.rows == a.rows);
        assert(c.cols == b.cols);
        assert(a.cols == b.rows);
    }

    for j: 0..c.cols-1 {
        for k: 0.. a.cols-1 {
            tmp := b[k][j];
            for i: 0..c.rows-1 {
                c[i][j] += a[i][k] * tmp;
            }
        }
    }
}















/*
######## ########     ###    ##    ##  ######  ########
   ##    ##     ##   ## ##   ###   ## ##    ## ##
   ##    ##     ##  ##   ##  ####  ## ##       ##
   ##    ########  ##     ## ## ## ##  ######  ######
   ##    ##   ##   ######### ##  ####       ## ##
   ##    ##    ##  ##     ## ##   ### ##    ## ##
   ##    ##     ## ##     ## ##    ##  ######  ##
*/



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
        (<<M)[it][i], (<<M)[it][j] = swap((<<M)[it][i], (<<M)[it][j]);
    }
}

swap_rows :: (M: *Matrix($T), i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0 && i != j);
        assert(i < M.rows && j < M.rows);
    }

    for 0..M.cols-1 {
        (<<M)[i][it], (<<M)[j][it] = swap((<<M)[i][it], (<<M)[j][it]);
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
tensor :: (a: Matrix($A), b: Matrix($B)) -> Matrix($C) #modify {
    C = upcast(A, B);
} {
    res := mat(C, a.rows * b.rows, a.cols * b.cols);
    for r: 0..res.rows-1 {
        ar := r/b.rows;
        br := r % b.rows;
        
        for c: 0..res.cols-1 {
            ac := c/b.cols;
            bc := c % b.cols;

            res[r][c] = a[ar][ac] * b[br][bc];
        }
    }
    return res;
}





/*
        ##        ######   ######     ###    ##          ###    ########
         ##      ##    ## ##    ##   ## ##   ##         ## ##   ##     ##
          ##     ##       ##        ##   ##  ##        ##   ##  ##     ##
#######    ##     ######  ##       ##     ## ##       ##     ## ########
          ##           ## ##       ######### ##       ######### ##   ##
         ##      ##    ## ##    ## ##     ## ##       ##     ## ##    ##
        ##        ######   ######  ##     ## ######## ##     ## ##     ##
*/



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





////////////////////////////////////////////////////////////////////////////////////////
/*
      ######## ########  ######  ########  ######
           ##    ##       ##    ##    ##    ##    ##
           ##    ##       ##          ##    ##
           ##    ######    ######     ##     ######
           ##    ##             ##    ##          ##
           ##    ##       ##    ##    ##    ##    ##
           ##    ########  ######     ##     ######
*/
////////////////////////////////////////////////////////////////////////////////////////

#scope_module

test_matrix :: () {
    print("\n\nMatrix(T) -------------------------------\n\n");
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
        print("a = \n%\n", str(a));
        print("b = \n%\n", str(b));
        add(*a, b);
        print("a = \n%\n", str(a));
        print("2 b = \n%\n", 2.0 * b);
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

    {
        print("%\n", MATRIX_FLAGS.UNIT_LOWER_TRIANGULAR);
        m : Matrix(float64) = ---;
        m.flags = .UNIT_LOWER_TRIANGULAR;
        assert( is(m, .TRAPEZOIDAL) );
    }
}