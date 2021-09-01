






/*

 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##

*/

str :: (m: Matrix($T, $N, $M), indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    append(*builder, "\n");
    for i: 0..N-1 {
        append(*builder, indent);
        for j: 0..M-1 {
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

######## ##          ###     ######    ######
##       ##         ## ##   ##    ##  ##    ##
##       ##        ##   ##  ##        ##
######   ##       ##     ## ##   ####  ######
##       ##       ######### ##    ##        ##
##       ##       ##     ## ##    ##  ##    ##
##       ######## ##     ##  ######    ######

*/

// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Fig. 3.5 Matrix class derivation tree, page 268
MATRIX_FLAGS :: enum_flags u32 {
    NONE                   :: 0;
    SQUARE                 :: (1 <<  0);
    SYMMETRIC              :: (1 <<  1) + xx SQUARE;
    SYMMETRIC_POSITIVE     :: (1 <<  2) + xx SYMMETRIC;
    UPPER_HESSENBERG       :: (1 <<  3);
    LOWER_HESSENBERG       :: (1 <<  4); 
    ORTHOGONAL             :: (1 <<  5);
    TRAPEZOIDAL            :: (1 <<  6);
    LOWER_TRAPEZOIDAL      :: (1 <<  7) + xx TRAPEZOIDAL;
    LOWER_TRIANGULAR       :: (1 <<  8) + xx LOWER_TRAPEZOIDAL;
    UNIT_LOWER_TRAPEZOIDAL :: (1 <<  9) + xx LOWER_TRAPEZOIDAL;
    UNIT_LOWER_TRIANGULAR  :: (1 << 10) + xx UNIT_LOWER_TRAPEZOIDAL;
    UPPER_TRAPEZOIDAL      :: (1 << 11) + xx TRAPEZOIDAL;
    UPPER_TRIANGULAR       :: (1 << 12) + xx UPPER_TRAPEZOIDAL;
    UNIT_UPPER_TRAPEZOIDAL :: (1 << 13) + xx UPPER_TRAPEZOIDAL;
    UNIT_UPPER_TRIANGULAR  :: (1 << 14) + xx UNIT_UPPER_TRAPEZOIDAL;    
    HERMITIAN              :: (1 << 15) + xx SQUARE; // not the same as symmetric for complex matrices
    DIAGONAL               :: (1 << 16) + xx BAND;
    BAND                   :: (1 << 17) + xx TRAPEZOIDAL;
}


is :: inline (m: Matrix($T, $N, $M), flag: MATRIX_FLAGS) -> bool {
    return (m.flags & flag) > 0;
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



Matrix :: struct(T: Type, N: int, M: int) {
    data: [N*M]T;
    flags: MATRIX_FLAGS = .NONE;
}

operator [] :: inline (m: Matrix($T, $N, $M), i: int) -> []T {
    #if CHECKS {
        assert(i >= 0, "Row-index negative! i = %", i);
        assert(i < N, "Row-index out of bounds: % rows", N);
    }

    return array_view_count(m.data, i * M, M);
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


mat_cast :: ($type: Type, m: Matrix($T, $N, $M)) -> Matrix(type, N, M) #must {
    res : Matrix(type, N, M) = ---;
    res.flags = m.flags;

    for m.data {
        res.data[it_index] = it;
    }

    return res;
}

copy :: (m: Matrix($T, N, M)) -> Matrix(T) #must {
    res := mat(T, m.rows, m.cols, m.flags);
    res.data = array_copy(m.data);
    return res;
}

copy :: (to: *Matrix($T, $N, $M), from: Matrix(T, N, M)) {
    to.data = array_copy(from.data);
}


mat :: ($rows: int, $cols: int, init: $T, $flags: MATRIX_FLAGS = .NONE) -> Matrix(T, rows, cols) #must {
    m : Matrix(T, rows, cols) = ---;
    m.flags = flags;
    for 0..N*M-1 {
        m.data[it] = init;
    }
    return m;
}



mat :: ($rows: int, $cols: int, values: [rows * cols]$T, flags: MATRIX_FLAGS = .NONE) -> Matrix(T, rows, cols) #must {
    m : Matrix(T, rows, cols);
    m.flags = flags;
    for values {
        m.data[it_index] = it;
    }
    return m;
}


mat :: inline ($dim: int, init: $T, flags: MATRIX_FLAGS = .NONE) -> Matrix(T, dim, dim) #must {
    return mat(dim, dim, init, flags);
}


mat_ones :: inline ($type: Type, $rows: int, $cols: int, flags: MATRIX_FLAGS = .NONE) -> Matrix(type, rows, cols) #must {
    return mat(rows, cols, one(type), flags);
}

mat_ones :: inline ($type: Type, $dim: int, flags: MATRIX_FLAGS = .NONE) -> Matrix(type, dim, dim) #must {
    return mat(dim, dim, one(type), flags);
}



mat_id :: ($type: Type, $rows: int, $cols: int) -> Matrix(type, rows, cols) #must {
    m : Matrix(type, rows, cols);
    m.flags = .DIAGONAL;

    #if rows <= cols {
        dim :: rows;
    } else {
        dim :: cols;
    }

    for 0..dim-1 {
        m[it][it] = one(type);
    }
    return m;
}

mat_id :: inline ($type: Type, $dim: int) -> Matrix(type, dim, dim) #must {
    return mat_id(type, dim, dim);
}



mat_hadamard_2 :: inline ($type: Type) -> Matrix(type, 2, 2) {
    // TODO: update flags
    return mat(2, 2, one(type), one(type), one(type), -one(type), .NONE);
}

mat_hadamard :: ($type: Type, $nestings: int) -> Matrix(type, D, D) #modify {
    D = 1 << nestings;
} {
    #if CHECKS {
        #assert(nestings >= 0, "Nestings must be >=0 but = %", nestings);
    }
    #if nestings == 0 {
        return mat_hadamard_2(type);
    } else {
        return tensor(mat_hadamard_2(type), mat_hadamard(type, nestings-1));
    }
}


mat_reflector :: (normal: Vector($T, $N)) -> Matrix(T, N, N) #must {
    tau := ncast(T, 2)/(conjugate(normal)*normal);

    res := mat_id(T, N);
    sub(res, tau * normal / normal);
    return res;
}






// inclusive first bound, exclusive last bound
mat_submatrix :: (M: Matrix($T, $N, $M), $row_from: int, $row_to: int, $column_from: int, $column_to: int) -> Matrix(T, row_to-row_from, column_to-column_from) #must {
    #if CHECKS {
        #assert(row_to > row_from, "\ncannot create submatrix with rows %->%\n", row_to, row_from);
        #assert(column_to > column_from, "\ncannot create submatrix with columns %->%\n", column_to, column_from);
    }
    res : Matrix(T, row_to-row_from, column_to-column_from) = ---;
    res.flags = .NONE;

    for i: 0..row_to-row_from-1 {
        for j: 0..column_to-column_from-1 {
            res[i][j] = M[i+row_from][j+column_from];
        }
    }

    return res;
}









/*
########   #######  ##      ##          ######   #######  ##       ##     ## ##     ## ##    ##
##     ## ##     ## ##  ##  ##         ##    ## ##     ## ##       ##     ## ###   ### ###   ##
##     ## ##     ## ##  ##  ##         ##       ##     ## ##       ##     ## #### #### ####  ##
########  ##     ## ##  ##  ##         ##       ##     ## ##       ##     ## ## ### ## ## ## ##
##   ##   ##     ## ##  ##  ##         ##       ##     ## ##       ##     ## ##     ## ##  ####
##    ##  ##     ## ##  ##  ##         ##    ## ##     ## ##       ##     ## ##     ## ##   ###
##     ##  #######   ###  ###           ######   #######  ########  #######  ##     ## ##    ##
*/


row :: (m: Matrix($T, $N, $M), i: int) -> Vector(T, M) #must {
    #if CHECKS {
        assert(i >= 0 && i < N);
    }
    res : Vector(T, M);
    for 0..M-1 {
        res.data[it] = m[i][it];
    }
    return res;
}

column :: (m: Matrix($T, $N, $M), j: int) -> Vector(T, N) #must {
    #if CHECKS {
        assert(i >= 0 && i < M);
    }
    res : Vector(T, N) = ---;
    for 0..N-1 {
        res.data[it] = m[it][i];
    }
    return res;
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




operator == :: (a: Matrix($A, $aN, $aM), b: Matrix($B, $bN, $bM)) -> bool {
    #if aN != bN || aM != bM {
        return false;
    } else {
        for a.data {
            if it != b.data[it_index] {
                return false;
            }
        }
        return true;
    }
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
operator * :: (a: Matrix($A, $aN, $aM), b: Matrix($B, $bN, $bM)) -> Matrix($C, aN, bM) #modify {
    C = dominant_type(A, B);
}{
    #if CHECKS {
        #assert(aM == bN, "Matrix dimensions don't match: a.cols = %, b.rows = %", aM, bN);
    }

    res : Matrix(C, aN, bM);
    for j: 0..bM-1 {
        for k: 0..aM-1 {
            tmp := b[k][j];
            for i: 0..aN-1 {
                res[i][j] += a[i][k] * tmp;
            }
        }
    }

    return res;
}

operator * :: (a: Matrix($A, $aN, $aM), b: Vector($B, $bN)) -> Vector($C, aN) #modify {
    C = dominant_type(A, B);
}{
    #if CHECKS {
        #assert(aM == bN, "Dimensions mismatch: a.cols = %, b.dim = %", aM, bN);
    }

    res : Vector(C, aN);
    for i: 0..aN-1 {
        for j: 0..aM-1 {
            res.data[i] += a[i][j] * b[j];
        }
    }

    return res;
}

operator * :: (a: Vector($A, $aN), b: Matrix($B, $bN, $bM)) -> Vector($C, bM) #modify {
    C = dominant_type(A, B);
}{
    #if CHECKS {
        #assert(aN == bN, "Dimensions mismatch: a.dim = %, b.rows = %", aN, bN);
    }

    res : Vector(C, bM);
    for j: 0..bN-1 {
        for i: 0..aN-1 {
            res.data[j] += a[i] * b[i][j];
        }
    }

    return res;
}

operator * :: (a: Matrix($A, $N, $M), b: $B) -> Matrix($C, N, M) #modify {
    if !is_scalar(B) { C = null; return; }
    C = dominant_type(A, B);
} #symmetric {

    res : Matrix(C, N, M) = ---;
    res.flags = a.flags;

    for 0..a.data.count-1 {
        res.data[it] = a.data[it] * b;
    }

    return res;
}


mul :: (a: *Matrix($A, $N, $M), b: $B) {
    #if CHECKS {
        #assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] *= b;
    }
}






operator / :: inline (a: Matrix($A, $aN, $aM), b: Matrix($B, $bN, $bM)) -> Matrix($C, aM, bM) #modify {
    C = dominant_type(A, B);
} { 
    return inverse(b) * a;
}

operator / :: (a: Vector($A, $aN), b: Matrix($B, $bN, $bM)) -> Matrix($C, bM, aN) #modify {
    C = dominant_type(A, B);
} { 
    return inverse(b) * a;
}

operator / :: inline (a: Matrix($A, $N, $M), b: $B) -> Matrix($C, N, M) #modify {
    if !is_scalar(B) { C = null; return; }
    C = dominant_type(A, B);
} { 
    res : Matrix(C, N, M) = ---;
    res.flags = a.flags;

    for 0..a.data.count-1 {
        res.data[it] = a.data[it] / b;
    }
    return res;
}



div :: (a: *Matrix($A, $N, $M), b: $B) {
    #if CHECKS {
        #assert(is_scalar(B));
        #assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] /= b;
    }
}




operator + :: (a: Matrix($A, $N, $M), b: Matrix($B, N, M)) -> Matrix($C, N, M) #modify {
    C = dominant_type(A, B);
}{
    res : Matrix(C, N, M) = ---;
    res.flags = .NONE;

    for 0..a.data.count-1 {
        res.data[it] = a.data[it] + b.data[it];
    }
    return res;
}





add :: (a: *Matrix($A, $N, $M), b:Matrix($B, N, M)) {
    #if CHECKS {
        assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] += b.data[it];
    }
}








operator - :: (a: Matrix($A, $N, $M), b: Matrix($B, N, M)) -> Matrix($C, N, M) #modify {
    C = dominant_type(A, B);
}{
    res : Matrix(C, N, M) = ---;
    res.flags = .NONE;

    for 0..a.data.count-1 {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}



sub :: (a: *Matrix($A, $N, $M), b:Matrix($B, N, M)) {
    #if CHECKS {
        #assert(is_minor_type(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] -= b.data[it];
    }
}





operator - :: (a: Matrix($A, $N, $M)) -> Matrix(A, N, M) {
    res : Matrix(A, N, M) = ---;
    res.flags = a.flags;

    for 0..a.data.count-1 {
        res.data[it] = -a.data[it];
    }
    return res;
}




neg :: (a: Matrix($T, $N, $M)) {
    for a.data {
        a.data[it_index] = -it;
    }
}















// Matrix Computations, 4th edition
// The Johns Hopkins University Press
// Algorithm 1.1.6 Dot Product Matrix Multiplication, page 10
// Algorithm 1.1.7 Saxpy Matrix Multiplication, page 10
// C = C + A*B
saxpy :: (c: *Matrix($C, $cN, $cM), a: Matrix($A, cN, $aM), b: Matrix($B, aM, cM)) {
    #if CHECKS {
        #assert(is_minor_type(C, A) && is_minor_type(C, B));
    }

    for j: 0..cM-1 {
        for k: 0..aM-1 {
            tmp := b[k][j];
            for i: 0..cN-1 {
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
permute_rows :: (m: Matrix($T, $N, $M), order: ..int) -> Matrix(T, N, M) #must {
    #if CHECKS {
        assert(order.count == N);
    }

    res : Matrix(T, N, M) = ---;
    res.flags = .NONE;

    for i: 0..N-1 {
        for j: 0..M-1 {
            res[i][j] = m[order[i]][j];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute_columns :: (m: Matrix($T, $N, $M), order: ..int) -> Matrix(T, N, M) #must {
    #if CHECKS {
        assert(order.count == M);
    }

    res : Matrix(T, N, M) = ---;
    res.flags = .NONE;

    for i: 0..N-1 {
        for j: 0..M-1 {
            res[i][j] = m[i][order[j]];
        }
    }

    return res;
}

// there are faster ways of doing this.
permute :: (m: Matrix($T, $N, $M), row_order: []int, col_order: []int) -> Matrix(T, N, M) #must {
    #if CHECKS {
        assert(row_order.count == N);
        assert(col_order.count == M);
    }

    res : Matrix(T, N, M) = ---;
    res.flags = .NONE;

    for i: 0..N-1 {
        for j: 0..M-1 {
            res[i][j] = m[row_order[i]][col_order[j]];
        }
    }

    return res;
}


swap_columns :: (m: *Matrix($T, $N, $M), i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0);
        assert(i < M && j < M);
    }

    if i == j then return;

    for 0..N-1 {
        (<<m)[it][i], (<<m)[it][j] = swap((<<m)[it][i], (<<m)[it][j]);
    }
}

swap_rows :: (m: *Matrix($T, $N, $M), i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && j >= 0);
        assert(i < m.rows && j < m.rows);
    }

    if i == j then return;

    for 0..M-1 {
        (<<m)[i][it], (<<m)[j][it] = swap((<<m)[i][it], (<<m)[j][it]);
    }
}





transpose :: (m: Matrix($T, $N, $M)) -> Matrix(T, M, N) #must {
    res : Matrix(T, M, N);
    for r: 0..N-1 {
        for c: 0..M-1 {
            res[c][r] = m[r][c];
        }
    }
    return res;
}

// There are better ways to do this!
transpose :: (m: *Matrix($T, $N, $M)) {
    <<m = transpose(<<m);
}





conjugate :: (m: Matrix($T, $N, $M)) -> Matrix(T, N, M) #must {
    #if T != Complex64 {
        return m;
    } else {
        res : Matrix(T, N, M);
        res.flags = m.flags;

        for m.data {
            res.data[it_index] = conjugate(it);
        }
        return res;
    }
}

conjugate :: (m: *Matrix($T, $N, $M)) {
    #if T != Complex64 {
        return;
    } else {
        for m.data {
            m.data[it_index] = conjugate(it);
        }
    }
}


conjugate_transpose :: inline (m: Matrix($T, $N, $M)) -> Matrix(T, M, N) #must {
    return conjugate(transpose(m));
}
conjugate_transpose :: inline (m: *Matrix($T, $N, $M)) {
    transpose(m);
    conjugate(m);
}
dagger :: conjugate_transpose;
// hermitian :: conjugate_transpose;





// https://www.wikiwand.com/en/Tensor_product
tensor :: (a: Matrix($A, $aN, $aM), b: Matrix($B, $bN, $bM)) -> Matrix($C, aN*bN, aM*bM) #modify {
    C = dominant_type(A, B);
} {
    res : Matrix(C, aN * bN, aM * bM) = ---;
    res.flags = .NONE;

    for r: 0..res.M-1 {
        ar := r/bN;
        br := r % bN;
        
        for c: 0..res.M-1 {
            ac := c/bM;
            bc := c % bM;

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

norm_1 :: (m: Matrix($T, $N, $M)) -> float64 {
    res : float64 = 0.0;

    for j: 0..M {
        sum : float64 = 0.0;
        for i: 0..N {
            sum += abs(m[i][j]);
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}

norm_inf :: (m: Matrix($T, $N, $M)) -> float64 {
    res : float64 = 0.0;

    for i: 0..N {
        sum : float64 = 0.0;
        for j: 0..M {
            sum += abs(m[i][j]);
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}

norm_frobenius :: (m: Matrix($T, $N, $M)) -> float64 {
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
        m : Matrix(float64, 3, 3);
        print("%\n", str(m));
    }
    {
        a :: 3;
        m : Matrix(float64, a, a);
        print("%\n", str(m));
    }



    {
        m := mat(3, 3, 42.0);
        print("%\n", str(m));
    }
    {
        a :: 3;
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
        m : Matrix(Complex64, 3, 3);
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
        m : Matrix(float64, 4, 4) = ---;
        m.flags = .UNIT_LOWER_TRIANGULAR;
        assert( is(m, .TRAPEZOIDAL) );
    }
}