






/*

 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##

*/

pstr :: inline (m: $V/MatrixType) -> string {
    return str(m, print_options.offset, print_options.shift, "│", print_options.indent_charactoer_every);
}

str :: (m: $M/MatrixType, offset: int = 0, indent: int = 4, indent_character: string = " ", indent_charactoer_every: int = 4) -> string {
    
    shift :: (ind: int = 0, additional_indent: int = 0) #expand {
        for 0..offset-1 {
            append(*builder, " ");
        }
        for 0..ind-1+additional_indent {
            if it % indent_charactoer_every == 0 {
                append(*builder, indent_character);
            } else {
                append(*builder, " ");
            }
        }
    }
    
    builder: String_Builder;
    defer free_buffers(*builder);
    init_string_builder(*builder);

    append(*builder, "\n");

    column_width : [M.M]int; // * initialize all to 0
    for m {
        column_width[it_column] = max(
            column_width[it_column], 
            str(it).count
        );
    }

    matrix_width := 0;
    for column_width {
        matrix_width += it + 2;
    }
    matrix_width -= 2;


    for i: 0..M.N-1 {
        shift(indent+4);

        if i == 0 {
            append(*builder, "⎛ ");
        } else if i == M.N-1 {
            append(*builder, "⎝ ");
        } else {
            append(*builder, "⎜ ");
        }

        for j: 0..M.M-1 {
            padding := column_width[j] - str(get(m,i,j)).count;
            for 0..padding-1 {
                append(*builder, " ");
            }
            append(*builder, str(get(m,i,j)) );
            if j != M.M-1 {
                append(*builder, ", ");
            }
        }

        if i == 0 {
            append(*builder, " ⎞");
        } else if i == M.N-1 {
            append(*builder, " ⎠");
        } else {
            append(*builder, " ⎟");
        }

        append(*builder, "\n");
    }

    shift(indent);
    append(*builder, "└╴");
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

// MATRIX_FLAGS :: enum_flags u32 {
//     NONE                   :: 0;
//     SQUARE                 :: (1 <<  0);
//     SYMMETRIC              :: (1 <<  1) + xx SQUARE;
//     SYMMETRIC_POSITIVE     :: (1 <<  2) + xx SYMMETRIC;
//     UPPER_HESSENBERG       :: (1 <<  3);
//     LOWER_HESSENBERG       :: (1 <<  4); 
//     ORTHOGONAL             :: (1 <<  5);
//     TRAPEZOIDAL            :: (1 <<  6);
//     LOWER_TRAPEZOIDAL      :: (1 <<  7) + xx TRAPEZOIDAL;
//     LOWER_TRIANGULAR       :: (1 <<  8) + xx LOWER_TRAPEZOIDAL;
//     UNIT_LOWER_TRAPEZOIDAL :: (1 <<  9) + xx LOWER_TRAPEZOIDAL;
//     UNIT_LOWER_TRIANGULAR  :: (1 << 10) + xx UNIT_LOWER_TRAPEZOIDAL;
//     UPPER_TRAPEZOIDAL      :: (1 << 11) + xx TRAPEZOIDAL;
//     UPPER_TRIANGULAR       :: (1 << 12) + xx UPPER_TRAPEZOIDAL;
//     UNIT_UPPER_TRAPEZOIDAL :: (1 << 13) + xx UPPER_TRAPEZOIDAL;
//     UNIT_UPPER_TRIANGULAR  :: (1 << 14) + xx UNIT_UPPER_TRAPEZOIDAL;    
//     HERMITIAN              :: (1 << 15) + xx SQUARE; // not the same as symmetric for complex matrices
//     DIAGONAL               :: (1 << 16) + xx BAND;
//     BAND                   :: (1 << 17) + xx TRAPEZOIDAL;
// }


is :: inline (m: $M/MatrixType, flag: MATRIX_FLAGS) -> bool {
    return (m.flags & flag) > 0;
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


mat_cast :: ($type: Type, m: $M/MatrixType) -> MatrixCastHelper(M, type).T #must {
    res := make(MatrixCastHelper(M, type).T, false);

    for raw(*m) {
        set_raw(*res, it_raw_index, ncast(type, it));
    }

    return res;
}

copy :: (m: $M/MatrixType) -> M #must {
    res := make(M, false);

    for raw(*m) {
        set_raw(*res, it_raw_index, it);
    }

    return res;
}

copy :: inline (to: *$M/MatrixType, from: M) {
    for from {
        set_lin(to, it_index, it);
    }
}

mat_init :: (m: *$M/MatrixType, vals: ..$T) {
    #if CHECKS {
        #run assert(is_minor_number_type(M.data_type, T));
        assert(vals.count <= M.N * M.M, "too many values (%) for matrix of % entries\n", vals.count, M.N * M.M);
    }
    // println("mat_init");
    // println("    vals: %", vals);

    for vals {
        j := it_index % M.M;
        i := (it_index - j) / M.M;
        // println("    original %, cast %", it, ncast(M.data_type, it));
        set(m, i, j, ncast(M.data_type, it));
    }
}

mat_init_all :: (m: *$M/MatrixType, value: $T) {
    for m {
        set_lin(m, it_index, ncast(M.data_type, value));
    }
}

mat_init_unit :: (m: *$M/MatrixType, $init_zeroes: bool = false) {
    #if init_zeroes {
        for m {
            if it_row == it_column {
                set_lin(m, it_index, zero(M.data_type));
            } else {
                set_lin(m, it_index, one(M.data_type));
            }
        }
    } else {
        for i: 0..min(M.N, M.M)-1 {
            set(m, i, i, one(M.data_type));
        }
    }
}

mat_init_diag :: (m: *$M/MatrixType, values: ..$T, $init_zeroes: bool = false) {
    #if CHECKS {
        assert(values.count == min(M.N, M.M), "not enough values (% values) for diagonal (% values)\n", values.count, min(M.N, M.M));
    }
    #if init_zeroes {
        for m {
            if it_column == it_row {
                set_lin(m, it_index, ncast(M.data_type, values[it_row]));
            } else {
                set_lin(m, it_index, zero(M.data_type));
            }
        }
    } else {
        for i: 0..min(M.N, M.M)-1 {
            set(m, i, i, ncast(M.data_type, values[i]));
        }
    }
}


mat_init_hadamard_2 :: inline (m: *$M/MatrixType) {
    #if CHECKS {
        #run assert(M.N == 2 && M.M == 2);
    }
    mat_init(m, 1, 1, 1, -1);
}


// TODO come back to that later!
// mat_hadamard :: ($type: Type, $nestings: int) -> Matrix(type, 1<<nestings, 1<<nestings) {
//     #if CHECKS {
//         #run assert(nestings >= 0, "Nestings must be >=0 but = %", nestings);
//     }
//     #if nestings == 0 {
//         return mat_hadamard_2(type);
//     } else {
//         return tensor(mat_hadamard_2(type), mat_hadamard(type, nestings-1));
//     }
// }


mat_reflector :: (normal: $V/VectorType) -> DenseMatrix(V.data_type, V.N, V.N) #must {
    tau := ncast(T, 2)/(conjugate(normal)*normal);

    res := make(DenseMatrix(V.data_type, V.N, V.N), false);
    mat_init_unit(*res, true);
    sub(*res, tau * normal / normal);
    return res;
}






// inclusive first bound, exclusive last bound
mat_submatrix :: (
    m: $M/MatrixType, 
    $row_from: int, 
    $row_to: int, 
    $column_from: int, 
    $column_to: int
) -> MatrixSubmatrixHelper(M, row_to-row_from, column_to-column_from).T #must {
    #if CHECKS {
        #run assert(row_to > row_from, "\ncannot create submatrix with rows %->%\n", row_to, row_from);
        #run assert(column_to > column_from, "\ncannot create submatrix with columns %->%\n", column_to, column_from);
        #run assert(row_to <= M.N);
        #run assert(column_to <= M.M);
        #run assert(row_from >= 0);
        #run assert(column_from >= 0);
    }
    res := make(MatrixSubmatrixHelper(M, row_to-row_from, column_to-column_from).T, false);

    for i: 0..row_to-row_from-1 {
        for j: 0..column_to-column_from-1 {
            set(*res, i, j, get(m, i+row_from, j+column_from));
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


// TODO non-allocating versions of these functions!

row :: (m: $M/MatrixType, i: int) -> MatrixVectorHelper(M, M.M).T #must {
    #if CHECKS {
        assert(i >= 0 && i < M.N);
    }
    res := make(MatrixVectorHelper(M, M.M).T, false);
    for 0..M.M-1 {
        set(*res, it, get(m,i,it));
    }
    return res;
}

column :: (m: $M/MatrixType, i: int) -> MatrixVectorHelper(M, M.N).T #must {
    #if CHECKS {
        assert(i >= 0 && i < M.M);
    }
    res := make(MatrixVectorHelper(M, M.N).T, false);
    for 0..M.N-1 {
        set(*res,it, get(m,it,i));
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




operator == :: (a: $A/MatrixType, b: $B/MatrixType) -> bool {
    #if A.N != B.N || A.M != B.M {
        return false;
    } else {
        // TODO be smart about this, e.g. for sparse matrices
        for a {
            if get_lin(b, it_index) != it then return false;
        }
        return true;
    }
}



// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 2.11 Case Study: Matrix-Matrix Multiplication, page 118
// ? According to Table 2.12, the best performing loop order is j/k/i with
// ? the least number of cache misses.
// ? 
// ? Similar discussion in:
// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.5 ijk Matrix Multiplication, page 8

// * Assuming res is initialized to zeroes.
// * User has to know, what type the result has to be!
mat_mul_mat :: (a: $A/MatrixType, b: $B/MatrixType, res: *$C/MatrixType) {
    #if CHECKS {
        #run assert(A.M == B.N, "Matrix dimensions don't match: a.cols = %, b.rows = %", A.M, B.N);
        #run assert(A.N == C.N, "Matrix dimensions don't match: a.rows = %, c.rows = %", A.N, C.N);
        #run assert(B.M == C.M, "Matrix dimensions don't match: b.cols = %, c.cols = %", B.M, C.M);
    }

    for j: 0..B.M-1 {
        for k: 0..A.M-1 {
            tmp := get(b,k,j);
            for i: 0..A.N-1 {
                set(res,i,j, get(res,i,j) + get(a,i,k) * tmp);
            }
        }
    }
}

operator * :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMultiplicationHelper(A,B,A.N, B.M).T {
    res := make(MatrixMultiplicationHelper(A,B,A.N, B.M).T);
    inline mat_mul_mat(a,b,*res);
    return res;
}

// * Assuming res is initialized to zeroes.
// * User has to know, what type the result has to be!
mat_mul_vec :: (a: $A/MatrixType, b: $B/VectorType, res: *$C/VectorType) {
    #if CHECKS {
        #run assert(A.M == B.N, "Dimensions mismatch: a.cols = %, b.dim   = %", A.M, B.N);
        #run assert(A.N == C.N, "Dimensions mismatch: a.rows = %, res.dim = %", A.N, C.N);
    }

    for a {
        set(res, it_row, get(res, it_row) + it * get(b, it_column));
    }
}
operator * :: (a: $A/MatrixType, b: $B/VectorType) -> MatrixVectorMultiplicationHelper(A,B,A.N).T {
    res := make(MatrixVectorMultiplicationHelper(A,B,A.N).T);
    inline mat_mul_vec(a,b,*res);
    return res;
}

// * Assuming res is initialized to zeroes.
// * User has to know, what type the result has to be!
vec_mul_mat :: (a: $A/VectorType, b: $B/MatrixType, res: *$C/VectorType) {
    #if CHECKS {
        #run assert(A.N == B.N, "Dimensions mismatch: a.dim  = %, b.rows  = %", A.N, B.N);
        #run assert(B.M == C.N, "Dimensions mismatch: b.cols = %, res.dim = %", B.M, C.N);
    }

    for j: 0..B.N-1 {
        for i: 0..A.N-1 {
            set(res,j, get(res,j) + get(a,i) * get(b,i,j));
        }
    }
}
operator * :: (a: $A/VectorType, b: $B/MatrixType) -> MatrixVectorMultiplicationHelper(A,B,B.M).T {
    res := make(MatrixVectorMultiplicationHelper(A,B,B.M).T);
    inline vec_mul_mat(a,b,*res);
    return res;
}

// * User has to know, what type the result has to be!
mul :: (a: $A/MatrixType, b: $B, res: *$C/MatrixType) #modify {
    if !is_number(B) {B = null; return;}
} {
    #if CHECKS {
        #run assert(A.N == C.N, "Dimensions mismatch: a.rows  = %, res.rows  = %", A.N, C.N);
        #run assert(A.M == C.M, "Dimensions mismatch: a.cols  = %, res.cols  = %", A.M, C.M);
    }

    for a {
        set_lin(res, it_index, it * b);
    }
}
operator * :: (a: $A/MatrixType, b: $B) -> MatrixCastHelper(A, B).T #modify {
    if !is_number(B) {B = null; return;}    
} #symmetric {
    res := make(MatrixCastHelper(A, B).T, false);
    inline mul(a,b,*res);
    return res;
}


mul :: (a: *$A/MatrixType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    #if CHECKS {
        #run assert(is_minor_number_type(A.data_type, B));
    }

    for raw(a) {
        set_raw(a, it_raw_index, it * b);
    }
}





// mat_div_mat :: (a: $A/MatrixType, b: $B/MatrixType, res: *$C/MatrixType) {
// }
operator / :: inline (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMultiplicationHelper(A,B,A.M,B.M).T { 
    return inverse(b) * a;
}

operator / :: (a: $A/VectorType, b: $B/MatrixType) -> MatrixVectorMultiplicationHelper(B, A, A.N).T { 
    return inverse(b) * a;
}

// * User has to know, what type the result has to be!
div :: (a: $A/MatrixType, b: $B, res: *$C/MatrixType) #modify {
    if !is_number(B) { B = null; return; }
} {
    for a {
        set_lin(res, it_index, it / b);
    }
}
operator / :: inline (a: $A/MatrixType, b: $B) -> MatrixCastHelper(A, B).T #modify {
    if !is_number(B) {B = null; return;}
} { 
    res := make(MatrixCastHelper(A, B).T, false);
    inline div(a,b,*res);
    return res;
}



div :: (a: *$A/MatrixType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    #if CHECKS {
        #run assert(is_minor_number_type(A.data_type, B));
    }

    for raw(a) {
        set_raw(a, it_raw_index, it / b);
    }
}


// * User has to know, what type the result has to be!
mat_add_mat :: (a: $A/MatrixType, b: $B/MatrixType, res: *$C/MatrixType) {
    #if CHECKS {
        #run assert(A.N == B.N, "Matrix dimensions don't match: a.rows = %, b.rows = %", A.N, B.N);
        #run assert(A.N == C.N, "Matrix dimensions don't match: a.rows = %, b.rows = %", A.N, C.N);
        #run assert(A.M == B.M, "Matrix dimensions don't match: a.cols = %, c.cols = %", A.M, B.M);
        #run assert(A.M == C.M, "Matrix dimensions don't match: a.cols = %, c.cols = %", A.M, C.M);
    }

    for a {
        set_lin(res, it_index, it + get_lin(b, it_index));
    }
}
operator + :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMatrixCastHelper(A,B).T {
    res := make(MatrixMatrixCastHelper(A,B).T, false);
    inline mat_add_mat(a,b,*res);
    return res;
}





add :: (a: *$A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(A.N == B.N && A.M == B.M, "matrix dimensions do not match A(% x %), B(% x %)", A.N, A.M, B.N, B.M);
        #run assert(is_minor_number_type(A.data_type, B.data_type), "cannot cast % to %.", B.data_type, A.data_type);
    }

    for a {
        set_lin(a, it_index, it + get_lin(b, it_index));
    }
}







// * User has to know, what type the result has to be!
mat_sub_mat :: (a: $A/MatrixType, b: $B/MatrixType, res: *$C/MatrixType) {
    #if CHECKS {
        #run assert(A.N == B.N, "Matrix dimensions don't match: a.rows = %, b.rows = %", A.N, B.N);
        #run assert(A.N == C.N, "Matrix dimensions don't match: a.rows = %, b.rows = %", A.N, C.N);
        #run assert(A.M == B.M, "Matrix dimensions don't match: a.cols = %, c.cols = %", A.M, B.M);
        #run assert(A.M == C.M, "Matrix dimensions don't match: a.cols = %, c.cols = %", A.M, C.M);
    }

    for a {
        set_lin(res, it_index, it - get_lin(b, it_index));
    }
}
operator - :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMatrixCastHelper(A,B).T {
    res := make(MatrixMatrixCastHelper(A,B).T, false);
    inline mat_sub_mat(a,b,*res);
    return res;
}


// TODO check for "UpCast" of Matrices
sub :: (a: *$A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(A.N == B.N && A.M == B.M, "matrix dimensions do not match A(% x %), B(% x %)", A.N, A.M, B.N, B.M);
        #run assert(is_minor_number_type(A.data_type, B.data_type), "cannot cast % to %.", B.data_type, A.data_type);
    }

    for a {
        set_lin(a, it_index, it - get_lin(b, it_index));
    }
}





operator - :: (a: $A/MatrixType) -> A {
    res := make(A, false);

    for a {
        set_lin(*res, it_index, -it);
    }

    return res;
}



// ? by iterating through the raw data, there is no double-negation issue!
neg :: (a: *$A/MatrixType) {
    for raw(a) {
        set_raw(a, it_raw_index, -it);
    }
}















// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.6 Dot Product Matrix Multiplication, page 10
// ? Algorithm 1.1.7 Saxpy Matrix Multiplication, page 10
// ? C = C + A*B
saxpy :: (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(A.N == C.N);
        #run assert(B.N == A.M && B.M == C.M);
        #run assert(is_minor_number_type(C.data_type, A.data_type) && is_minor_number_type(C.data_type, B.data_type));
    }

    for j: 0..C.M-1 {
        for k: 0..A.M-1 {
            tmp := get(b,k,j);
            for i: 0..C.N-1 {
                set(c,i,j, get(c,i,j) + get(a,i,k) * tmp);
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
// check for compatible matrix types: e.g. this is not possible with triangular matrices
permute_rows :: (m: $M/MatrixType, order: ..int) -> M #must {
    #if CHECKS {
        #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(order.count == N);
    }

    res := make(M, false);

    for res {
        set(*res, it_row, it_column, get(m, order[it_row], it_column));
    }

    return res;
}

// there are faster ways of doing this.
permute_columns :: (m: $M/MatrixType, order: ..int) -> M #must {
    #if CHECKS {
        #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(order.count == M.M);
    }

    res := make(M, false);

    for res {
        set(*res, it_row, it_column, get(m, it_row, order[it_column]));
    }

    return res;
}

// there are faster ways of doing this.
permute :: (m: $M/MatrixType, row_order: [M.N]int, col_order: [M.M]int) -> M #must {
    #if CHECKS {
        #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
    }

    res := make(M, false);

    for res {
        set(*res, it_row, it_column, get(m, row_order[it_row], col_order[it_column]));
    }

    return res;
}


swap_columns :: (m: *$M/MatrixType, i: int, j: int) {
    #if CHECKS {
        #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(i >= 0 && j >= 0);
        assert(i < M.M && j < M.M);
    }

    if i == j then return;

    for r: 0..M.N-1 {
        tmp := get(m, r, i);
        set(m, r, i, get(m, r, j));
        set(m, r, j, tmp);
    }
}

swap_rows :: (m: *$M/MatrixType, i: int, j: int) {
    #if CHECKS {
        #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense) : %\n", M.flags);
        assert(i >= 0 && j >= 0);
        assert(i < M.N && j < M.N);
    }

    if i == j then return;

    for c: 0..M.M-1 {
        tmp := get(m, i, c);
        set(m, i, c, get(m, j, c));
        set(m, j, c, tmp);
    }
}


transpose :: (m: $M/MatrixType) -> MatrixTransposeHelper(M).T #must {
    #if CHECKS {
        #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)");
    }

    // println_push("transpose", color = .FG_BLUE);

    res := make(MatrixTransposeHelper(M).T);

    for m {
        set(*res, it_column, it_row, it);
    }

    // println("res = %", pstr(m));
    return res;
}

// There are better ways to do this!
transpose :: (m: *$M/MatrixType) {
    #if CHECKS {
        #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)");
    }

    for i: 0..M.N-1 {
        for j: i+1..M.M-1 {
            tmp := get(m, i, j);
            set(m, i, j, get(m, j, i));
            set(m, j, i, tmp);
        }
    }
}





conjugate :: (m: $M/MatrixType) -> M #must {
    res := make(M, false);
    #if #run !is_non_real_number(M.data_type) {
        copy(*res, m);
    } else {
        for raw(*m) {
            set_raw(*res, it_raw_index, conjugate(it));
        }
    }
    return res;
}

conjugate :: inline (m: *$M/MatrixType) {
    #if #run !is_non_real_number(M.data_type) {
        return;
    } else {
        for raw(m) {
            set_raw(m, it_raw_index, conjugate(it));
        }
    }
}


conjugate_transpose :: inline (m: $M/MatrixType) -> MatrixTransposeHelper(M).T #must {
    return conjugate(transpose(m));
}
conjugate_transpose :: inline (m: *$M/MatrixType) {
    transpose(m);
    conjugate(m);
}
dagger :: conjugate_transpose;
// hermitian :: conjugate_transpose;





// https://www.wikiwand.com/en/Tensor_product
tensor :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixHelper(A, B, A.N*B.N, A.M*B.M).T {
    res := make(MatrixHelper(A, B, A.N*B.N, A.M*B.M).T, false);

    for r: 0..res.M-1 {
        ar := r/B.N;
        br := r % B.N;
        
        for c: 0..res.M-1 {
            ac := c/B.M;
            bc := c % B.M;

            set(*res, r, c, get(a, ar, ac) * get(b, br, bc));
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



// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Lemma 3.5.4 Matrix Norms, page 176

// this needs the largest eigenvalue
// norm_2 :: inline (v: Vector) -> float64 {
//     return norm(v, 2.0);
// }

norm_1 :: (m: $M/MatrixType) -> float64 {
    res : float64 = 0.0;

    for j: 0..M.M-1 {
        sum : float64 = 0.0;
        for i: 0..M.N-1 {
            sum += abs(get(m,i,j));
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}

norm_inf :: (m: $M/MatrixType) -> float64 {
    res : float64 = 0.0;

    for i: 0..M.N-1 {
        sum : float64 = 0.0;
        for j: 0..M.M-1 {
            sum += abs(get(m,i,j));
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}

norm_frobenius :: (m: $M/MatrixType) -> float64 {
    sum : float64 = 0.0;

    for m {
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

    println_push("Matrix", color = .FG_WHITE);

    {
        println_push("DenseMatrix", color = .FG_GREEN);
        {
            m := make(DenseMatrix(float64, 3, 3));
            println("%", pstr(m));
        }
        {
            a :: 3;
            m := make(DenseMatrix(float64, a, a));
            println("%", pstr(m));
        }



        {
            m := make(DenseMatrix(float64, 3, 3));
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }
        {
            a :: 3;
            m := make(DenseMatrix(float64, a, a));
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }


        {
            m := make(DenseMatrix(C64, 3, 3));
            mat_init_all(*m, cpx(42.0));
            set(*m,1,2,R);
            assert( get(m,1,2) == cpx(1.0, 0.0) );
            println("m = %", pstr(m));
        }

        {
            m := make(DenseMatrix(C64, 4, 4));
            mat_init_unit(*m);
            println("m = %", pstr(m));
            assert(get(m,0,1) == C_ZERO);
            assert(get(m,1,1) == R);
            assert(get(m,2,1) == C_ZERO);
            assert(get(m,3,3) == R);
        }

        {
            m := make(DenseMatrix(C64, 3, 3));
            set(*m,0,2,R);
            m = transpose(m);
            assert(get(m,2,0) == R);
        }

        {
            a := make(DenseMatrix(float64, 2, 2));
            b := make(DenseMatrix(float64, 2, 2));
            mat_init_unit(*a);
            mat_init_unit(*b);
            println("a = %", pstr(a));
            println("b = %", pstr(b));
            add(*a, b);
            println("a = %", pstr(a));
            println("2 b = %", pstr(2.0 * b));
            assert(a == 2.0 * b);

            a = a + b;
            assert(a == 3.0 * b);
        }

        {
            m := make(DenseMatrix(C64, 2, 2));
            mat_init(*m, R, 2*R, 3*R, 4*R);
            println("m = %", pstr(m));
            m = m * m;
            println("m = %", pstr(m));
        }
        {
            m := make(DenseMatrix(C64, 2, 2));
            mat_init(*m, R, 2.0*R, 3.0*R, 4.0*R);
            println("m = %", pstr(m));
            m = m * m;
            println("m = %", pstr(m));
        }
        // {
        //     m := mat(1, 3, 1.0, 2.0, 3.0);
        //     print("m = %\n", m);
        //     mt := transpose(m);
        //     m = m * mt;
        //     print("m = %\n", m);
        // }
        {
            m := make(DenseMatrix(C64, 2, 2));
            mat_init(*m, 
                C_ZERO, -I, 
                I, C_ZERO
            );
            assert( dagger(m) == m );
        }


        {
            m := make(DenseMatrix(C64, 2, 2));
            mat_init(*m, R, I, 2*R-I, 2*R + 3*I);
            println("m = %", pstr(m));
            t1 := tensor(m, m);
            println("t1 = %", pstr(t1));
            t2 := make(DenseMatrix(C64, 4, 4));
            mat_init(*t2,
                R            ,             I,             I,   -R           ,
                2.0 -       I, 2.0 + 3.0 * I, 1.0 + 2.0 * I, -3.0 +  2.0 * I,
                2.0 -       I, 1.0 + 2.0 * I, 2.0 + 3.0 * I, -3.0 +  2.0 * I,
                3.0 - 4.0 * I, 7.0 + 4.0 * I, 7.0 + 4.0 * I, -5.0 + 12.0 * I
            );
            println("t2 = %", pstr(t2));
            assert(t1 == t2);
        }

        {
            M := make(DenseMatrix(C64, 3, 3));
            mat_init(*M,
                C_ZERO, R, 2*R,
                3*R, 4*R, 5*R,
                6*R, 7*R, 8*R
            );

            S := mat_submatrix(M, 1, 3, 1, 3);
            
            res := make(DenseMatrix(C64, 2, 2));
            mat_init(*res, 4*R, 5*R, 7*R, 8*R);
            
            assert(S == res);
        }

        // {
        //     print("%\n", MATRIX_FLAGS.UNIT_LOWER_TRIANGULAR);
        //     m := make(DenseMatrix(float64, 4, 4)) = ---;
        //     m.flags = .UNIT_LOWER_TRIANGULAR;
        //     assert( is(m, .TRAPEZOIDAL) );
        // }
    }



    {
        println_push("DenseHeapMatrix", color = .FG_GREEN);
        {
            println_push("test 1", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(float64, 3, 3));
            defer free(m);
            println("%", pstr(m));
        }
        {
            println_push("test 2", color = .FG_YELLOW);
            a :: 3;
            m := make(DenseHeapMatrix(float64, a, a));
            defer free(m);
            println("%", pstr(m));
        }



        {
            println_push("test 3", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(float64, 3, 3));
            defer free(m);
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }
        {
            println_push("test 4", color = .FG_YELLOW);
            a :: 3;
            m := make(DenseHeapMatrix(float64, a, a));
            defer free(m);
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }


        {
            println_push("test 5", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(C64, 3, 3));
            defer free(m);
            mat_init_all(*m, cpx(42.0));
            set(*m,1,2,R);
            assert( get(m,1,2) == cpx(1.0, 0.0) );
            println("m = %", pstr(m));
        }

        {
            println_push("test 6", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(C64, 4, 4));
            defer free(m);
            mat_init_unit(*m);
            println("m = %", pstr(m));
            assert(get(m,0,1) == C_ZERO);
            assert(get(m,1,1) == R);
            assert(get(m,2,1) == C_ZERO);
            assert(get(m,3,3) == R);
        }

        {
            println_push("test 7", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(C64, 3, 3));
            defer free(m);
            set(*m,0,2,R);
            println("transposing: ");
            m = transpose(m);
            assert(get(m,2,0) == R);
        }

        {
            println_push("test 8", color = .FG_YELLOW);
            a := make(DenseHeapMatrix(float64, 2, 2));
            defer free(a);
            b := make(DenseHeapMatrix(float64, 2, 2));
            defer free(b);
            mat_init_unit(*a);
            mat_init_unit(*b);
            println("a = %", pstr(a));
            println("b = %", pstr(b));
            add(*a, b);
            println("a = %", pstr(a));
            println("2 b = %", pstr(2.0 * b));
            assert(a == 2.0 * b);

            a = a + b;
            assert(a == 3.0 * b);
        }

        {
            println_push("test 9", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(C64, 2, 2));
            defer free(m);
            mat_init(*m, R, 2*R, 3*R, 4*R);
            println("m = %", pstr(m));
            m = m * m;
            println("m = %", pstr(m));
        }
        {
            println_push("test 10", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(C64, 2, 2));
            defer free(m);
            mat_init(*m, R, 2.0*R, 3.0*R, 4.0*R);
            println("m = %", pstr(m));
            m = m * m;
            println("m = %", pstr(m));
        }
        {
            println_push("test 11", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(C64, 2, 2));
            defer free(m);
            mat_init(*m, 
                C_ZERO, -I, 
                I, C_ZERO
            );
            assert( dagger(m) == m );
        }


        {
            println_push("test 12", color = .FG_YELLOW);
            m := make(DenseHeapMatrix(C64, 2, 2));
            defer free(m);
            mat_init(*m, R, I, 2*R-I, 2*R + 3*I);
            println("m = %", pstr(m));
            t1 := tensor(m, m);
            println("t1 = %", pstr(t1));
            t2 := make(DenseHeapMatrix(C64, 4, 4));
            mat_init(*t2,
                R            ,             I,             I,   -R           ,
                2.0 -       I, 2.0 + 3.0 * I, 1.0 + 2.0 * I, -3.0 +  2.0 * I,
                2.0 -       I, 1.0 + 2.0 * I, 2.0 + 3.0 * I, -3.0 +  2.0 * I,
                3.0 - 4.0 * I, 7.0 + 4.0 * I, 7.0 + 4.0 * I, -5.0 + 12.0 * I
            );
            println("t2 = %", pstr(t2));
            assert(t1 == t2);
        }

        {
            println_push("test 13", color = .FG_YELLOW);
            M := make(DenseHeapMatrix(C64, 3, 3));
            defer free(M);
            mat_init(*M,
                C_ZERO, R, 2*R,
                3*R, 4*R, 5*R,
                6*R, 7*R, 8*R
            );

            S := mat_submatrix(M, 1, 3, 1, 3);
            
            res := make(DenseHeapMatrix(C64, 2, 2));
            mat_init(*res, 4*R, 5*R, 7*R, 8*R);
            
            assert(S == res);
        }
    }
}