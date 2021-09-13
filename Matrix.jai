






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
    return str(m, p_offset, p_shift + 4, "|");
}

str :: (m: $M/MatrixType, offset: int = 0, indent: int = 4, indent_character: string = " ", indent_every: int = 4) -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    append(*builder, "\n");
    for i: 0..M.N-1 {
        for 0..offset-1 {
            append(*builder, " ");
        }
        for 0..indent-1 {
            if it % indent_every == 0 {
                append(*builder, indent_character);
            } else {
                append(*builder, " ");
            }
        }
        for j: 0..M.M-1 {
            #if M.data_type == Complex64 {
                print_to_builder(*builder, "%, ", str(get(m,i,j)) );
            } else {
                print_to_builder(*builder, "%, ", get(m,i,j) );
            }
        }
        if i != M.N-1 then
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
    res : MatrixCastHelper(M, type).T = ---;

    for i: 0..M.N-1 {
        for j: 0..M.M-1 {
            set(*res, i, j, ncast(type, get(m, i, j)));
        }
    }

    return res;
}

copy :: (m: $M/MatrixType) -> M #must {
    res : M = ---;
    for i: 0..M.N-1 {
        for j: 0..M.M-1 {
            set(*res, i, j, get(m, i, j));
        }
    }
    return res;
}

copy :: inline (to: *$M/MatrixType, from: M) {
    for i: 0..M.N-1 {
        for j: 0..M.M-1 {
            set(to, i, j, get(from, i, j));
        }
    }
}

mat_init :: (m: *$M/MatrixType, vals: ..$T) {
    #if CHECKS {
        #run assert(is_minor_number(M.data_type, T));
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
    for i: 0..M.N-1 {
        for j: 0..M.M-1 {
            set(m, i, j, ncast(M.data_type, value));
        }
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
        for i: 0..M.N-1 {
            for j: 0..M.M-1 {
                if i == j { 
                    set(m, i, j, ncast(M.data_type, values[i]));
                } else {
                    set(m, i, j, zero(M.data_type));
                }
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

    res : DenseMatrix(V.data_type, V.N, V.N) = ---;
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
    res : MatrixSubmatrixHelper(M, row_to-row_from, column_to-column_from).T = ---;

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


row :: (m: $M/MatrixType, i: int) -> MatrixVectorHelper(M, M.M).T #must {
    #if CHECKS {
        assert(i >= 0 && i < M.N);
    }
    res : MatrixVectorHelper(M, M.M).T = ---;
    for 0..M.M-1 {
        set(*res, it, get(m,i,it));
    }
    return res;
}

column :: (m: $M/MatrixType, i: int) -> MatrixVectorHelper(M, M.N).T #must {
    #if CHECKS {
        assert(i >= 0 && i < M.M);
    }
    res : MatrixVectorHelper(M, M.N).T = ---;
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
        for i: 0..A.N-1 {
            for j: 0..A.M-1 {
                if get(a, i, j) != get(b, i, j) then return false;
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
operator * :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMultiplicationHelper(A,B,A.N, B.M).T {
    #if CHECKS {
        #run assert(A.M == B.N, "Matrix dimensions don't match: a.cols = %, b.rows = %", A.M, B.N);
    }

    res : MatrixMultiplicationHelper(A,B,A.N, B.M).T;
    for j: 0..B.M-1 {
        for k: 0..A.M-1 {
            tmp := get(b,k,j);
            for i: 0..A.N-1 {
                set(*res,i,j, get(res,i,j) + get(a,i,k) * tmp);
            }
        }
    }

    return res;
}


#scope_file 
matrix_vector_multiplication_helper :: ($a: Type, $b: Type, $N: int) -> Type {
    // TODO other types
    return DenseVector(UpCast(a.data_type, b.data_type).T, N);
}
#scope_export
MatrixVectorMultiplicationHelper :: struct (A: Type, B: Type, N: int) {
    T :: #run matrix_vector_multiplication_helper(A,B,N);
}
operator * :: (a: $A/MatrixType, b: $B/VectorType) -> MatrixVectorMultiplicationHelper(A,B,A.N).T {
    #if CHECKS {
        #run assert(A.M == B.N, "Dimensions mismatch: a.cols = %, b.dim = %", A.M, B.N);
    }

    res : MatrixVectorMultiplicationHelper(A,B,A.N).T; // ! need to initialize values to 0!!!!

    for a {
        set(*res, it_row, get(res, it_row) + it * get(b, it_column));
    }

    return res;
}

operator * :: (a: $A/VectorType, b: $B/MatrixType) -> MatrixVectorMultiplicationHelper(A,B,B.M).T {
    #if CHECKS {
        #run assert(A.N == B.N, "Dimensions mismatch: a.dim = %, b.rows = %", A.N, B.N);
    }

    res : MatrixVectorMultiplicationHelper(A,B,B.M).T;
    for j: 0..B.N-1 {
        for i: 0..A.N-1 {
            set(*res,j, get(res,j) + get(a,i) * get(b,i,j));
        }
    }

    return res;
}


operator * :: (a: $A/MatrixType, b: Scalar($B)) -> MatrixCastHelper(A, B).T #symmetric {
    res : MatrixCastHelper(A, B).T = ---;

    for a {
        set_lin(*res, it_index, it * b.value);
    }

    return res;
}


mul :: (a: *$A/MatrixType, b: Scalar($B)) {
    #if CHECKS {
        #run assert(is_minor_number(A.data_type, B));
    }

    for i: 0..A.N-1 {
        for j: 0..A.M-1 {
            set(a, i, j, get(<<a, i, j) * b.value);
        }
    }
}






operator / :: inline (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMultiplicationHelper(A,B,A.M,B.M).T { 
    return inverse(b) * a;
}

operator / :: (a: $A/VectorType, b: $B/MatrixType) -> MatrixVectorMultiplicationHelper(B, A, A.N).T { 
    return inverse(b) * a;
}

operator / :: inline (a: $A/MatrixType, b: Scalar($B)) -> MatrixCastHelper(A, B).T { 
    res : MatrixCastHelper(A, B).T = ---;

    for a {
        set(*res, it_row, it_column, it / b.value);
    }

    return res;
}



div :: (a: *$A/MatrixType, b: Scalar($B)) {
    #if CHECKS {
        #run assert(is_minor_number(A.data_type, B));
    }

    for i: 0..A.N-1 {
        for j: 0..A.M-1 {
            set(a, i, j, get(<<a, i, j) / b.value);
        }
    }
}




operator + :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMatrixCastHelper(A,B).T {
    #if CHECKS {
        #run assert(A.N == B.N && A.M == B.M, "matrix dimensions do not match A(% x %), B(% x %)", A.N, A.M, B.N, B.M);
    }
    
    res : MatrixMatrixCastHelper(A,B).T = ---;

    for i: 0..A.N-1 {
        for j: 0..A.M-1 {
            set(*res, i, j, get(a, i, j) + get(b, i, j));
        }
    }
    return res;
}





add :: (a: *$A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(A.N == B.N && A.M == B.M, "matrix dimensions do not match A(% x %), B(% x %)", A.N, A.M, B.N, B.M);
        #run assert(is_minor_number(A.data_type, B.data_type), "cannot cast % to %.", B.data_type, A.data_type);
    }

    for i: 0..A.N-1 {
        for j: 0..A.M-1 {
            set(a, i, j, get(<<a, i, j) + get(b, i, j));
        }
    }
}








operator - :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMatrixCastHelper(A,B).T {
    #if CHECKS {
        #run assert(A.N == B.N && A.M == B.M, "matrix dimensions do not match A(% x %), B(% x %)", A.N, A.M, B.N, B.M);
    }
    
    res : MatrixMatrixCastHelper(A,B).T = ---;

    for i: 0..A.N-1 {
        for j: 0..A.M-1 {
            set(*res, i, j, get(a, i, j) - get(b, i, j));
        }
    }
    return res;
}


// TODO check for "UpCast" of Matrices
sub :: (a: *$A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(A.N == B.N && A.M == B.M, "matrix dimensions do not match A(% x %), B(% x %)", A.N, A.M, B.N, B.M);
        #run assert(is_minor_number(A.data_type, B.data_type), "cannot cast % to %.", B.data_type, A.data_type);
    }

    for i: 0..A.N-1 {
        for j: 0..A.M-1 {
            set(a, i, j, get(<<a, i, j) - get(b, i, j));
        }
    }
}





operator - :: (a: $A/MatrixType) -> A {
    res : A = ---;

    for a {
        set_lin(*res, it_index, -it);
    }

    return res;
}



// TODO be careful with e.g. symmetric matrices:
// if I first set Aij = -Aji and then Aji = -Aij = -(-Aji) = Aji
// so there's no negation happening at all!
neg :: (a: *$A/MatrixType) {
    for a {
        set_lin(a, it_index, -it);
    }
}















// Matrix Computations, 4th edition
// The Johns Hopkins University Press
// Algorithm 1.1.6 Dot Product Matrix Multiplication, page 10
// Algorithm 1.1.7 Saxpy Matrix Multiplication, page 10
// C = C + A*B
saxpy :: (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(A.N == C.N);
        #run assert(B.N == A.M && B.M == C.M);
        #run assert(is_minor_number(C.data_type, A.data_type) && is_minor_number(C.data_type, B.data_type));
    }

    for j: 0..C.M-1 {
        for k: 0..A.M-1 {
            tmp := get(b,k,j);
            for i: 0..C.N-1 {
                set(c,i,j, get(<<c,i,j) + get(a,i,k) * tmp);
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
        #run assert(M.flags == .DENSE, "not a compatible matrix type (not dense)\n");
        assert(order.count == N);
    }

    res : M = ---;

    for res {
        set(*res, it_row, it_column, get(m, order[it_row], it_column));
    }

    return res;
}

// there are faster ways of doing this.
permute_columns :: (m: $M/MatrixType, order: ..int) -> M #must {
    #if CHECKS {
        #run assert(M.flags == .DENSE, "not a compatible matrix type (not dense)\n");
        assert(order.count == M.M);
    }

    res : M = ---;

    for res {
        set(*res, it_row, it_column, get(m, it_row, order[it_column]));
    }

    return res;
}

// there are faster ways of doing this.
permute :: (m: $M/MatrixType, row_order: [M.N]int, col_order: [M.M]int) -> M #must {
    #if CHECKS {
        #run assert(M.flags == .DENSE, "not a compatible matrix type (not dense)\n");
    }

    res : M = ---;

    for res {
        set(*res, it_row, it_column, get(m, row_order[it_row], col_order[it_column]));
    }

    return res;
}


swap_columns :: (m: *$M/MatrixType, i: int, j: int) {
    #if CHECKS {
        #run assert(M.flags == .DENSE, "not a compatible matrix type (not dense)\n");
        assert(i >= 0 && j >= 0);
        assert(i < M.M && j < M.M);
    }

    if i == j then return;

    for r: 0..M.N-1 {
        tmp := get(<<m, r, i);
        set(m, r, i, get(<<m, r, j));
        set(m, r, j, tmp);
    }
}

swap_rows :: (m: *$M/MatrixType, i: int, j: int) {
    #if CHECKS {
        #run assert(M.flags == .DENSE, "not a compatible matrix type (not dense) : %\n", M.flags);
        assert(i >= 0 && j >= 0);
        assert(i < M.N && j < M.N);
    }

    if i == j then return;

    for c: 0..M.M-1 {
        tmp := get(<<m, i, c);
        set(m, i, c, get(<<m, j, c));
        set(m, j, c, tmp);
    }
}



#scope_file 
MatrixTransposeHelper :: struct(A: Type) {
    T :: #run matrix_transpose_helper(A);
}
matrix_transpose_helper :: ($a: Type) -> Type {
    return DenseMatrix(a.data_type, a.M, a.N);
}
#scope_export


transpose :: (m: $M/MatrixType) -> MatrixTransposeHelper(M).T #must {
    #if CHECKS {
        #run assert(M.flags == .DENSE, "not a compatible matrix type (not dense)");
    }

    res : MatrixTransposeHelper(M).T = ---;

    for m {
        set(*res, it_column, it_row, it);
    }
    return res;
}

// There are better ways to do this!
transpose :: (m: *$M/MatrixType) {
    #if CHECKS {
        #run assert(M.flags == .DENSE, "not a compatible matrix type (not dense)");
    }

    for i: 0..M.N-1 {
        for j: i+1..M.M-1 {
            tmp := get(<<m, i, j);
            set(m, i, j, get(<<m, j, i));
            set(m, j, i, tmp);
        }
    }
}





conjugate :: (m: $M/MatrixType) -> M #must {
    #if M.data_type != Complex64 {
        return m;
    } else {
        res : M = ---;

        for m {
            set_lin(*res, it_index, conjugate(it));
        }

        return res;
    }
}

conjugate :: inline (m: *$M/MatrixType) {
    #if M.data_type != Complex64 {
        return;
    } else {
        for m {
            set(m, it_index, conjugate(it));
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
    res : MatrixHelper(A, B, A.N*B.N, A.M*B.M).T = ---;

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



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Lemma 3.5.4 Matrix Norms, page 176

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

    for i: 0..M.N-1 {
        for j: 0..M.M-1 {
            sum += abs(get(m, i, j));
        }
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
    print("\n\nDenseMatrix(T) -------------------------------\n\n");
    {
        m : DenseMatrix(float64, 3, 3);
        print("%\n", str(m));
    }
    {
        a :: 3;
        m : DenseMatrix(float64, a, a);
        print("%\n", str(m));
    }



    {
        m : DenseMatrix(float64, 3, 3);
        mat_init_all(*m, 42.0);
        print("%\n", str(m));
    }
    {
        a :: 3;
        m : DenseMatrix(float64, a, a);
        mat_init_all(*m, 42.0);
        print("%\n", str(m));
    }


    {
        m : DenseMatrix(Complex64, 3, 3);
        mat_init_all(*m, c64(42.0));
        set(*m,1,2,R);
        assert( get(m,1,2) == c64(1.0, 0.0) );
        print("m = \n%\n", str(m));
    }

    {
        m : DenseMatrix(Complex64, 4, 4);
        mat_init_unit(*m);
        print("m = \n%\n", str(m));
        assert(get(m,0,1) == CZERO);
        assert(get(m,1,1) == R);
        assert(get(m,2,1) == CZERO);
        assert(get(m,3,3) == R);
    }

    {
        m : DenseMatrix(Complex64, 3, 3);
        set(*m,0,2,R);
        m = transpose(m);
        assert(get(m,2,0) == R);
    }

    {
        a : DenseMatrix(float64, 2, 2);
        b : DenseMatrix(float64, 2, 2);
        mat_init_unit(*a);
        mat_init_unit(*b);
        print("a = \n%\n", str(a));
        print("b = \n%\n", str(b));
        add(*a, b);
        print("a = \n%\n", str(a));
        print("2 b = \n%\n", scl(2.0) * b);
        assert(a == scl(2.0) * b);

        a = a + b;
        assert(a == scl(3.0) * b);
    }

    {
        m : DenseMatrix(Complex64, 2, 2);
        mat_init(*m, R, c64(2.0), c64(3.0), c64(4.0));
        print("m = \n%\n", str(m));
        m = m * m;
        print("m = \n%\n", str(m));
    }
    {
        m : DenseMatrix(Complex64, 2, 2);
        mat_init(*m, R, 2.0*R, 3.0*R, 4.0*R);
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
        m : DenseMatrix(Complex64, 2, 2);
        mat_init(*m, 
            CZERO, -I, 
            I, CZERO
        );
        assert( dagger(m) == m );
    }


    {
        m : DenseMatrix(Complex64, 2, 2);
        mat_init(*m, R, I, c64(2.0, -1.0), 2.0 + 3.0*I);
        print("m = \n%\n", str(m));
        t1 := tensor(m, m);
        print("t1 = \n%\n", str(t1));
        t2 : DenseMatrix(Complex64, 4, 4);
        mat_init(*t2,
            R            ,             I,             I,   -R           ,
            2.0 -       I, 2.0 + 3.0 * I, 1.0 + 2.0 * I, -3.0 +  2.0 * I,
            2.0 -       I, 1.0 + 2.0 * I, 2.0 + 3.0 * I, -3.0 +  2.0 * I,
            3.0 - 4.0 * I, 7.0 + 4.0 * I, 7.0 + 4.0 * I, -5.0 + 12.0 * I
        );
        print("t2 = \n%\n", str(t2));
        assert(t1 == t2);
    }

    {
        M : DenseMatrix(Complex64, 3, 3);
        mat_init(*M,
            CZERO, R, 2*R,
            3*R, 4*R, 5*R,
            6*R, 7*R, 8*R
        );

        S := mat_submatrix(M, 1, 3, 1, 3);
        
        res : DenseMatrix(Complex64, 2, 2);
        mat_init(*res, 4*R, 5*R, 7*R, 8*R);
        
        assert(S == res);
    }

    // {
    //     print("%\n", MATRIX_FLAGS.UNIT_LOWER_TRIANGULAR);
    //     m : DenseMatrix(float64, 4, 4) = ---;
    //     m.flags = .UNIT_LOWER_TRIANGULAR;
    //     assert( is(m, .TRAPEZOIDAL) );
    // }
}