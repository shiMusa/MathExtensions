






/*

 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##

*/

str :: (m: $M/MatrixType, indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    append(*builder, "\n");
    for i: 0..M.N-1 {
        append(*builder, indent);
        for j: 0..M.M-1 {
            #if T == Complex64 {
                print_to_builder(*builder, "%, ", str(get(m,i,j)) );
            } else {
                print_to_builder(*builder, "%, ", get(m,i,j) );
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

    for i: M.N {
        for j: M.M {
            set(*res, i, j, ncast(type, get(m, i, j)));
        }
    }

    return res;
}

copy :: (m: $M/MatrixType) -> M #must {
    res : M = ---;
    for i: 0..M.N {
        for j: 0..M.M {
            set(*res, i, j, get(m, i, j));
        }
    }
    return res;
}

copy :: inline (to: *M/MatrixType, from: M) {
    for i: 0..M.N {
        for j: 0..M.M {
            set(to, i, j, get(from, i, j));
        }
    }
}

mat_init :: (m: *$M/MatrixType, vals: ..$T) {
    #if CHECKS {
        assert(vals.count <= M.N * M.M, "too many values (%) for matrix of % entries\n", vals.count, M.N * M.M);
    }
    for vals {
        j := it_index % M.M;
        i := (it_index - j) / M.M;
        set(m, i, j, ncast(M.data_type, it));
    }
}

mat_init_all :: (m: *$M/MatrixType, value: $T) {
    for i: 0..M.N {
        for j: 0..M.M {
            set(m, i, j, ncast(M.data_type, value));
        }
    }
}

mat_init_unit :: (m: *$M/MatrixType, $init_zeroes: bool = false) {
    #if init_zeroes {
        for i: 0..M.N {
            for j: 0..M.M {
                if i == j { 
                    set(m, i, j, one(M.data_type));
                } else {
                    set(m, i, j, zero(M.data_type));
                }
            }
        }
    } else {
        k := min(M.N, M.M);
        for i: 0..k {
            set(m, i, i, one(M.data_type));
        }
    }
}

mat_init_diag :: (m: *$M/MatrixType, values: ..$T, $init_zeroes: bool = false) {
    #if CHECKS {
        assert(values.count == min(M.N, M.M), "not enough values (% values) for diagonal (% values)\n", values.count, min(M.N, M.M));
    }
    #if init_zeroes {
        for i: 0..M.N {
            for j: 0..M.M {
                if i == j { 
                    set(m, i, j, ncast(M.data_type, values[i]));
                } else {
                    set(m, i, j, zero(M.data_type));
                }
            }
        }
    } else {
        for i: 0..#run min(M.N, M.M) {
            set(m, i, i, ncast(M.data_type, values[i]));
        }
    }
}


mat_init_hadamard_2 :: inline (m: *$M/MatrixType) {
    #if CHECKS {
        #assert(M.N == 2 && M.M == 2);
    }
    mat_init(m, 1, 1, 1, -1);
}


// TODO come back to that later!
// mat_hadamard :: ($type: Type, $nestings: int) -> Matrix(type, 1<<nestings, 1<<nestings) {
//     #if CHECKS {
//         #assert(nestings >= 0, "Nestings must be >=0 but = %", nestings);
//     }
//     #if nestings == 0 {
//         return mat_hadamard_2(type);
//     } else {
//         return tensor(mat_hadamard_2(type), mat_hadamard(type, nestings-1));
//     }
// }


// TODO Vector(...) -> $V/VectorType
mat_reflector :: (normal: Vector($T, $N)) -> DenseMatrix(T, N, N) #must {
    tau := ncast(T, 2)/(conjugate(normal)*normal);

    res : DenseMatrix(T, N, N) = ---;
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
        #run assert(row_to < M.N);
        #run assert(column_to < M.M);
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
        assert(i >= 0 && i < M);
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
        for i: 0..A.N {
            for j: 0..A.M {
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

    res : MatrixMultiplicationHelper(A,B,A.N, B.M).T = ---;
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


#scope_file matrix_vector_multiplication_helper :: ($a: Type, $b: Type, $N: int) -> Type {
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

    res : MatrixVectorMultiplicationHelper(A,B,A.N).T = ---;
    for i: 0..A.N-1 {
        for j: 0..A.M-1 {
            set(*res,i, get(res,i) + get(a,i,j) * get(b,j));
        }
    }

    return res;
}

operator * :: (a: $A/VectorType, b: $B/MatrixType) -> MatrixVectorMultiplicationHelper(A,B,B.M).T {
    #if CHECKS {
        #run assert(A.N == B.N, "Dimensions mismatch: a.dim = %, b.rows = %", A.N, B.N);
    }

    res : MatrixVectorMultiplicationHelper(A,B,B.M).T = ---;
    for j: 0..B.N-1 {
        for i: 0..A.N-1 {
            set(*res,j, get(res,j) + get(a,i) * get(b,i,j));
        }
    }

    return res;
}


operator * :: (a: $A/MatrixType, b: Scalar($B)) -> MatrixCastHelper(A, UpCast(A.data_type, B)).T #symmetric {
    res : MatrixCastHelper(A, UpCast(A.data_type, B)).T = ---;

    for i: 0..A.N {
        for j: 0..A.M {
            set(*res, i, j, get(a, i, j) * b.value);
        }
    }

    return res;
}


mul :: (a: *$A/MatrixType, b: Scalar($B)) {
    #if CHECKS {
        #run assert(is_minor_number(A.data_type, B));
    }

    for i: 0..A.N {
        for j: 0..A.M {
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

operator / :: inline (a: $A/MatrixType, b: Scalar($B)) -> MatrixCastHelper(A, UpCast(A.data_type, B).T).T { 
    res : MatrixCastHelper(A, UpCast(A.data_type, B).T).T = ---;

    for i: 0..A.N {
        for j: 0..A.M {
            set(*res, i, j, get(a, i, j)/b.value);
        }
    }

    return res;
}



div :: (a: *$A/MatrixType, b: Scalar($B)) {
    #if CHECKS {
        #assert(is_minor_number(A.data_type, B));
    }

    for i: 0..A.N {
        for j: 0..A.M {
            set(a, i, j, get(<<a, i, j) / b.value);
        }
    }
}




operator + :: (a: $A/MatrixType, b: $B/MatrixType) -> MatrixHelper(A,B,N,M).T {
    #if CHECKS {
        #run assert(A.N == B.N && A.M == B.M, "matrix dimensions do not match A(% x %), B(% x %)", A.N, A.M, B.N, B.M);
    }
    
    res : MatrixHelper(A,B,N,M).T = ---;

    for i: 0..A.N {
        for j: 0..A.M {
            set(*res, i, j, get(a, i, j) + get(b, i, j));
        }
    }
    return res;
}





add :: (a: *Matrix($A, $N, $M), b:Matrix($B, N, M)) {
    #if CHECKS {
        assert(is_minor_number(A, B));
    }
    for 0..a.data.count-1 {
        a.data[it] += b.data[it];
    }
}








operator - :: (a: Matrix($A, $N, $M), b: Matrix($B, N, M)) -> Matrix(UpCast(A,B).T, N, M) {
    res : Matrix(UpCast(A,B).T, N, M) = ---;
    res.flags = .NONE;

    for 0..a.data.count-1 {
        res.data[it] = a.data[it] - b.data[it];
    }
    return res;
}



sub :: (a: *Matrix($A, $N, $M), b:Matrix($B, N, M)) {
    #if CHECKS {
        #assert(is_minor_number(A, B));
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
        #assert(is_minor_number(C, A) && is_minor_number(C, B));
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
        assert(i < N && j < N);
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
tensor :: (a: Matrix($A, $aN, $aM), b: Matrix($B, $bN, $bM)) -> Matrix(UpCast(A,B).T, aN*bN, aM*bM) {
    res : Matrix(UpCast(A,B).T, aN * bN, aM * bM) = ---;
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
        print("2 b = \n%\n", scl(2.0) * b);
        assert(a == scl(2.0) * b);

        a = a + b;
        assert(a == scl(3.0) * b);
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