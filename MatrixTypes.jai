





// ######## ##          ###     ######    ######
// ##       ##         ## ##   ##    ##  ##    ##
// ##       ##        ##   ##  ##        ##
// ######   ##       ##     ## ##   ####  ######
// ##       ##       ######### ##    ##        ##
// ##       ##       ##     ## ##    ##  ##    ##
// ##       ######## ##     ##  ######    ######

MATRIX_FLAGS :: enum_flags u32 {
    DENSE :: 0;
    // band
    // symmetric
    // hermitian
    // triangular
    // sparse
    // etc.
}



// ######## ########     ###    #### ########
//    ##    ##     ##   ## ##    ##     ##
//    ##    ##     ##  ##   ##   ##     ##
//    ##    ########  ##     ##  ##     ##
//    ##    ##   ##   #########  ##     ##
//    ##    ##    ##  ##     ##  ##     ##
//    ##    ##     ## ##     ## ####    ##

MatrixType :: struct(
    matrix_type: Type, 
    data_type: Type,
    N: int,
    M: int,
    N_lin: int,
    flags: MATRIX_FLAGS,
    get: (m: matrix_type, int, int) -> data_type,
    set: (m: *matrix_type, int, int, data_type) -> (),
    get_lin: (m: matrix_type, int) -> data_type,
    set_lin: (m: *matrix_type, int, data_type) -> (),
    from_lin_index: (m: matrix_type, int) -> (int, int),
    to_lin_index: (m: matrix_type, int, int) -> int
){}
from_lin_index :: inline (m: $M/MatrixType, i: int) -> int, int {
    r, c := m.from_lin_index(m, i);
    return r, c;
}
to_lin_index :: inline (m: $M/MatrixType, i: int, j: int) -> int {
    return m.to_lin_index(m, i,j);
}
get :: inline (m: $M/MatrixType, i: int, j: int) -> M.data_type {
    return m.get(m, i, j);
}
set :: inline (m: *$M/MatrixType, i: int, j: int, value: M.data_type) {
    m.set(m, i, j, value);
}
get_lin :: inline (m: $M/MatrixType, i: int) -> M.data_type {
    return m.get_lin(m, i);
}
set_lin :: inline (m: *$M/MatrixType, i: int, value: M.data_type) {
    m.set_lin(m, i, value);
}
for_expansion :: (m: *$M/MatrixType, body: Code, flags: For_Flags) #expand {
    i := 0;
    while i < M.N_lin {
        `it_index := i;
        `it := get_lin(<<m, i);
        `it_row, `it_column := from_lin_index(<<m, i);
        #insert body;
        i += 1;
    }
}




// ##     ##    ###    ########    ######## ##    ## ########  ########  ######
// ###   ###   ## ##      ##          ##     ##  ##  ##     ## ##       ##    ##
// #### ####  ##   ##     ##          ##      ####   ##     ## ##       ##
// ## ### ## ##     ##    ##          ##       ##    ########  ######    ######
// ##     ## #########    ##          ##       ##    ##        ##             ##
// ##     ## ##     ##    ##          ##       ##    ##        ##       ##    ##
// ##     ## ##     ##    ##          ##       ##    ##        ########  ######

// The only data that actually exists in memory during runtime is "data", the ONELY field of the struct.
DenseMatrix :: struct(_T: Type, _N: int, _M: int) {
    using _m: MatrixType(
        DenseMatrix(_T, _N, _M),
        _T,
        _N,
        _M,
        _N * _M,
        .DENSE,
        dense_matrix_get,
        dense_matrix_set,
        dense_matrix_get_lin,
        dense_matrix_set_lin,
        dense_matrix_from_lin_index,
        dense_matrix_to_lin_index
    ); 
    data: [_N*_M]_T;
}

dmat :: ($T: Type, $N: int, $M: int, values: ..$U) -> DenseMatrix(T,N,M) {
    #if CHECKS {
        #run assert(is_minor_number(T, U));
        assert(values.count <= N*M);
    }
    res : DenseMatrix(T, N, M);
    for values {
        set_lin(*res, it_index, ncast(T,it));
    }
    return res;
}

dmat :: inline ($N: int, $M: int, values: ..$T) -> DenseMatrix(T,N,M) {
    #if CHECKS {
        assert(values.count <= N*M);
    }
    res : DenseMatrix(T,N,M);
    // println("dmat:");
    // println("    res %", res);
    // println("    values %", values);
    mat_init(*res, ..values);
    return res;
}

#scope_file
dense_matrix_get :: inline (m: DenseMatrix($T, $N, $M), i: int, j: int) -> T {
    return m.data[i*M + j];
}
dense_matrix_set :: inline (m: *DenseMatrix($T, $N, $M), i: int, j: int, value: T) {
    m.data[i*M + j] = value;
}
dense_matrix_from_lin_index :: inline (m: DenseMatrix($T, $N, $M), i: int) -> int, int {
    return (i- i%M)/M, i%M;
}
dense_matrix_to_lin_index :: inline (m: DenseMatrix($T, $N, $M), i: int, j: int) -> int {
    return i*M+j;
}
dense_matrix_get_lin :: inline (m: DenseMatrix($T, $N, $M), i: int) -> T {
    return m.data[i];
}
dense_matrix_set_lin :: inline (m: *DenseMatrix($T, $N, $M), i: int, value: T) {
    m.data[i] = value;
}
#scope_export




// ##     ## ######## ##       ########  ######## ########
// ##     ## ##       ##       ##     ## ##       ##     ##
// ##     ## ##       ##       ##     ## ##       ##     ##
// ######### ######   ##       ########  ######   ########
// ##     ## ##       ##       ##        ##       ##   ##
// ##     ## ##       ##       ##        ##       ##    ##
// ##     ## ######## ######## ##        ######## ##     ##

MatrixHelper :: struct(A: Type, B: Type, N: int, M: int) {
    T :: #run matrix_helper(A,B,N,M);
}

// TODO: maybe be smart about this at some point?
#scope_file
matrix_helper :: ($a: Type, $b: Type, $N: int, $M: int) -> Type {
    return DenseMatrix(UpCast(a.data_type, b.data_type).T, N, M);
}

#scope_export
MatrixCastHelper :: struct(A: Type, data_type: Type) {
    T :: #run matrix_cast_helper(A,data_type);
}
#scope_file
matrix_cast_helper :: ($a: Type, $b: Type) -> Type {
    // TODO other types!
    return DenseMatrix(UpCast(a.data_type, b).T, a.N, a.M);
}

#scope_export
MatrixMatrixCastHelper :: struct(A: Type, B: Type) {
    T :: #run matrix_matrix_cast_helper(A,B);
}
#scope_file
matrix_matrix_cast_helper :: ($a: Type, $b: Type) -> Type {
    // TODO other types!
    return DenseMatrix(UpCast(a.data_type, b.data_type).T, a.N, a.M);
}


#scope_export
MatrixSubmatrixHelper :: struct(A: Type, N: int, M: int) {
    T :: #run matrix_submatrix_helper(A, N, M);
}
#scope_file
matrix_submatrix_helper :: ($a: Type, $N: int, $M: int) -> Type {
    // TODO other types!
    return DenseMatrix(a.data_type, N, M);
}


#scope_export 
MatrixVectorHelper :: struct(A: Type, N: int) {
    T :: #run matrix_vector_helper(A, N);
}
#scope_file
matrix_vector_helper :: ($a: Type, $N: int) -> Type {
    // TODO other types!
    return DenseVector(a.data_type, N);
}

#scope_export 
MatrixMultiplicationHelper :: struct(A: Type, B: Type, N: int, M: int) {
    T :: #run matrix_multiplication_helper(A,B,N,M);
}
#scope_file
matrix_multiplication_helper :: ($a: Type, $b: Type, $N: int, $M: int) -> Type {
    return DenseMatrix(UpCast(a.data_type, b.data_type).T, N, M);
}