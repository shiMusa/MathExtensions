





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
    flags: MATRIX_FLAGS,
    get: (m: matrix_type, int, int) -> data_type,
    set: (m: *matrix_type, int, int, data_type) -> ()
){}

get :: inline (m: $M/MatrixType, i: int, j: int) -> M.data_type {
    return m.get(m, i, j);
}

set :: inline (m: *$M/MatrixType, i: int, j: int, value: M.data_type) {
    m.set(m, i, j, value);
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
        .NONE,
        dense_matrix_get,
        dense_matrix_set
    ); 
    data: [_N*_M]_T;
}

#scope_file
dense_matrix_get :: inline (m: DenseMatrix($T, $N, $M), i: int, j: int) -> T {
    return m.data[i*N + j];
}
dense_matrix_set :: inline (m: *DenseMatrix($T, $N, $M), i: int, j: int, value: T) {
    m.data[i*N + j] = value;
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
    return DenseMatrix(b, a.N, a.M);
}


#scope_export
MatrixSubmatrixHelper :: struct(A: Type, N: int, M: int) {
    T :: #run matrix_submatrix_helper(A, N, M);
}
#scope_file
matrix_submatrix_helper :: ($a: Type, N: int, M: int) -> Type {
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