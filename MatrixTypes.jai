
// ######## ##          ###     ######    ######
// ##       ##         ## ##   ##    ##  ##    ##
// ##       ##        ##   ##  ##        ##
// ######   ##       ##     ## ##   ####  ######
// ##       ##       ######### ##    ##        ##
// ##       ##       ##     ## ##    ##  ##    ##
// ##       ######## ##     ##  ######    ######


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


MATRIX_FLAGS :: enum_flags u32 {
    VIEW  :: 0;
    DENSE :: (1 << 0);
    HEAP  :: (1 << 1);
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
    rows: (m: matrix_type) -> int,
    cols: (m: matrix_type) -> int,
    flags: MATRIX_FLAGS,

    // * (row,column) access
    get: (m: matrix_type, int, int) -> data_type,
    set: (m: *matrix_type, int, int, data_type) -> (),
    
    // * linear access to all elements of the matrix
    get_lin: (m: matrix_type, int) -> data_type,
    set_lin: (m: *matrix_type, int, data_type) -> (),

    // * raw data access
    n_raw : (m: matrix_type) -> int,
    get_raw : (m: matrix_type, int) -> data_type,
    set_raw: (m: *matrix_type, int, data_type) -> (),
    from_raw_index: (m: matrix_type, int) -> int
){}


n :: inline (m: $M/MatrixType) -> int {
    return M.rows(m) * M.cols(m);
}
n :: inline (m: *$M/MatrixType) -> int {
    return M.rows(<<m) * M.cols(<<m);
}

rows :: inline (m: $M/MatrixType) -> int {
    return M.rows(m);
}
rows :: inline (m: *$M/MatrixType) -> int {
    return M.rows(<<m);
}
cols :: inline (m: $M/MatrixType) -> int {
    return M.cols(m);
}
cols :: inline (m: *$M/MatrixType) -> int {
    return M.cols(<<m);
}


get :: inline (m: $M/MatrixType, i: int, j: int) -> M.data_type {
    return M.get(m, i, j);
}
get :: inline (m: *$M/MatrixType, i: int, j: int) -> M.data_type {
    return M.get(<<m, i, j);
}
set :: inline (m: *$M/MatrixType, i: int, j: int, value: M.data_type) {
    M.set(m, i, j, value);
}

get_lin :: inline (m: $M/MatrixType, i: int) -> M.data_type {
    return M.get_lin(m, i);
}
get_lin :: inline (m: *$M/MatrixType, i: int) -> M.data_type {
    return M.get_lin(<<m, i);
}
set_lin :: inline (m: *$M/MatrixType, i: int, value: M.data_type) {
    M.set_lin(m, i, value);
}
from_lin_index :: inline (m: $M/MatrixType, i: int) -> int, int {
    c := cols(m);
    return (i - i%c)/c, i%c;
}
to_lin_index :: inline (m: $M/MatrixType, i: int, j: int) -> int {
    return i * cols(m) + j;
}

n_raw :: inline (m: $M/MatrixType) -> int {
    return M.n_raw(m);
}
get_raw :: inline (m: $M/MatrixType, i: int) -> M.data_type {
    return M.get_raw(m, i);
}
set_raw :: inline (m: *$M/MatrixType, i: int, value: M.data_type) {
    M.set_raw(m, i, value);
}
from_raw_index :: inline (m: $M/MatrixType, i: int) -> int {
    return M.from_raw_index(m, i);
}



// ########  #######  ########          ######## ##     ## ########  ##    ##
// ##       ##     ## ##     ##         ##        ##   ##  ##     ## ###   ##
// ##       ##     ## ##     ##         ##         ## ##   ##     ## ####  ##
// ######   ##     ## ########  ####### ######      ###    ########  ## ## ##
// ##       ##     ## ##   ##           ##         ## ##   ##        ##  ####
// ##       ##     ## ##    ##          ##        ##   ##  ##        ##   ###
// ##        #######  ##     ##         ######## ##     ## ##        ##    ##


// * standard for expansion of all elements of the matrix
for_expansion :: (m: *$M/MatrixType, body: Code, flags: For_Flags) #expand {
    i := 0;
    while i < n(m) {
        `it_index := i;
        `it := get_lin(<<m, i);
        `it_row, `it_column := from_lin_index(<<m, i);
        #insert body;
        i += 1;
    }
}



raw :: inline (m: *$M/MatrixType) -> MatrixRawForExpansion(M) {
    res : MatrixRawForExpansion(M) = ---;
    res.matrix = m;
    return res;
}
raw :: inline (m: $M/MatrixType) -> MatrixRawForExpansion(M) {
    res : MatrixRawForExpansion(M) = ---;
    res.matrix = *m;
    return res;
}

#scope_file
MatrixRawForExpansion :: struct(T: Type) {
    matrix: *T;
}
#scope_export

// * for expansion for only the really existing data
for_expansion :: (me: *MatrixRawForExpansion($M), body: Code, flags: For_Flags) #expand {
    m := <<me.matrix;
    i := 0;
    while i < n_raw(m) {
        `it_raw_index := i;
        `it_index := from_raw_index(m, i);
        `it := get_raw(m, i);
        `it_row, `it_column := from_lin_index(m, from_raw_index(m, i));
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


// ##     ## #### ######## ##      ##
// ##     ##  ##  ##       ##  ##  ##
// ##     ##  ##  ##       ##  ##  ##
// ##     ##  ##  ######   ##  ##  ##
//  ##   ##   ##  ##       ##  ##  ##
//   ## ##    ##  ##       ##  ##  ##
//    ###    #### ########  ###  ###


MatrixView :: struct(V: Type, N: int, M: int) {
    using _m: MatrixType(
        #this,
        V.data_type,
        matrix_view_rows,
        matrix_view_cols,
        .VIEW & V.flags,

        matrix_view_get,
        matrix_view_set,

        matrix_view_get_lin,
        matrix_view_set_lin,

        matrix_view_n_raw,
        matrix_view_get_lin, // ? ther is no difference between raw and lin for DENSE matrices
        matrix_view_set_lin, // ? ther is no difference between raw and lin for DENSE matrices
        matrix_view_from_raw_index
    );

    start_row: int = 0;
    start_col: int = 0;
    stride_row: int = 1;
    stride_col: int = 1;
    mat: *V;
}


matview :: (mat: *$M/MatrixType, $num_rows: int, $num_cols: int, start_row: int = 0, start_col: int = 0, stride_row: int = 1, stride_col: int = 1) -> MatrixView(M, num_rows, num_cols) {
    #if CHECKS {
        assert(start_row + (num_rows-1)*stride_row <=  rows(mat));
        assert(start_col + (num_cols-1)*stride_col <=  cols(mat));
    }
    res : MatrixView(M, num_rows, num_cols);
    res.start_row = start_row;
    res.start_col = start_col;
    res.stride_row = stride_row;
    res.stride_col = stride_col;
    res.mat = mat;
    return res;
}


#scope_file
matrix_view_rows :: inline (m: MatrixView($V, $N, $M)) -> int {
    return N;
}
matrix_view_cols :: inline (m: MatrixView($V, $N, $M)) -> int {
    return M;
}

matrix_view_get :: inline (m: MatrixView, i: int, j: int) -> m.data_type {
    return get(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col);
}
matrix_view_set :: inline (m: *MatrixView, i: int, j: int, val: m.data_type) {
    set(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col, val);
}

matrix_view_get_lin :: inline (m: MatrixView, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return matrix_view_get(m, r, c);
}
matrix_view_set_lin :: inline (m: *MatrixView, i: int, val: m.data_type) {
    r, c := from_lin_index(m, i);
    matrix_view_set(m, r, c, val);
}

matrix_view_n_raw :: inline (m: MatrixView($T, $N, $M)) -> int {
    return N*M;
}
matrix_view_from_raw_index :: inline(m: MatrixView, i: int) -> int {
    return i;
}
#scope_export


// ##     ## ########    ###    ########     ##     ## #### ######## ##      ##
// ##     ## ##         ## ##   ##     ##    ##     ##  ##  ##       ##  ##  ##
// ##     ## ##        ##   ##  ##     ##    ##     ##  ##  ##       ##  ##  ##
// ######### ######   ##     ## ########     ##     ##  ##  ######   ##  ##  ##
// ##     ## ##       ######### ##            ##   ##   ##  ##       ##  ##  ##
// ##     ## ##       ##     ## ##             ## ##    ##  ##       ##  ##  ##
// ##     ## ######## ##     ## ##              ###    #### ########  ###  ###

HeapMatrixView :: struct(V: Type) {
    using _m: MatrixType(
        #this,
        V.data_type,
        heap_matrix_view_rows,
        heap_matrix_view_cols,
        .VIEW & V.flags,

        heap_matrix_view_get,
        heap_matrix_view_set,

        heap_matrix_view_get_lin,
        heap_matrix_view_set_lin,

        heap_matrix_view_n_raw,
        heap_matrix_view_get_lin, // ? ther is no difference between raw and lin for DENSE matrices
        heap_matrix_view_set_lin, // ? ther is no difference between raw and lin for DENSE matrices
        heap_matrix_view_from_raw_index
    );

    rows: int;
    cols : int;
    start_row: int = 0;
    start_col: int = 0;
    stride_row: int = 1;
    stride_col: int = 1;
    mat: *V;
}


hmatview :: (mat: *$M/MatrixType, num_rows: int, num_cols: int, start_row: int = 0, start_col: int = 0, stride_row: int = 1, stride_col: int = 1) -> HeapMatrixView(M) {
    #if CHECKS {
        assert(start_row + (num_rows-1)*stride_row <=  rows(mat));
        assert(start_col + (num_cols-1)*stride_col <=  cols(mat));
    }
    res : HeapMatrixView(M);
    res.rows = num_rows;
    res.cols = num_cols;
    res.start_row = start_row;
    res.start_col = start_col;
    res.stride_row = stride_row;
    res.stride_col = stride_col;
    res.mat = mat;
    return res;
}


#scope_file
heap_matrix_view_rows :: inline (m: HeapMatrixView) -> int {
    return m.rows;
}
heap_matrix_view_cols :: inline (m: HeapMatrixView) -> int {
    return m.cols;
}

heap_matrix_view_get :: inline (m: HeapMatrixView, i: int, j: int) -> m.data_type {
    return get(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col);
}
heap_matrix_view_set :: inline (m: *HeapMatrixView, i: int, j: int, val: m.data_type) {
    set(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col, val);
}

heap_matrix_view_get_lin :: inline (m: HeapMatrixView, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return heap_matrix_view_get(m, r, c);
}
heap_matrix_view_set_lin :: inline (m: *HeapMatrixView, i: int, val: m.data_type) {
    r, c := from_lin_index(m, i);
    heap_matrix_view_set(m, r, c, val);
}

heap_matrix_view_n_raw :: inline (m: HeapMatrixView) -> int {
    return m.rows * m.cols;
}
heap_matrix_view_from_raw_index :: inline(m: HeapMatrixView, i: int) -> int {
    return i;
}
#scope_export



// ########  ######## ##    ##  ######  ########
// ##     ## ##       ###   ## ##    ## ##
// ##     ## ##       ####  ## ##       ##
// ##     ## ######   ## ## ##  ######  ######
// ##     ## ##       ##  ####       ## ##
// ##     ## ##       ##   ### ##    ## ##
// ########  ######## ##    ##  ######  ########




// The only data that actually exists in memory during runtime is "data", the ONELY field of the struct.
DenseMatrix :: struct(_T: Type, _N: int, _M: int) {
    using _m: MatrixType(
        #this,
        _T,
        dense_matrix_rows,
        dense_matrix_cols,
        .DENSE,

        dense_matrix_get,
        dense_matrix_set,
        
        dense_matrix_get_lin,
        dense_matrix_set_lin,

        dense_matrix_n_raw,
        dense_matrix_get_lin, // ? ther is no difference between raw and lin for DENSE matrices
        dense_matrix_set_lin, // ? ther is no difference between raw and lin for DENSE matrices
        dense_matrix_from_raw_index
    ); 
    data: [_N*_M]_T;
}

dmat :: ($T: Type, $N: int, $M: int, values: ..$U) -> DenseMatrix(T,N,M) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        assert(values.count <= N*M);
    }
    res : DenseMatrix(T, N, M) = ---;
    for res {
        if it_index < values.count {
            set_lin(*res, it_index, ncast(T, values[it_index]));
        } else {
            set_lin(*res, it_index, zero(T));
        }
    }
    return res;
}

dmat :: inline ($N: int, $M: int, values: ..$T) -> DenseMatrix(T,N,M) {
    #if CHECKS {
        assert(values.count <= N*M);
    }
    res : DenseMatrix(T, N, M) = ---;
    for res {
        if it_index < values.count {
            set_lin(*res, it_index, values[it_index]);
        } else {
            set_lin(*res, it_index, zero(T));
        }
    }
    return res;
}

#scope_file
dense_matrix_rows :: inline (m: DenseMatrix($T, $N, $M)) -> int {
    return N;
}
dense_matrix_cols :: inline (m: DenseMatrix($T, $N, $M)) -> int {
    return M;
}

dense_matrix_get :: inline (m: DenseMatrix($T, $N, $M), i: int, j: int) -> T {
    return m.data[i*M + j];
}
dense_matrix_set :: inline (m: *DenseMatrix($T, $N, $M), i: int, j: int, value: T) {
    m.data[i*M + j] = value;
}


dense_matrix_get_lin :: inline (m: DenseMatrix($T, $N, $M), i: int) -> T {
    return m.data[i];
}
dense_matrix_set_lin :: inline (m: *DenseMatrix($T, $N, $M), i: int, value: T) {
    m.data[i] = value;
}

dense_matrix_n_raw :: inline (m: DenseMatrix($T, $N, $M)) -> int {
    return N*M;
}
dense_matrix_from_raw_index :: inline (m: DenseMatrix($T, $N, $M), i: int) -> int {
    return i;
}
#scope_export


// ########  ######## ##    ##  ######  ########    ##     ## ########    ###    ########
// ##     ## ##       ###   ## ##    ## ##          ##     ## ##         ## ##   ##     ##
// ##     ## ##       ####  ## ##       ##          ##     ## ##        ##   ##  ##     ##
// ##     ## ######   ## ## ##  ######  ######      ######### ######   ##     ## ########
// ##     ## ##       ##  ####       ## ##          ##     ## ##       ######### ##
// ##     ## ##       ##   ### ##    ## ##          ##     ## ##       ##     ## ##
// ########  ######## ##    ##  ######  ########    ##     ## ######## ##     ## ##


DenseHeapMatrix :: struct(_T: Type) {
    using _m: MatrixType(
        #this,
        _T,
        dense_heap_matrix_rows,
        dense_heap_matrix_cols,
        .DENSE & .HEAP,

        dense_heap_matrix_get,
        dense_heap_matrix_set,
        
        dense_heap_matrix_get_lin,
        dense_heap_matrix_set_lin,

        dense_heap_matrix_n_raw,
        dense_heap_matrix_get_lin, // ? ther is no difference between raw and lin for DENSE matrices
        dense_heap_matrix_set_lin, // ? ther is no difference between raw and lin for DENSE matrices
        dense_heap_matrix_from_raw_index
    ); 
    rows: int;
    cols: int;
    data : []_T; // * heap-allocated data!
}

free :: inline (m: $M/DenseHeapMatrix) {
    array_free(m.data);
}

dhmat :: ($T: Type, N: int, M: int, init_zeroes: bool = true) -> DenseHeapMatrix(T) {
    res : DenseHeapMatrix(T) = ---;
    res.rows = N;
    res.cols = M;
    res.data = NewArray(N*M, T, init_zeroes);
    return res;
}

// TODO simplify this to only take `values` of type `T`?
dhmat :: ($T: Type, N: int, M: int, values: ..$U) -> DenseHeapMatrix(T) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        assert(values.count <= N*M);
    }
    res : DenseHeapMatrix(T) = ---;
    res.rows = N;
    res.cols = M;
    res.data = NewArray(N*M, T);

    for values {
        set_lin(*res, it_index, ncast(T,it));
    }
    return res;
}

dhmat :: inline (N: int, M: int, values: ..$T) -> DenseHeapMatrix(T) {
    #if CHECKS {
        assert(values.count <= N*M);
    }
    res : DenseHeapMatrix(T);
    res.rows = N;
    res.cols = M;
    res.data = NewArray(N*M, T);

    for values {
        set_lin(*res, it_index, it);
    }

    return res;
}

#scope_file
dense_heap_matrix_rows :: inline (m: DenseHeapMatrix($T)) -> int {
    return m.rows;
}
dense_heap_matrix_cols :: inline (m: DenseHeapMatrix($T)) -> int {
    return m.cols;
}

dense_heap_matrix_get :: inline (m: DenseHeapMatrix($T), i: int, j: int) -> T {
    return m.data[i*cols(m) + j];
}
dense_heap_matrix_set :: inline (m: *DenseHeapMatrix($T), i: int, j: int, value: T) {
    m.data[i*cols(m) + j] = value;
}


dense_heap_matrix_get_lin :: inline (m: DenseHeapMatrix($T), i: int) -> T {
    return m.data[i];
}
dense_heap_matrix_set_lin :: inline (m: *DenseHeapMatrix($T), i: int, value: T) {
    m.data[i] = value;
}

dense_heap_matrix_n_raw :: inline (m: DenseHeapMatrix($T)) -> int {
    return rows(m)*cols(m);
}
dense_heap_matrix_from_raw_index :: inline (m: DenseHeapMatrix($T), i: int) -> int {
    return i;
}
#scope_export