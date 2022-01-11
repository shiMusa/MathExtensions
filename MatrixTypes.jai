
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
    BAND  :: (1 << 2);
    TRIANGULAR :: (1 << 3);
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


/*
* All concrete types M(T) of ArrayType2D(T) must implement the following functions:
* get: (M(T), int, int) -> T
* set: (*M(T), int, int, T)
*/
ArrayType2D :: struct(data_type: Type){};

get :: inline (a: *$A/ArrayType2D, i: int, j: int) -> A.data_type {
    return get(<<a, i, j);
}

// * implementation for simple arrays
Array2D :: struct(N: int, M: int, _data_type: Type) {
    using _a2d: ArrayType2D(_data_type);
    data: [N][M]_data_type;
}
array2d :: (values: [$N][$M]$T) -> Array2D(N,M,T) {
    res : Array2D(N,M,T) = ---;
    res.data = values;
    return res;
}
get :: inline (a: Array2D, i: int, j: int) -> a.data_type {
    return a.data[i][j];
}
set :: inline (a: *Array2D, i: int, j: int, val: a.data_type) {
    a.data[i][j] = val;
}


/*
* All concrete types M(T) of MatrixType(T, flags) must implement the following functions:
* 
* rows: (m: M(T)) -> int
* cols: (m: M(T)) -> int
* 
* get_lin: (m: M(T), int) -> T
* set_lin: (m: *M(T), int, T) -> ()
* 
* n_raw : (m: M(T)) -> int
* get_raw : (m: M(T), int) -> T
* set_raw: (m: *M(T), int, T) -> ()
* from_raw_index: (m: M(T), int) -> int
*/
MatrixType :: struct(
    __data_type: Type,
    flags: MATRIX_FLAGS
){
    using _a2d: ArrayType2D(__data_type);
}


get_lin :: inline (m: *$M/MatrixType, i: int) -> m.data_type {
    return get_lin(<<m, i);
}
get_raw :: inline (m: *$M/MatrixType, i: int) -> m.data_type {
    return get_raw(<<m, i);
}
rows :: inline (m: *$M/MatrixType) -> int {
    return rows(<<m);
}
cols :: inline (m: *$M/MatrixType) -> int {
    return cols(<<m);
}

from_raw_index :: inline (m: *$M/MatrixType, i: int) -> int {
    return from_raw_index(<<m, i);
}

n :: inline (m: $M/MatrixType) -> int {
    return rows(m) * M.cols(m);
}
n :: inline (m: *$M/MatrixType) -> int {
    return rows(<<m) * cols(<<m);
}

from_lin_index :: inline (m: $M/MatrixType, i: int) -> int, int {
    c := cols(m);
    return (i - i%c)/c, i%c;
}
to_lin_index :: inline (m: $M/MatrixType, i: int, j: int) -> int {
    return i * cols(m) + j;
}
from_lin_index :: inline (m: *$M/MatrixType, i: int) -> int, int {
    r, c := inline from_lin_index(<<m, i);
    return r, c;
}
to_lin_index :: inline (m: *$M/MatrixType, i: int, j: int) -> int {
    return inline to_lin_index(<<m, i, j);
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
        #insert (continue={i += 1; continue;}, break={break;}) body;
        i += 1;
    }
}



raw :: inline (m: *$M/MatrixType) -> *MatrixRawForExpansion(M) {
    res := New(MatrixRawForExpansion(M), false);
    res.matrix = m;
    return res;
}

#scope_file
MatrixRawForExpansion :: struct(T: Type) {
    matrix: *T;
}
#scope_export

// * for expansion for only the really existing data
for_expansion :: (me: *MatrixRawForExpansion, body: Code, flags: For_Flags) #expand {
    m := <<(me.matrix);
    i := 0;
    while i < n_raw(m) {
        `it_raw_index := i;
        `it_index := from_raw_index(m, i);
        `it := get_raw(m, i);
        `it_row, `it_column := from_lin_index(m, from_raw_index(m, i));
        #insert (continue={i += 1; continue;}, break={break;}) body;
        i += 1;
    }
    free(me);
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
    using _m: MatrixType(V.data_type, MATRIX_FLAGS.VIEW & V.flags);
    start_row: int = 0;
    start_col: int = 0;
    stride_row: int = 1;
    stride_col: int = 1;
    mat: *V;
}


mat_view :: (mat: *$M/MatrixType, $num_rows: int, $num_cols: int, start_row: int = 0, start_col: int = 0, stride_row: int = 1, stride_col: int = 1) -> MatrixView(M, num_rows, num_cols) {
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

rows :: inline (m: MatrixView) -> int {
    return m.N;
}
cols :: inline (m: MatrixView) -> int {
    return m.M;
}

get :: inline (m: MatrixView, i: int, j: int) -> m.data_type {
    return get(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col);
}
set :: inline (m: *MatrixView, i: int, j: int, val: m.data_type) {
    set(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col, val);
}

get_lin :: inline (m: MatrixView, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_lin :: inline (m: *MatrixView, i: int, val: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, val);
}

get_raw :: inline (m: MatrixView, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_raw :: inline (m: *MatrixView, i: int, val: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, val);
}

n_raw :: inline (m: MatrixView) -> int {
    return m.N*m.M;
}
from_raw_index :: inline(m: MatrixView, i: int) -> int {
    return i;
}


// ##     ## ########    ###    ########     ##     ## #### ######## ##      ##
// ##     ## ##         ## ##   ##     ##    ##     ##  ##  ##       ##  ##  ##
// ##     ## ##        ##   ##  ##     ##    ##     ##  ##  ##       ##  ##  ##
// ######### ######   ##     ## ########     ##     ##  ##  ######   ##  ##  ##
// ##     ## ##       ######### ##            ##   ##   ##  ##       ##  ##  ##
// ##     ## ##       ##     ## ##             ## ##    ##  ##       ##  ##  ##
// ##     ## ######## ##     ## ##              ###    #### ########  ###  ###

MatrixHeapView :: struct(V: Type) {
    using _m: MatrixType(V.data_type, MATRIX_FLAGS.VIEW & V.flags);
    rows: int;
    cols : int;
    start_row: int = 0;
    start_col: int = 0;
    stride_row: int = 1;
    stride_col: int = 1;
    mat: *V;
}


mat_heap_view :: (mat: *$M/MatrixType, num_rows: int, num_cols: int, start_row: int = 0, start_col: int = 0, stride_row: int = 1, stride_col: int = 1) -> MatrixHeapView(M) {
    #if CHECKS {
        assert(start_row + (num_rows-1)*stride_row <=  rows(mat));
        assert(start_col + (num_cols-1)*stride_col <=  cols(mat));
    }
    res : MatrixHeapView(M);
    res.rows = num_rows;
    res.cols = num_cols;
    res.start_row = start_row;
    res.start_col = start_col;
    res.stride_row = stride_row;
    res.stride_col = stride_col;
    res.mat = mat;
    return res;
}


rows :: inline (m: MatrixHeapView) -> int {
    return m.rows;
}
cols :: inline (m: MatrixHeapView) -> int {
    return m.cols;
}

get :: inline (m: MatrixHeapView, i: int, j: int) -> m.data_type {
    return get(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col);
}
set :: inline (m: *MatrixHeapView, i: int, j: int, val: m.data_type) {
    set(m.mat, m.start_row + i * m.stride_row, m.start_col + j * m.stride_col, val);
}

get_lin :: inline (m: MatrixHeapView, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_lin :: inline (m: *MatrixHeapView, i: int, val: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, val);
}

get_raw :: inline (m: MatrixHeapView, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_raw :: inline (m: *MatrixHeapView, i: int, val: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, val);
}
n_raw :: inline (m: MatrixHeapView) -> int {
    return m.rows * m.cols;
}
from_raw_index :: inline(m: MatrixHeapView, i: int) -> int {
    return i;
}



// ########  ######## ##    ##  ######  ########
// ##     ## ##       ###   ## ##    ## ##
// ##     ## ##       ####  ## ##       ##
// ##     ## ######   ## ## ##  ######  ######
// ##     ## ##       ##  ####       ## ##
// ##     ## ##       ##   ### ##    ## ##
// ########  ######## ##    ##  ######  ########




// The only data that actually exists in memory during runtime is "data", the ONELY field of the struct.
DenseMatrix :: struct(_T: Type, N: int, M: int) {
    using _m: MatrixType(_T, MATRIX_FLAGS.DENSE); 
    data: [N*M]_T;
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

rows :: inline (m: DenseMatrix) -> int {
    return m.N;
}
cols :: inline (m: DenseMatrix) -> int {
    return m.M;
}

get :: inline (m: DenseMatrix, i: int, j: int) -> m.data_type {
    return m.data[i*m.M + j];
}
set :: inline (m: *DenseMatrix, i: int, j: int, value: m.data_type) {
    m.data[i*m.M + j] = value;
}

get_lin :: inline (m: DenseMatrix, i: int) -> m.data_type {
    return m.data[i];
}
set_lin :: inline (m: *DenseMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}

get_raw :: inline (m: DenseMatrix, i: int) -> m.data_type {
    return m.data[i];
}
set_raw :: inline (m: *DenseMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}
n_raw :: inline (m: DenseMatrix) -> int {
    return m.N*m.M;
}
from_raw_index :: inline (m: DenseMatrix, i: int) -> int {
    return i;
}


// ########  ######## ##    ##  ######  ########    ##     ## ########    ###    ########
// ##     ## ##       ###   ## ##    ## ##          ##     ## ##         ## ##   ##     ##
// ##     ## ##       ####  ## ##       ##          ##     ## ##        ##   ##  ##     ##
// ##     ## ######   ## ## ##  ######  ######      ######### ######   ##     ## ########
// ##     ## ##       ##  ####       ## ##          ##     ## ##       ######### ##
// ##     ## ##       ##   ### ##    ## ##          ##     ## ##       ##     ## ##
// ########  ######## ##    ##  ######  ########    ##     ## ######## ##     ## ##


DenseHeapMatrix :: struct(_T: Type) {
    using _m: MatrixType(_T, MATRIX_FLAGS.DENSE & .HEAP); 
    rows: int;
    cols: int;
    data : []_T; // * heap-allocated data!
}

free :: inline (m: DenseHeapMatrix) {
    array_free(m.data);
}

dhmat :: ($T: Type, rows: int, cols: int, $init_zeroes: bool = true) -> DenseHeapMatrix(T) {
    res : DenseHeapMatrix(T) = ---;
    res.rows = rows;
    res.cols = cols;
    res.data = NewArray(rows*cols, T, init_zeroes);
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

rows :: inline (m: DenseHeapMatrix) -> int {
    return m.rows;
}
cols :: inline (m: DenseHeapMatrix) -> int {
    return m.cols;
}

get :: inline (m: DenseHeapMatrix, i: int, j: int) -> m.data_type {
    return m.data[i*m.cols + j];
}
set :: inline (m: *DenseHeapMatrix, i: int, j: int, value: m.data_type) {
    m.data[i*m.cols + j] = value;
}


get_lin :: inline (m: DenseHeapMatrix, i: int) -> m.data_type {
    return m.data[i];
}
set_lin :: inline (m: *DenseHeapMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}

get_raw :: inline (m: DenseHeapMatrix, i: int) -> m.data_type {
    return m.data[i];
}
set_raw :: inline (m: *DenseHeapMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}
n_raw :: inline (m: DenseHeapMatrix) -> int {
    return m.rows * m.cols;
}
from_raw_index :: inline (m: DenseHeapMatrix, i: int) -> int {
    return i;
}



// ########     ###    ##    ## ########
// ##     ##   ## ##   ###   ## ##     ##
// ##     ##  ##   ##  ####  ## ##     ##
// ########  ##     ## ## ## ## ##     ##
// ##     ## ######### ##  #### ##     ##
// ##     ## ##     ## ##   ### ##     ##
// ########  ##     ## ##    ## ########

// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? 1.2.1 Band Matrices, page 15 & 17

// ? Band matrices have entries in the diagonal, to the upper bandwidth `q` up from the diagonal, and `p` down from the diagonal.

// ?      |a a a . . . .|        |a a a . . . . . . . . . . . . . . . . .|
// ?      |x x x x . . .|        |x x x x . . . . . . . . . . . . . . . .|
// ?      |. x x x x . .|        |. x x x x . . . . . . . . . . . . . . .|
// ?      |. . x x x x .|        |. . x x x x . . . . . . . . . . . . . .|
// ?      |. . . x x x x|        |. . . x x x x . . . . . . . . . . . . .|
// ?      |. . . . b b b| 
// ?      |. . . . . b b|   
// ?      |. . . . . . b|     
// ?      |. . . . . . .|       
// ? 

// ? P is the lower bandwidth
// ? Q is the upper bandwidth
BandMatrix :: struct(_T: Type, N: int, M: int, P: int, Q: int) {
    using _m: MatrixType(_T, MATRIX_FLAGS.DENSE | MATRIX_FLAGS.BAND); 
    ROWS : int : P+Q+1; // ? the BandMatrix is of size ROWS × M
    data: [ROWS*M]_T;
}

bmat :: ($T: Type, $N: int, $M: int, $P: int, $Q: int, values: ..$U) -> BandMatrix(T,N,M,P,Q) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        #run assert(N>=0 && M>=0 && P>=0 && Q>=0);
        assert(values.count <= (P+Q+1)*M); // ! this is not the correct limit. See notes "Nov. 1st/2nd, 2021" (handwritten) for correct numbers
        // TODO update above forumla with correct value!
    }
    res : BandMatrix(T, N, M, P, Q) = ---;
    for res {
        if it_index < values.count {
            set_lin(*res, it_index, ncast(T, values[it_index]));
        } else {
            set_lin(*res, it_index, zero(T));
        }
    }
    return res;
}

bmat :: inline ($N: int, $M: int, $P: int, $Q: int, values: ..$T) -> BandMatrix(T,N,M,P,Q) {
    #if CHECKS {
        assert(values.count <= N*M);
        assert(N>=0 && M>=0 && P>=0 && Q>=0);
        assert(values.count <= (P+Q+1)*M); // ! this is not the correct limit. See notes "Nov. 1st/2nd, 2021" (handwritten) for correct numbers
        // TODO update above forumla with correct value!
    }
    res : BandMatrix(T, N, M, P, Q) = ---;
    for res {
        if it_index < values.count {
            set_lin(*res, it_index, values[it_index]);
        } else {
            set_lin(*res, it_index, zero(T));
        }
    }
    return res;
}

// ? these are the numbers for rows and columns for the user.
rows :: inline (m: BandMatrix) -> int {
    return m.N;
}
cols :: inline (m: BandMatrix) -> int {
    return m.M;
}

get :: inline (m: BandMatrix, i: int, j: int) -> m.data_type {
    if j < i-m.P || j > i+m.Q then return zero(m.data_type);
    return get_internal(m, i-j+m.Q, j);
}
set :: inline (m: *BandMatrix, i: int, j: int, value: m.data_type) {
    #if CHECKS {
        assert(j >= i-m.P && j <= i+m.Q);
    }
    if j < i-m.P || j > i+m.Q then return;
    set_internal(m, i-j+m.Q, j, value);
}


// * specific to BandMatrix
get_internal :: inline (m: BandMatrix, i: int, j: int) -> m.data_type {
    return m.data[i*m.M + j];
}
set_internal :: inline (m: *BandMatrix, i: int, j: int, value: m.data_type) {
    m.data[i*m.M + j] = value;
}


// TODO there should be a better way to do this
get_lin :: inline (m: BandMatrix, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_lin :: inline (m: *BandMatrix, i: int, value: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, value);
}

get_raw :: inline (m: BandMatrix, i: int) -> m.data_type {
    // ? this is necessary since some values in the data are nonsense!
    r_int := i / m.M;
    c_int := i % m.M;
    r := r_int + c_int - m.Q;
    c := c_int;
    if c < r-m.P || c > r+m.Q then return zero(m.data_type);
    return m.data[i];
}
set_raw :: inline (m: *BandMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}
n_raw :: inline (m: BandMatrix) -> int {
    return m.ROWS * m.M; // ! this is not the correct number!
}
from_raw_index :: inline (m: BandMatrix, i: int) -> int {
    r_int := i / m.M;
    c_int := i % m.M;
    r := r_int + c_int - m.Q;
    c := c_int;
    return to_lin_index(m, r, c);
}


test_BandMatrix :: () {
    println_push("BandMatrix", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        bm : BandMatrix(int, 15, 17, 1, 2) = ---;
        for raw(*bm) {
            set_raw(*bm, it_raw_index, it_raw_index);
        }
        println(bm);
        println(pstr(bm));
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        bm : BandMatrix(int, 8, 5, 2, 0);
        for raw(*bm) {
            set_raw(*bm, it_raw_index, it_raw_index);
        }
        set(*bm, 3, 3, -42);
        set(*bm, 2, 1, 313);
        println(pstr(bm));
    }
}




// ########     ###    ##    ## ########     ##     ## ########    ###    ########
// ##     ##   ## ##   ###   ## ##     ##    ##     ## ##         ## ##   ##     ##
// ##     ##  ##   ##  ####  ## ##     ##    ##     ## ##        ##   ##  ##     ##
// ########  ##     ## ## ## ## ##     ##    ######### ######   ##     ## ########
// ##     ## ######### ##  #### ##     ##    ##     ## ##       ######### ##
// ##     ## ##     ## ##   ### ##     ##    ##     ## ##       ##     ## ##
// ########  ##     ## ##    ## ########     ##     ## ######## ##     ## ##

// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? 1.2.1 Band Matrices, page 15 & 17

// ? P is the lower bandwidth
// ? Q is the upper bandwidth
BandHeapMatrix :: struct(_T: Type) {
    using _m: MatrixType(_T, MATRIX_FLAGS.DENSE | MATRIX_FLAGS.BAND | MATRIX_FLAGS.HEAP); 
    rows: int;
    cols: int;
    irows: int;
    lower_bandwidth: int;
    upper_bandwidth: int;
    data: []_T;
}

free :: inline (m: BandHeapMatrix) {
    array_free(m.data);
}

bhmat :: ($T: Type, N: int, M: int, P: int, Q: int, $init_zeroes: bool = true) -> BandHeapMatrix(T) {
    res : BandHeapMatrix(T) = ---;
    res.rows = N;
    res.cols = M;
    res.lower_bandwidth = P;
    res.upper_bandwidth = Q;
    res.irows = P+Q+1;
    res.data = NewArray(res.irows * M, T, init_zeroes);
    return res;
}


bhmat :: ($T: Type, N: int, M: int, P: int, Q: int, values: ..$U) -> BandHeapMatrix(T) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        assert(N>=0 && M>=0 && P>=0 && Q>=0);
        assert(values.count <= (P+Q+1)*M); // ! this is not the correct limit. See notes "Nov. 1st/2nd, 2021" (handwritten) for correct numbers
        // TODO update above forumla with correct value!
    }
    res := dhmat(T, N, M, P, Q);
    for i: 0..values.count-1 {
        set_lin(*res, i, ncast(T, values[i]));
    }
    return res;
}

bmat :: inline (N: int, M: int, P: int, Q: int, values: ..$T) -> BandHeapMatrix(T) {
    #if CHECKS {
        assert(values.count <= N*M);
        assert(N>=0 && M>=0 && P>=0 && Q>=0);
        assert(values.count <= (P+Q+1)*M); // ! this is not the correct limit. See notes "Nov. 1st/2nd, 2021" (handwritten) for correct numbers
        // TODO update above forumla with correct value!
    }
    res := bhmat(T, N, M, P, Q);
    for i: 0..values.count-1 {
        set_lin(*res, i, values[i]);
    }
    return res;
}

// ? these are the numbers for rows and columns for the user.
rows :: inline (m: BandHeapMatrix) -> int {
    return m.rows;
}
cols :: inline (m: BandHeapMatrix) -> int {
    return m.cols;
}

get :: inline (m: BandHeapMatrix, i: int, j: int) -> m.data_type {
    if j < i-m.lower_bandwidth || j > i+m.upper_bandwidth then return zero(m.data_type);
    return get_internal(m, i-j+m.upper_bandwidth, j);
}
set :: inline (m: *BandHeapMatrix, i: int, j: int, value: m.data_type) {
    P := m.lower_bandwidth;
    Q := m.upper_bandwidth;
    #if CHECKS {
        assert(j >= i-P && j <= i+Q);
    }
    if j < i-P || j > i+Q then return;
    set_internal(m, i-j+Q, j, value);
}


// * specific to BandHeapMatrix
get_internal :: inline (m: BandHeapMatrix, i: int, j: int) -> m.data_type {
    return m.data[i*m.cols + j];
}
set_internal :: inline (m: *BandHeapMatrix, i: int, j: int, value: m.data_type) {
    m.data[i*m.cols + j] = value;
}


// TODO there should be a better way to do this
get_lin :: inline (m: BandHeapMatrix, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_lin :: inline (m: *BandHeapMatrix, i: int, value: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, value);
}

get_raw :: inline (m: BandHeapMatrix, i: int) -> m.data_type {
    // TODO the whole indexing math can be skipped if you know, which values in the BandHeapMatrix are garbage.
    r_int := i / m.cols;
    c_int := i % m.cols;
    P := m.lower_bandwidth;
    Q := m.upper_bandwidth;
    r := r_int + c_int - Q;
    c := c_int;
    if c < r-P || c > r+Q then return zero(m.data_type);
    return m.data[i];
}
set_raw :: inline (m: *BandHeapMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}
n_raw :: inline (m: BandHeapMatrix) -> int {
    return m.irows * m.cols; // ! this is not the correct number!
}
from_raw_index :: inline (m: BandHeapMatrix, i: int) -> int {
    r_int := i / m.cols;
    c_int := i % m.cols;
    r := r_int + c_int - m.upper_bandwidth;
    c := c_int;
    return to_lin_index(m, r, c);
}


test_BandHeapMatrix :: () {
    println_push("BandHeapMatrix", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        bm := bhmat(int, 15, 17, 1, 2, false);
        defer free(bm);
        for raw(*bm) {
            set_raw(*bm, it_raw_index, it_raw_index);
        }
        println(bm);
        println(pstr(bm));
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        bm := bhmat(int, 8, 5, 2, 0, false);
        defer free(bm);
        for raw(*bm) {
            set_raw(*bm, it_raw_index, it_raw_index);
        }
        set(*bm, 3, 3, -42);
        set(*bm, 2, 1, 313);
        println(pstr(bm));
    }
}


// ######## ########  ####    ###    ##    ##  ######   ##     ## ##          ###    ########
//    ##    ##     ##  ##    ## ##   ###   ## ##    ##  ##     ## ##         ## ##   ##     ##
//    ##    ##     ##  ##   ##   ##  ####  ## ##        ##     ## ##        ##   ##  ##     ##
//    ##    ########   ##  ##     ## ## ## ## ##   #### ##     ## ##       ##     ## ########
//    ##    ##   ##    ##  ######### ##  #### ##    ##  ##     ## ##       ######### ##   ##
//    ##    ##    ##   ##  ##     ## ##   ### ##    ##  ##     ## ##       ##     ## ##    ##
//    ##    ##     ## #### ##     ## ##    ##  ######    #######  ######## ##     ## ##     ##

#scope_file
triangluar_data_size :: (n: int, m: int) -> int {
    if n >= m {
        return ((m+1)*m)/2;
    }
    k := m-n;
    return ((m+1)*m)/2 - ((k+1)*k)/2;
}
triangular_first_index_of_row :: inline (cols: int, i: int) -> int {
    k := cols - i;
    return ((cols+1)*cols)/2 - ((k+1)*k)/2;
}
#scope_export

TriangularMatrix :: struct(_T: Type, N: int, M: int) {
    using _m: MatrixType(_T, MATRIX_FLAGS.DENSE | MATRIX_FLAGS.TRIANGULAR);
    transposed: bool = false;
    data: [#run triangluar_data_size(N,M)]_T;
}



tmat :: ($T: Type, $N: int, $M: int, values: ..$U) -> TriangularMatrix(T,N,M) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        #run assert(N >= 0 && M >= 0);
        assert(values.count <= triangluar_data_size(N,M));
    }
    res : TriangularMatrix(T, N, M);
    for values {
        res.data[it_index] = ncast(T, it);
    }
    return res;
}

tmat :: ($N: int, $M: int, values: ..$T) -> TriangularMatrix(T,N,M) {
    #if CHECKS {
        #run assert(N >= 0 && M >= 0);
        assert(values.count <= triangluar_data_size(N,M));
    }
    res : TriangularMatrix(T, N, M);
    for values {
        res.data[it_index] = it;
    }
    return res;
}


rows :: inline (m: TriangularMatrix, $check_transposed: bool = true) -> int {
    #if check_transposed {
        return ifx !m.transposed then m.N else m.M;
    } else {
        return m.N;
    }
}
cols :: inline (m: TriangularMatrix, $check_transposed: bool = true) -> int {
    #if check_transposed {
        return ifx !m.transposed then m.M else m.N;
    } else {
        return m.M;
    }
}

get :: inline (m: *TriangularMatrix, i: int, j: int, $check_transposed: bool = true) -> m.data_type {
    return inline get(<<m, i, j, check_transposed);
}
get :: inline (m: TriangularMatrix, i: int, j: int, $check_transposed: bool = true) -> m.data_type {
    #if check_transposed {
        if m.transposed {
            swap_code(i, j);
        } 
    }
    if i > j then return zero(m.data_type);
    return m.data[triangular_first_index_of_row(m.M, i) + (j-i)];
}
set :: inline (m: *TriangularMatrix, i: int, j: int, value: m.data_type, $check_transposed: bool = true) {
    #if check_transposed {
        if m.transposed {
            swap_code(i,j);
        }
    }
    #if CHECKS {
        assert( i <= j );
    }
    m.data[triangular_first_index_of_row(m.M,i)+(j-i)] = value;
}

get_lin :: inline (m: TriangularMatrix, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_lin :: inline (m: *TriangularMatrix, i: int, value: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, value);
}

get_raw :: inline (m: TriangularMatrix, i: int) -> m.data_type {
    return m.data[i];
}
set_raw :: inline (m: *TriangularMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}
n_raw :: inline (m: TriangularMatrix) -> int {
    return m.data.count;
}
from_raw_index :: inline (m: TriangularMatrix, i: int) -> int {
    #if CHECKS {
        assert(i >= 0 && i < m.data.count);
    }
    count := 0;
    r := 0;
    c := 0;
    for 0..m.N-1 {
        next := m.M - it;
        if count + next > i {
            r = it;
            c = i - count;
            break;
        }
    }
    return ifx !m.transposed then to_lin_index(m,r,c) else to_lin_index(m,c,r);
}




test_TriangularMatrix :: () {
    println_push("TriangularMatrix", color = .FG_GREEN);

    {
        println_push("test utils 1", color = .FG_YELLOW);
        assert( triangular_first_index_of_row(4, 2) == 7 );
        for 0..6-1 {
            println("%: %", it, triangular_first_index_of_row(6, it));
        }
    }

    {
        println_push("test utils 2", color = .FG_YELLOW);
        assert( triangluar_data_size(3, 5) == 12 );
        assert( triangluar_data_size(5, 3) == 6);
    }

    {
        println_push("test 1", color = .FG_YELLOW);
        tm : TriangularMatrix(int, 8, 5) = ---;
        for raw(*tm) {
            set_raw(*tm, it_raw_index, it_raw_index);
        }
        tm.transposed = false;
        println("TriangularMatrix: %", pstr(tm));
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        tm : TriangularMatrix(int, 8, 5) = ---;
        for raw(*tm) {
            set_raw(*tm, it_raw_index, it_raw_index);
        }
        tm.transposed = true;
        println("TriangularMatrixᵀ: %", pstr(tm));
    }

    {
        println_push("test 3", color = .FG_YELLOW);
        tm : TriangularMatrix(int, 8, 5);
        set(*tm, 1, 3, -42);
        set(*tm, 0, 4, 313);
        println("TriangularMatrix: %", pstr(tm));
        assert(get(tm, 1, 3) == -42);
        assert(get(tm, 0, 4) == 313);
    }

    {
        println_push("test 4", color = .FG_YELLOW);
        tm : TriangularMatrix(int, 8, 5);
        tm.transposed = true;
        set(*tm, 3, 1, -42);
        set(*tm, 4, 0, 313);
        println("TriangularMatrixᵀ : %", pstr(tm));
        assert(get(tm, 3, 1) == -42);
        assert(get(tm, 4, 0) == 313);
    }
}



// ######## ########  ####    ##     ## ########    ###    ########
//    ##    ##     ##  ##     ##     ## ##         ## ##   ##     ##
//    ##    ##     ##  ##     ##     ## ##        ##   ##  ##     ##
//    ##    ########   ##     ######### ######   ##     ## ########
//    ##    ##   ##    ##     ##     ## ##       ######### ##
//    ##    ##    ##   ##     ##     ## ##       ##     ## ##
//    ##    ##     ## ####    ##     ## ######## ##     ## ##


TriangularHeapMatrix :: struct(_T: Type) {
    using _m: MatrixType(_T, MATRIX_FLAGS.DENSE | MATRIX_FLAGS.TRIANGULAR | MATRIX_FLAGS.HEAP);
    transposed: bool = false;
    rows: int;
    cols: int;
    data: []_T;
}

free :: inline (m: TriangularHeapMatrix) {
    array_free(m.data);
}

thmat :: ($T: Type, N: int, M: int, transposed: bool = false, $init_zeroes: bool = true) -> TriangularHeapMatrix(T) {
    res : TriangularHeapMatrix(T) = ---;
    res.rows = N;
    res.cols = M;
    res.transposed = transposed;
    res.data = NewArray(triangluar_data_size(N, M), T, init_zeroes);
    return res;
}

thmat :: ($T: Type, N: int, M: int, values: ..$U) -> TriangularHeapMatrix(T) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        assert(N >= 0 && M >= 0);
        assert(values.count <= triangluar_data_size(N,M));
    }
    res := thmat(T, N, M);
    for values {
        res.data[it_index] = ncast(T, it);
    }
    return res;
}

thmat :: (N: int, M: int, values: ..$T) -> TriangularHeapMatrix(T) {
    #if CHECKS {
        assert(N >= 0 && M >= 0);
        assert(values.count <= triangluar_data_size(N,M));
    }
    res := thmat(T,N,M);
    for values {
        res.data[it_index] = it;
    }
    return res;
}

rows :: inline (m: TriangularHeapMatrix, $check_transposed: bool = true) -> int {
    #if check_transposed {
        return ifx !m.transposed then m.rows else m.cols;
    } else {
        return m.rows;
    }
}
cols :: inline (m: TriangularHeapMatrix, $check_transposed: bool = true) -> int {
    #if check_transposed {
        return ifx !m.transposed then m.cols else m.rows;
    } else {
        return m.cols;
    }
}

get :: inline (m: *TriangularHeapMatrix, i: int, j: int, $check_transposed: bool = true) -> m.data_type {
    return inline get(<<m, i, j, check_transposed);
}
get :: inline (m: TriangularHeapMatrix, i: int, j: int, $check_transposed: bool = true) -> m.data_type {
    #if check_transposed {
        if m.transposed {
            swap_code(i, j);
        } 
    }
    if i > j then return zero(m.data_type);
    return m.data[triangular_first_index_of_row(m.cols, i) + (j-i)];
}
set :: inline (m: *TriangularHeapMatrix, i: int, j: int, value: m.data_type, $check_transposed: bool = true) {
    #if check_transposed {
        if m.transposed {
            swap_code(i,j);
        }
    }
    #if CHECKS {
        assert( i <= j );
    }
    m.data[triangular_first_index_of_row(m.cols,i)+(j-i)] = value;
}

get_lin :: inline (m: TriangularHeapMatrix, i: int) -> m.data_type {
    r, c := from_lin_index(m, i);
    return get(m, r, c);
}
set_lin :: inline (m: *TriangularHeapMatrix, i: int, value: m.data_type) {
    r, c := from_lin_index(m, i);
    set(m, r, c, value);
}

get_raw :: inline (m: TriangularHeapMatrix, i: int) -> m.data_type {
    return m.data[i];
}
set_raw :: inline (m: *TriangularHeapMatrix, i: int, value: m.data_type) {
    m.data[i] = value;
}
n_raw :: inline (m: TriangularHeapMatrix) -> int {
    return m.data.count;
}
from_raw_index :: inline (m: TriangularHeapMatrix, i: int) -> int {
    #if CHECKS {
        assert(i >= 0 && i < m.data.count);
    }
    count := 0;
    r := 0;
    c := 0;
    for 0..m.rows-1 {
        next := m.cols - it;
        if count + next > i {
            r = it;
            c = i - count;
            break;
        }
    }
    return ifx !m.transposed then to_lin_index(m,r,c) else to_lin_index(m,c,r);
}


// ######## ########  ######  ########  ######
//    ##    ##       ##    ##    ##    ##    ##
//    ##    ##       ##          ##    ##
//    ##    ######    ######     ##     ######
//    ##    ##             ##    ##          ##
//    ##    ##       ##    ##    ##    ##    ##
//    ##    ########  ######     ##     ######

test_TriangularHeapMatrix :: () {
    println_push("TriangularHeapMatrix", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        tm := thmat(int, 8, 5);
        defer free(tm);
        for raw(*tm) {
            set_raw(*tm, it_raw_index, it_raw_index);
        }
        tm.transposed = false;
        println("TriangularHeapMatrix: %", pstr(tm));
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        tm := thmat(int, 8, 5);
        defer free(tm);
        for raw(*tm) {
            set_raw(*tm, it_raw_index, it_raw_index);
        }
        tm.transposed = true;
        println("TriangularHeapMatrix: %", pstr(tm));
    }

    {
        println_push("test 3", color = .FG_YELLOW);
        tm := thmat(int, 8, 5);
        defer free(tm);
        set(*tm, 1, 3, -42);
        set(*tm, 0, 4, 313);
        println("TriangularHeapMatrix: %", pstr(tm));
        assert(get(tm, 1, 3) == -42);
        assert(get(tm, 0, 4) == 313);
    }

    {
        println_push("test 4", color = .FG_YELLOW);
        tm := thmat(int, 8, 5);
        defer free(tm);
        tm.transposed = true;
        set(*tm, 3, 1, -42);
        set(*tm, 4, 0, 313);
        println("TriangularHeapMatrix : %", pstr(tm));
        assert(get(tm, 3, 1) == -42);
        assert(get(tm, 4, 0) == 313);
    }
}

test_MatrixTypes :: () {
    println_push("MatrixTypes", color = .FG_WHITE);
    test_BandMatrix();
    // test_BandHeapMatrix();
    // test_TriangularMatrix();
    // test_TriangularHeapMatrix();
}