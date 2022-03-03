
// ######## ##          ###     ######    ######
// ##       ##         ## ##   ##    ##  ##    ##
// ##       ##        ##   ##  ##        ##
// ######   ##       ##     ## ##   ####  ######
// ##       ##       ######### ##    ##        ##
// ##       ##       ##     ## ##    ##  ##    ##
// ##       ######## ##     ##  ######    ######

VECTOR_FLAGS :: enum_flags u8 {
    VIEW  :: 0;
    DENSE :: (1 << 0);
    HEAP  :: (1 << 1);
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
* All concrete types A(T) of ArrayType1D(T) have to implement the following functions:
* get: (A(T), int) -> T
* set: (*A(T), int, T)
* dim: (A) -> int
*/
ArrayType1D :: struct(data_type: Type) {}

// get :: inline (a: *$A/ArrayType1D, i: int) -> a.data_type {
//     return get(<<a, i);
// }

// dim :: inline (a: *$A/ArrayType1D) -> int {
//     return dim(<<a);
// }

// operator [] :: inline (a: $A/ArrayType1D, i: int) -> A.data_type {
//     return get(a, i);
// }
// operator [] :: inline (a: *$A/ArrayType1D, i: int) -> A.data_type {
//     return get(<<a, i);
// }

//    ###    ########  ########     ###    ##    ##    ##   ########
//   ## ##   ##     ## ##     ##   ## ##    ##  ##   ####   ##     ##
//  ##   ##  ##     ## ##     ##  ##   ##    ####      ##   ##     ##
// ##     ## ########  ########  ##     ##    ##       ##   ##     ##
// ######### ##   ##   ##   ##   #########    ##       ##   ##     ##
// ##     ## ##    ##  ##    ##  ##     ##    ##       ##   ##     ##
// ##     ## ##     ## ##     ## ##     ##    ##     ###### ########

// * implementation for simple arrays
Array1D :: struct(N: int, _T: Type) {
    using _a1d: ArrayType1D(_T);
    data: [N]_T;
}
array1d :: (vals: [$N]$T) -> Array1D(N,T) {
    res : Array1D(N,T) = ---;
    for vals {
        res.data[it_index] = it;
    }
    return res;
}
array1d :: ($N: int, vals: ..$T) -> Array1D(N,T) {
    res : Array1D(N,T) = ---;
    for vals {
        res.data[it_index] = it;
    }
    return res;
}
get :: inline (a: Array1D, i: int) -> a.data_type {
    return a.data[i];
}
set :: inline (a: *Array1D, i: int, val: a.data_type) {
    a.data[i] = val;
}

operator == :: (a: Array1D, b: Array1D) -> bool {
    #if a.N != b.N { 
        return false; 
    } else {
        for a.data {
            if it != b.data[it_index] then return false;
        }
        return true;
    }
}

//    ###          ##   ########     ##     ## #### ######## ##      ##
//   ## ##       ####   ##     ##    ##     ##  ##  ##       ##  ##  ##
//  ##   ##        ##   ##     ##    ##     ##  ##  ##       ##  ##  ##
// ##     ##       ##   ##     ##    ##     ##  ##  ######   ##  ##  ##
// #########       ##   ##     ##     ##   ##   ##  ##       ##  ##  ##
// ##     ##       ##   ##     ##      ## ##    ##  ##       ##  ##  ##
// ##     ##     ###### ########        ###    #### ########  ###  ###

Array1DView :: struct(_A: Type) {
    using _a1d: ArrayType1D(_A.data_type);
    a2d : *_A;
    cols: int;
}

array1d_view :: (a2d: *$A/ArrayType2D, cols: int) -> Array1DView(A) {
    res : Array1DView(A) = ---;
    res.a2d = a2d;
    res.cols = cols;
    return res;
}

get :: inline (a: Array1DView, i: int) -> a.data_type {
    col := i % a.cols;
    row := i / a.cols;
    return get(a.a2d, row, col);
}
set :: inline (a: *Array1DView, i: int, val: a.data_type) {
    col := i % a.cols;
    row := i / a.cols;
    set(a.a2d, row, col, val);
}


// ##     ## ########    ###    ########     ###    ########  ########     ###    ##    ##    ##   ########
// ##     ## ##         ## ##   ##     ##   ## ##   ##     ## ##     ##   ## ##    ##  ##   ####   ##     ##
// ##     ## ##        ##   ##  ##     ##  ##   ##  ##     ## ##     ##  ##   ##    ####      ##   ##     ##
// ######### ######   ##     ## ########  ##     ## ########  ########  ##     ##    ##       ##   ##     ##
// ##     ## ##       ######### ##        ######### ##   ##   ##   ##   #########    ##       ##   ##     ##
// ##     ## ##       ##     ## ##        ##     ## ##    ##  ##    ##  ##     ##    ##       ##   ##     ##
// ##     ## ######## ##     ## ##        ##     ## ##     ## ##     ## ##     ##    ##     ###### ########

// * implementation for simple arrays
HeapArray1D :: struct(_T: Type) {
    using _a1d: ArrayType1D(_T);
    data: []_T;
}
harray1d :: (vals: []$T) -> HeapArray1D(T) {
    res : HeapArray1D(T) = ---;
    res.data = NewArray(vals.count, T, false); // ! alloc
    for vals {
        res.data[it_index] = it;
    }
    return res;
}
harray1d :: (val: ..$T) -> HeapArray1D(T) {
    res : HeapArray1D(T) = ---;
    res.data = NewArray(val.count, T, false); // ! alloc
    for val {
        res.data[it_index] = it;
    }
    return res;
}
get :: inline (a: HeapArray1D, i: int) -> a.data_type {
    return a.data[i];
}
set :: inline (a: *HeapArray1D, i: int, val: a.data_type) {
    a.data[i] = val;
}

operator == :: (a: HeapArray1D, b: HeapArray1D) -> bool {
    if a.data.count != b.data.count then return false;
    for a.data {
        if it != b.data[it_index] then return false;
    }
    return true;
}


// ##     ## ########  ######     ######## ##    ## ########  ########
// ##     ## ##       ##    ##       ##     ##  ##  ##     ## ##
// ##     ## ##       ##             ##      ####   ##     ## ##
// ##     ## ######   ##             ##       ##    ########  ######
//  ##   ##  ##       ##             ##       ##    ##        ##
//   ## ##   ##       ##    ##       ##       ##    ##        ##
//    ###    ########  ######        ##       ##    ##        ########


/*
* All concrete types V(T) of VectorType(T,flags) have to implement the following functions:
* 
* dim: (V(T)) -> int
* 
* get: (V(T), int) -> T
* set: (*V(T), int, T) -> ()
* 
* n_raw: (V(T)) -> int
* get_raw: (V(T), int) -> T
* set_raw: (*V(T), int, T) -> ()
* from_raw_index: (V(T), int) -> int 
*/
VectorType :: struct(
    _data_type: Type,
    flags: VECTOR_FLAGS
){
    using _a1d: ArrayType1D(_data_type);
}

from_raw_index :: inline (v: *$V/VectorType, i: int) -> int {
    return from_raw_index(<<v, i);
}

// ########  #######  ########          ######## ##     ## ########  ##    ##
// ##       ##     ## ##     ##         ##        ##   ##  ##     ## ###   ##
// ##       ##     ## ##     ##         ##         ## ##   ##     ## ####  ##
// ######   ##     ## ########  ####### ######      ###    ########  ## ## ##
// ##       ##     ## ##   ##           ##         ## ##   ##        ##  ####
// ##       ##     ## ##    ##          ##        ##   ##  ##        ##   ###
// ##        #######  ##     ##         ######## ##     ## ##        ##    ##


// * linear for expansion over all elements
for_expansion :: (v: *$V/VectorType, body: Code, flags: For_Flags) #expand {
    #assert(flags == 0);
    i := 0;
    while i < dim(v) {
        `it_index := i;
        `it := get(v, i);
        #insert (continue={i += 1; continue;}, break={break;}) body;
        i += 1;
    }
}


raw :: inline (v: *$V/VectorType) -> *VectorRawForExpansion(V) {
    res := New(VectorRawForExpansion(V), false);
    res.vector = v;
    return res;
}
#scope_file
VectorRawForExpansion :: struct(T: Type) {
    vector: *T;
}
#scope_export
// * for expansion only over the really existing data
for_expansion :: (ve: *VectorRawForExpansion, body: Code, flags: For_Flags) #expand {
    defer free(ve);
    v := <<(ve.vector);
    i := 0;
    while i < n_raw(v) {
        `it_index := from_raw_index(v, i);
        `it := get_raw(v, i);
        `it_raw_index := i;
        #insert (continue={i += 1; continue;}, break={break;}) body;
        i += 1;
    }
}

// ##     ## ########  ######     ######## ##    ## ########  ########  ######
// ##     ## ##       ##    ##       ##     ##  ##  ##     ## ##       ##    ##
// ##     ## ##       ##             ##      ####   ##     ## ##       ##
// ##     ## ######   ##             ##       ##    ########  ######    ######
//  ##   ##  ##       ##             ##       ##    ##        ##             ##
//   ## ##   ##       ##    ##       ##       ##    ##        ##       ##    ##
//    ###    ########  ######        ##       ##    ##        ########  ######


// ##     ##    ###    ########    ########   #######  ##      ##    ##     ## #### ######## ##      ##
// ###   ###   ## ##      ##       ##     ## ##     ## ##  ##  ##    ##     ##  ##  ##       ##  ##  ##
// #### ####  ##   ##     ##       ##     ## ##     ## ##  ##  ##    ##     ##  ##  ##       ##  ##  ##
// ## ### ## ##     ##    ##       ########  ##     ## ##  ##  ##    ##     ##  ##  ######   ##  ##  ##
// ##     ## #########    ##       ##   ##   ##     ## ##  ##  ##     ##   ##   ##  ##       ##  ##  ##
// ##     ## ##     ##    ##       ##    ##  ##     ## ##  ##  ##      ## ##    ##  ##       ##  ##  ##
// ##     ## ##     ##    ##       ##     ##  #######   ###  ###        ###    #### ########  ###  ###

MatrixRowView :: struct(M: Type, N: int) {
    // using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW & M.flags); // TODO implement flags conversion
    using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW);
    row: int;
    start_index: int = 0;
    stride: int = 1;
    mat: *M;
}


mat_row_view :: (mat: *$M/MatrixType, row: int, $num_elements: int, start_index: int = 0, stride: int = 1) -> MatrixRowView(M, num_elements) {
    #if CHECKS {
        assert(row >= 0 && row < rows(mat));
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : MatrixRowView(M, num_elements);
    res.row = row;
    res.start_index = start_index;
    res.stride = stride;
    res.mat = mat;
    return res;
}

dim :: inline (v: MatrixRowView) -> int {
    return v.N;
}
get :: inline (v: MatrixRowView, i: int) -> v.data_type {
    return get(v.mat, v.row, v.start_index + v.stride * i);
}
set :: inline (v: *MatrixRowView, i: int, value: v.data_type) {
    set(v.mat, v.row, v.start_index + v.stride * i, value);
}
// TODO be smarter about this?
get_raw :: inline (v: MatrixRowView, i: int) -> v.data_type {
    return get(v.mat, v.row, v.start_index + v.stride * i);
}
set_raw :: inline (v: *MatrixRowView, i: int, value: v.data_type) {
    set(v.mat, v.row, v.start_index + v.stride * i, value);
}
n_raw :: inline (v: MatrixRowView) -> int {
    return v.N; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
}
from_raw_index :: inline (v: MatrixRowView, i: int) -> int {
    return i;
}

// ##     ##    ###    ########     ######   #######  ##          ##     ## #### ######## ##      ##
// ###   ###   ## ##      ##       ##    ## ##     ## ##          ##     ##  ##  ##       ##  ##  ##
// #### ####  ##   ##     ##       ##       ##     ## ##          ##     ##  ##  ##       ##  ##  ##
// ## ### ## ##     ##    ##       ##       ##     ## ##          ##     ##  ##  ######   ##  ##  ##
// ##     ## #########    ##       ##       ##     ## ##           ##   ##   ##  ##       ##  ##  ##
// ##     ## ##     ##    ##       ##    ## ##     ## ##            ## ##    ##  ##       ##  ##  ##
// ##     ## ##     ##    ##        ######   #######  ########       ###    #### ########  ###  ###

MatrixColumnView :: struct(M: Type, N: int) {
    // using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW & M.flags); // TODO implement flags conversion
    using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW);
    col: int;
    start_index: int = 0;
    stride: int = 1;
    mat: *M;
}


mat_col_view :: (mat: *$M/MatrixType, col: int, $num_elements: int, start_index: int = 0, stride: int = 1) -> MatrixColumnView(M, num_elements) {
    #if CHECKS {
        assert(col >= 0 && col < cols(mat));
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : MatrixColumnView(M, num_elements);
    res.col = col;
    res.start_index = start_index;
    res.stride = stride;
    res.mat = mat;
    return res;
}

dim :: inline (v: MatrixColumnView) -> int {
    return v.N;
}
get :: inline (v: MatrixColumnView, i: int) -> v.data_type {
    return get(v.mat, v.start_index + v.stride * i, v.col);
}
set :: inline (v: *MatrixColumnView, i: int, value: v.data_type) {
    set(v.mat, v.start_index + v.stride * i, v.col, value);
}
// TODO be smarter about this?
get_raw :: inline (v: MatrixColumnView, i: int) -> v.data_type {
    return get(v.mat, v.start_index + v.stride * i, v.col);
}
set_raw :: inline (v: *MatrixColumnView, i: int, value: v.data_type) {
    set(v.mat, v.start_index + v.stride * i, value, v.col);
}
n_raw :: inline (v: MatrixColumnView) -> int {
    return v.N; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
}
from_raw_index :: inline (v: MatrixColumnView, i: int) -> int {
    return i;
}


// ##     ##    ###    ########    ########   #######  ##      ##    ##     ## ########    ###    ########     ##     ## #### ######## ##      ##
// ###   ###   ## ##      ##       ##     ## ##     ## ##  ##  ##    ##     ## ##         ## ##   ##     ##    ##     ##  ##  ##       ##  ##  ##
// #### ####  ##   ##     ##       ##     ## ##     ## ##  ##  ##    ##     ## ##        ##   ##  ##     ##    ##     ##  ##  ##       ##  ##  ##
// ## ### ## ##     ##    ##       ########  ##     ## ##  ##  ##    ######### ######   ##     ## ########     ##     ##  ##  ######   ##  ##  ##
// ##     ## #########    ##       ##   ##   ##     ## ##  ##  ##    ##     ## ##       ######### ##            ##   ##   ##  ##       ##  ##  ##
// ##     ## ##     ##    ##       ##    ##  ##     ## ##  ##  ##    ##     ## ##       ##     ## ##             ## ##    ##  ##       ##  ##  ##
// ##     ## ##     ##    ##       ##     ##  #######   ###  ###     ##     ## ######## ##     ## ##              ###    #### ########  ###  ###

MatrixRowHeapView :: struct(M: Type) {
    // using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW & M.flags); // TODO implement flags conversion
    using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW);
    row: int;
    cols: int;
    start_index: int = 0;
    stride: int = 1;
    mat: *M;
}


mat_row_heap_view :: (mat: *$M/MatrixType, row: int, num_elements: int, start_index: int = 0, stride: int = 1) -> MatrixRowHeapView(M) {
    #if CHECKS {
        assert(row >= 0 && row < rows(mat));
        assert(start_index + (num_elements-1) * stride <= cols(mat));
    }
    res : MatrixRowHeapView(M);
    res.row = row;
    res.cols = num_elements;
    res.start_index = start_index;
    res.stride = stride;
    res.mat = mat;
    return res;
}

dim :: inline (v: MatrixRowHeapView) -> int {
    return v.cols;
}
get :: inline (v: MatrixRowHeapView, i: int) -> v.data_type {
    return get(v.mat, v.row, v.start_index + v.stride * i);
}
set :: inline (v: *MatrixRowHeapView, i: int, value: v.data_type) {
    set(v.mat, v.row, v.start_index + v.stride * i, value);
}
// TODO be smarter about this?
get_raw :: inline (v: MatrixRowHeapView, i: int) -> v.data_type {
    return get(v.mat, v.row, v.start_index + v.stride * i);
}
set_raw :: inline (v: *MatrixRowHeapView, i: int, value: v.data_type) {
    set(v.mat, v.row, v.start_index + v.stride * i, value);
}
n_raw :: inline (v: MatrixRowHeapView) -> int {
    return v.cols; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
}
from_raw_index :: inline (v: MatrixRowHeapView, i: int) -> int {
    return i;
}

// ##     ##    ###    ########     ######   #######  ##          ##     ## ########    ###    ########     ##     ## #### ######## ##      ##
// ###   ###   ## ##      ##       ##    ## ##     ## ##          ##     ## ##         ## ##   ##     ##    ##     ##  ##  ##       ##  ##  ##
// #### ####  ##   ##     ##       ##       ##     ## ##          ##     ## ##        ##   ##  ##     ##    ##     ##  ##  ##       ##  ##  ##
// ## ### ## ##     ##    ##       ##       ##     ## ##          ######### ######   ##     ## ########     ##     ##  ##  ######   ##  ##  ##
// ##     ## #########    ##       ##       ##     ## ##          ##     ## ##       ######### ##            ##   ##   ##  ##       ##  ##  ##
// ##     ## ##     ##    ##       ##    ## ##     ## ##          ##     ## ##       ##     ## ##             ## ##    ##  ##       ##  ##  ##
// ##     ## ##     ##    ##        ######   #######  ########    ##     ## ######## ##     ## ##              ###    #### ########  ###  ###


MatrixColumnHeapView :: struct(M: Type) {
    // using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW & M.flags); // TODO implement flags conversion
    using _v: VectorType(M.data_type, VECTOR_FLAGS.VIEW);
    col: int;
    rows: int;
    start_index: int = 0;
    stride: int = 1;
    mat: *M;
}


mat_col_heap_view :: (mat: *$M/MatrixType, col: int, num_elements: int, start_index: int = 0, stride: int = 1) -> MatrixColumnHeapView(M) {
    #if CHECKS {
        assert(col >= 0 && col < cols(mat));
        assert(start_index + (num_elements-1) * stride <= rows(mat));
    }
    res : MatrixColumnHeapView(M);
    res.col = col;
    res.rows = num_elements;
    res.start_index = start_index;
    res.stride = stride;
    res.mat = mat;
    return res;
}

dim :: inline (v: MatrixColumnHeapView) -> int {
    return v.rows;
}
get :: inline (v: MatrixColumnHeapView, i: int) -> v.data_type {
    return get(v.mat, v.start_index + v.stride * i, v.col);
}
set :: inline (v: *MatrixColumnHeapView, i: int, value: v.data_type) {
    set(v.mat, v.start_index + v.stride * i, v.col, value);
}
// TODO be smarter about this?
get_raw :: inline (v: MatrixColumnHeapView, i: int) -> v.data_type {
    return get(v.mat, v.start_index + v.stride * i, v.col);
}
set_raw :: inline (v: *MatrixColumnHeapView, i: int, value: v.data_type) {
    set(v.mat, v.start_index + v.stride * i, value, v.col);
}
n_raw :: inline (v: MatrixColumnHeapView) -> int {
    return v.rows; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
}
from_raw_index :: inline (v: MatrixColumnHeapView, i: int) -> int {
    return i;
}


// ##     ## #### ######## ##      ##
// ##     ##  ##  ##       ##  ##  ##
// ##     ##  ##  ##       ##  ##  ##
// ##     ##  ##  ######   ##  ##  ##
//  ##   ##   ##  ##       ##  ##  ##
//   ## ##    ##  ##       ##  ##  ##
//    ###    #### ########  ###  ###

VectorView :: struct(V: Type, N: int) {
    using _v: VectorType(V.data_type, VECTOR_FLAGS.VIEW & V.flags);
    start_index: int = 0;
    stride: int = 1;
    vec: *V;
}


vec_view :: (vec: *$V/VectorType, $num_elements: int, start_index: int = 0, stride: int = 1) -> VectorView(V, num_elements) {
    #if CHECKS {
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : VectorView(V, num_elements);
    res.start_index = start_index;
    res.stride = stride;
    res.vec = vec;
    return res;
}

dim :: inline (v: VectorView) -> int {
    return v.N;
}
get :: inline (v: VectorView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
set :: inline (v: *VectorView, i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
// TODO be smarter about this?
get_raw :: inline (v: VectorView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
set_raw :: inline (v: *VectorView, i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
n_raw :: inline (v: VectorView) -> int {
    return v.N; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
}
from_raw_index :: inline (v: VectorView, i: int) -> int {
    return i;
}



// ##     ## ########    ###    ########     ##     ## #### ######## ##      ##
// ##     ## ##         ## ##   ##     ##    ##     ##  ##  ##       ##  ##  ##
// ##     ## ##        ##   ##  ##     ##    ##     ##  ##  ##       ##  ##  ##
// ######### ######   ##     ## ########     ##     ##  ##  ######   ##  ##  ##
// ##     ## ##       ######### ##            ##   ##   ##  ##       ##  ##  ##
// ##     ## ##       ##     ## ##             ## ##    ##  ##       ##  ##  ##
// ##     ## ######## ##     ## ##              ###    #### ########  ###  ###

VectorHeapView :: struct(V: Type) {
    using _v: VectorType(V.data_type, VECTOR_FLAGS.VIEW & V.flags);
    dims: int;
    start_index: int = 0;
    stride: int = 1;
    vec: *V;
}


vec_heap_view :: inline(vec: *$V/VectorType, num_elements: int, start_index: int = 0, stride: int = 1) -> VectorHeapView(V) {
    #if CHECKS {
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : VectorHeapView(V);
    res.dims = num_elements;
    res.start_index = start_index;
    res.stride = stride;
    res.vec = vec;
    return res;
}


dim :: inline (v: VectorHeapView) -> int {
    return v.dims;
}
get :: inline (v: VectorHeapView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
set :: inline (v: *VectorHeapView, i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
get_raw :: inline (v: VectorHeapView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
set_raw :: inline (v: *VectorHeapView, i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
n_raw :: inline (v: VectorHeapView) -> int {
    return v.dims;
}
from_raw_index :: inline (v: VectorHeapView, i: int) -> int {
    return i;
}


// ########  ######## ##    ##  ######  ########
// ##     ## ##       ###   ## ##    ## ##
// ##     ## ##       ####  ## ##       ##
// ##     ## ######   ## ## ##  ######  ######
// ##     ## ##       ##  ####       ## ##
// ##     ## ##       ##   ### ##    ## ##
// ########  ######## ##    ##  ######  ########

DenseVector :: struct(_T: Type, N: int) {
    using _v: VectorType(_T, VECTOR_FLAGS.DENSE);
    data: [N]_T;
}



dvec :: inline($N: int, values: ..$T) -> DenseVector(T,N) {
    #if CHECKS {
        assert(values.count <= N);
    }
    res : DenseVector(T,N);
    for values {
        set(*res, it_index, it);
    }
    return res;
}

dvec :: ($T: Type, $N: int, values: ..$U) -> DenseVector(T,N) {
    #if CHECKS {
        #assert(is_minor_number_type(T, U));
        assert(values.count <= N);
    }
    res : DenseVector(T,N);
    for values {
        set(*res, it_index, ncast(T, it));
    }
    return res;
}

dim :: inline (v: DenseVector) -> int {
    return v.N;
}
get :: inline (v: DenseVector, i: int) -> v.data_type {
    return v.data[i];
}
set :: inline (v: *DenseVector, i: int, value: v.data_type) {
    v.data[i] = value;
}
get_raw :: inline (v: DenseVector, i: int) -> v.data_type {
    return v.data[i];
}
set_raw :: inline (v: *DenseVector, i: int, value: v.data_type) {
    v.data[i] = value;
}
n_raw :: inline (v: DenseVector) -> int {
    return v.N;
}
from_raw_index :: inline (v: DenseVector, i: int) -> int {
    return i;
}




// ########  ######## ##    ##  ######  ########    ##     ## ########    ###    ########
// ##     ## ##       ###   ## ##    ## ##          ##     ## ##         ## ##   ##     ##
// ##     ## ##       ####  ## ##       ##          ##     ## ##        ##   ##  ##     ##
// ##     ## ######   ## ## ##  ######  ######      ######### ######   ##     ## ########
// ##     ## ##       ##  ####       ## ##          ##     ## ##       ######### ##
// ##     ## ##       ##   ### ##    ## ##          ##     ## ##       ##     ## ##
// ########  ######## ##    ##  ######  ########    ##     ## ######## ##     ## ##


DenseHeapVector :: struct(_T: Type) {
    using _v: VectorType(_T, VECTOR_FLAGS.DENSE & VECTOR_FLAGS.HEAP);
    dims : int;
    data: []_T;
}

free :: inline (v: $V/DenseHeapVector) {
    array_free(v.data);
}

dhvec :: inline ($T: Type, N: int, $initialize_variables: bool = true) -> DenseHeapVector(T) {
    res : DenseHeapVector(T);
    res.dims = N;
    init(*res, initialize_variables);
    return res;
}

dhvec :: inline(N: int, values: ..$T) -> DenseHeapVector(T) {
    #if CHECKS {
        assert(values.count <= N);
    }
    res := make(DenseHeapVector(T), false);
    vec_init(*res, ..values);
    return res;
}

dhvec :: ($T: Type, N: int, values: ..$U) -> DenseHeapVector(T) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        assert(values.count <= N);
    }
    res : DenseHeapVector(T);
    res.dims = N;
    res.data = NewArray(N, T, false);
    for 0..N-1 {
        if it < values.count {
            res.data[it] = ncast(T, values[it]);
        } else {
            res.data[it] = zero(T);
        }
    }
    return res;
}

init :: inline (v: *DenseHeapVector, $initialize_variables: bool = true) {
    v.data = NewArray(v.dims, v.data_type, false);
    #if initialize_variables {
        for 0..v.dims-1 {
            v.data[it] = zero(v.data_type);
        }
    }
}

dim :: inline (v: DenseHeapVector) -> int {
    return v.dims;
}
get :: inline (v: DenseHeapVector, i: int) -> v.data_type {
    return v.data[i];
}
set :: inline (v: *DenseHeapVector, i: int, value: v.data_type) {
    v.data[i] = value;
}
get_raw :: inline (v: DenseHeapVector, i: int) -> v.data_type {
    return v.data[i];
}
set_raw :: inline (v: *DenseHeapVector, i: int, value: v.data_type) {
    v.data[i] = value;
}
n_raw :: inline (v: DenseHeapVector) -> int {
    return v.dims;
}
from_raw_index :: inline (v: DenseHeapVector, i: int) -> int {
    return i;
}