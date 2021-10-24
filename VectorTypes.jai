
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
All concrete types of VectorType have to implement the following functions:

dim: (v: vector_type) -> int

get: (v: vector_type, int) -> data_type
set: (v: *vector_type, int, data_type) -> ()

n_raw: (v: vector_type) -> int
get_raw: (v: vector_type, int) -> data_type
set_raw: (v: *vector_type, int, data_type) -> ()
from_raw_index: (v: vector_type, int) -> int

*/
VectorType :: struct(
    data_type: Type,
    flags: VECTOR_FLAGS
){}

get :: inline (v: *$V/VectorType, i: int) -> V.data_type {
    return get(<<v, i);
}
dim :: inline (v: *$V/VectorType) -> int {
    return dim(<<v);
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
    i := 0;
    while i < dim(v) {
        `it_index := i;
        `it := get(v, i);
        #insert body;
        i += 1;
    }
}


raw :: inline (v: *$V/VectorType) -> VectorRawForExpansion(V) {
    res : VectorRawForExpansion(V) = ---;
    res.vector = v;
    return res;
}
raw :: inline (v: $V/VectorType) -> VectorRawForExpansion(V) {
    res : VectorRawForExpansion(V) = ---;
    res.vector = *v;
    return res;
}
#scope_file
VectorRawForExpansion :: struct(T: Type) {
    vector: *T;
}
#scope_export
// * for expansion only over the really existing data
for_expansion :: (ve: *VectorRawForExpansion($V), body: Code, flags: For_Flags) #expand {
    v := <<(ve.vector);
    i := 0;
    while i < n_raw(v) {
        `it_index := from_raw_index(v, i);
        `it := get_raw(v, i);
        `it_raw_index := i;
        #insert body;
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

dim :: inline (v: MatrixRowView($M, $N)) -> int {
    return N;
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
n_raw :: inline (v: MatrixRowView($T, $N)) -> int {
    return N; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
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

dim :: inline (v: MatrixColumnView($M, $N)) -> int {
    return N;
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
n_raw :: inline (v: MatrixColumnView($T, $N)) -> int {
    return N; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
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
        assert(start_index + (num_elements-1) * stride <= dim(vec));
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
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : MatrixColumnHeapView(M, num_elements);
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

dim :: inline (v: VectorView($V, $N)) -> int {
    return N;
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
n_raw :: inline (v: VectorView($T, $N)) -> int {
    return N; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
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

HeapVectorView :: struct(V: Type) {
    using _v: VectorType(V.data_type, VECTOR_FLAGS.VIEW & V.flags);
    dims: int;
    start_index: int = 0;
    stride: int = 1;
    vec: *V;
}


vec_heap_view :: inline(vec: *$V/VectorType, num_elements: int, start_index: int = 0, stride: int = 1) -> HeapVectorView(V) {
    #if CHECKS {
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : HeapVectorView(V);
    res.dims = num_elements;
    res.start_index = start_index;
    res.stride = stride;
    res.vec = vec;
    return res;
}


dim :: inline (v: HeapVectorView) -> int {
    return v.dims;
}
get :: inline (v: HeapVectorView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
set :: inline (v: *HeapVectorView($T), i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
get_raw :: inline (v: HeapVectorView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
set_raw :: inline (v: *HeapVectorView($T), i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
n_raw :: inline (v: HeapVectorView) -> int {
    return v.dims;
}
from_raw_index :: inline (v: HeapVectorView, i: int) -> int {
    return i;
}


// TODO flags for array views?
// //    ###    ########  ########     ###    ##    ##    ##     ## #### ######## ##      ##
// //   ## ##   ##     ## ##     ##   ## ##    ##  ##     ##     ##  ##  ##       ##  ##  ##
// //  ##   ##  ##     ## ##     ##  ##   ##    ####      ##     ##  ##  ##       ##  ##  ##
// // ##     ## ########  ########  ##     ##    ##       ##     ##  ##  ######   ##  ##  ##
// // ######### ##   ##   ##   ##   #########    ##        ##   ##   ##  ##       ##  ##  ##
// // ##     ## ##    ##  ##    ##  ##     ##    ##         ## ##    ##  ##       ##  ##  ##
// // ##     ## ##     ## ##     ## ##     ##    ##          ###    #### ########  ###  ###

// dim :: inline (v: []$T) -> int {
//     return v.count;
// }

// get :: inline (v: []$T, i: int) -> T {
//     return v[i];
// }
// set :: inline (v: *[]$T, i: int, value: T) {
//     v[i] = value;
// }

// get_raw :: inline (v: []$T, i: int) -> T {
//     return v[i];
// }
// set_raw :: inline (v: *[]$T, i: int, value: T) {
//     v[i] = value;
// }
// n_raw :: inline (v: []$T) -> int {
//     return v.count;
// }
// from_raw_index :: inline (v: []$T, i: int) -> int {
//     return i;
// }


// ########  ######## ##    ##  ######  ########
// ##     ## ##       ###   ## ##    ## ##
// ##     ## ##       ####  ## ##       ##
// ##     ## ######   ## ## ##  ######  ######
// ##     ## ##       ##  ####       ## ##
// ##     ## ##       ##   ### ##    ## ##
// ########  ######## ##    ##  ######  ########

DenseVector :: struct(_T: Type, _N: int) {
    using _v: VectorType(_T, VECTOR_FLAGS.DENSE);
    data: [_N]_T;
}

dvec :: inline($N: int, values: ..$T) -> DenseVector(T,N) {
    #if CHECKS {
        assert(values.count <= N);
    }
    res : DenseVector(T,N);
    vec_init(*res, ..values);
    return res;
}

dvec :: ($T: Type, $N: int, values: ..$U) -> DenseVector(T,N) {
    #if CHECKS {
        #run assert(is_minor_number_type(T, U));
        assert(values.count <= N);
    }
    res : DenseVector(T,N);
    for values {
        set(*res, it_index, ncast(T, it));
    }
    return res;
}

dim :: inline (v: DenseVector($T, $N)) -> int {
    return N;
}
get :: inline (v: DenseVector($T, $N), i: int) -> T {
    return v.data[i];
}
set :: inline (v: *DenseVector($T, $N), i: int, value: T) {
    v.data[i] = value;
}
get_raw :: inline (v: DenseVector($T, $N), i: int) -> T {
    return v.data[i];
}
set_raw :: inline (v: *DenseVector($T, $N), i: int, value: T) {
    v.data[i] = value;
}
n_raw :: inline (v: DenseVector($T, $N)) -> int {
    return N;
}
from_raw_index :: inline (v: DenseVector($T, $N), i: int) -> int {
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

init :: inline (v: *DenseHeapVector($T), $initialize_variables: bool = true) {
    v.data = NewArray(v.dims, T, false);
    #if initialize_variables {
        for 0..v.dims-1 {
            v.data[it] = zero(T);
        }
    }
}

dim :: inline (v: DenseHeapVector($T)) -> int {
    return v.dims;
}
get :: inline (v: DenseHeapVector($T), i: int) -> T {
    return v.data[i];
}
set :: inline (v: *DenseHeapVector($T), i: int, value: T) {
    v.data[i] = value;
}
get_raw :: inline (v: DenseHeapVector($T), i: int) -> T {
    return v.data[i];
}
set_raw :: inline (v: *DenseHeapVector($T), i: int, value: T) {
    v.data[i] = value;
}
n_raw :: inline (v: DenseHeapVector($T)) -> int {
    return v.dims;
}
from_raw_index :: inline (v: DenseHeapVector($T), i: int) -> int {
    return i;
}