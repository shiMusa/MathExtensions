
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

VectorType :: struct(
    vector_type: Type, 
    data_type: Type,
    dim: (v: vector_type) -> int,
    flags: VECTOR_FLAGS,

    get: (v: vector_type, int) -> data_type,
    set: (v: *vector_type, int, data_type) -> (),

    n_raw: (v: vector_type) -> int,
    get_raw: (v: vector_type, int) -> data_type,
    set_raw: (v: *vector_type, int, data_type) -> (),
    from_raw_index: (v: vector_type, int) -> int
){}

dim :: inline (v: $V/VectorType) -> int {
    return V.dim(v);
}
dim :: inline (v: *$V/VectorType) -> int {
    return V.dim(<<v);
}

get :: inline (v: $V/VectorType, i: int) -> V.data_type {
    return V.get(v, i);
}
get :: inline (v: *$V/VectorType, i: int) -> V.data_type {
    return V.get(<<v, i);
}
set :: inline (v: *$V/VectorType, i: int, value: V.data_type) {
    V.set(v, i, value);
}


n_raw :: inline (v: $V/VectorType) -> int {
    return V.n_raw(v);
}
get_raw :: inline (v: $V/VectorType, i: int) -> V.data_type {
    return V.get_raw(v, i);
}
set_raw :: inline (v: *$V/VectorType, i: int, value: V.data_type) {
    V.set_raw(v, i, value);
}
from_raw_index :: inline (v: $V/VectorType, i: int) -> int {
    return V.from_raw_index(v, i);
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
    v := <<ve.vector;
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


// ##     ## #### ######## ##      ##
// ##     ##  ##  ##       ##  ##  ##
// ##     ##  ##  ##       ##  ##  ##
// ##     ##  ##  ######   ##  ##  ##
//  ##   ##   ##  ##       ##  ##  ##
//   ## ##    ##  ##       ##  ##  ##
//    ###    #### ########  ###  ###

VectorView :: struct(V: Type, N: int) {
    using _v: VectorType(
        #this,
        V.data_type,
        vector_view_dim,
        VECTOR_FLAGS.VIEW & V.flags,

        vector_view_init,

        vector_view_get,
        vector_view_set,

        vector_view_n_raw,
        vector_view_get, // same as above for dense vectors
        vector_view_set, // same as above for dense vectors
        vector_view_raw_index
    );

    start_index: int = 0;
    stride: int = 1;
    vec: *V;
}


vecview :: (vec: *$V/VectorType, start_index: int = 0, $num_elements: int, stride: int = 1) -> VectorView(V, num_elements) {
    #if CHECKS {
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : VectorView(V, num_elements);
    res.start_index = start_index;
    res.stride = stride;
    res.vec = vec;
    return res;
}

#scope_file
vector_view_dim :: inline (v: VectorView($V, $N)) -> int {
    return N;
}
vector_view_n_raw :: inline (v: VectorView($T, $N)) -> int {
    return N; // ! not correct since e.g. a view on a sparse vector probably has fewer elements
}
vector_view_raw_index :: inline (v: VectorView, i: int) -> int {
    return i;
}
vector_view_get :: inline (v: VectorView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
vector_view_set :: inline (v: *VectorView, i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
#scope_export



// ##     ## ########    ###    ########     ##     ## #### ######## ##      ##
// ##     ## ##         ## ##   ##     ##    ##     ##  ##  ##       ##  ##  ##
// ##     ## ##        ##   ##  ##     ##    ##     ##  ##  ##       ##  ##  ##
// ######### ######   ##     ## ########     ##     ##  ##  ######   ##  ##  ##
// ##     ## ##       ######### ##            ##   ##   ##  ##       ##  ##  ##
// ##     ## ##       ##     ## ##             ## ##    ##  ##       ##  ##  ##
// ##     ## ######## ##     ## ##              ###    #### ########  ###  ###

HeapVectorView :: struct(V: Type) {
    using _v: VectorType(
        #this,
        V.data_type,
        heap_vector_view_dim,
        VECTOR_FLAGS.VIEW & V.flags,

        heap_vector_view_init,

        heap_vector_view_get,
        heap_vector_view_set,

        heap_vector_view_n_raw,
        heap_vector_view_get, // same as above for dense vectors
        heap_vector_view_set, // same as above for dense vectors
        heap_vector_view_raw_index
    );

    dims: int;
    start_index: int = 0;
    stride: int = 1;
    vec: *V;
}


hvecview :: inline(vec: *$V/VectorType, num_elements: int, start_index: int = 0, stride: int = 1) -> HeapVectorView(V) {
    #if CHECKS {
        assert(start_index + (num_elements-1) * stride <= dim(vec));
    }
    res : HeapVectorView(V);
    res.dim = num_elements;
    res.start_index = start_index;
    res.stride = stride;
    res.vec = vec;
    return res;
}


#scope_file
heap_vector_view_dim :: inline (v: HeapVectorView) -> int {
    return v.dims;
}
heap_vector_view_n_raw :: inline (v: HeapVectorView) -> int {
    return v.dims;
}
heap_vector_view_raw_index :: inline (v: HeapVectorView, i: int) -> int {
    return i;
}
heap_vector_view_get :: inline (v: HeapVectorView, i: int) -> v.data_type {
    return get(v.vec, v.start_index + v.stride * i);
}
heap_vector_view_set :: inline (v: *HeapVectorView($T), i: int, value: v.data_type) {
    set(v.vec, v.start_index + v.stride * i, value);
}
#scope_export




// ########  ######## ##    ##  ######  ########
// ##     ## ##       ###   ## ##    ## ##
// ##     ## ##       ####  ## ##       ##
// ##     ## ######   ## ## ##  ######  ######
// ##     ## ##       ##  ####       ## ##
// ##     ## ##       ##   ### ##    ## ##
// ########  ######## ##    ##  ######  ########

DenseVector :: struct(_T: Type, _N: int) {
    using _v: VectorType(
        #this,
        _T,
        dense_vector_dim,
        VECTOR_FLAGS.DENSE,

        dense_vector_get,
        dense_vector_set,

        dense_vector_n_raw,
        dense_vector_get, // same as above for dense vectors
        dense_vector_set, // same as above for dense vectors
        dense_vector_raw_index
    );

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

#scope_file
dense_vector_dim :: inline (v: DenseVector($T, $N)) -> int {
    return N;
}
dense_vector_n_raw :: inline (v: DenseVector($T, $N)) -> int {
    return N;
}
dense_vector_raw_index :: inline (v: DenseVector($T, $N), i: int) -> int {
    return i;
}
dense_vector_get :: inline (v: DenseVector($T, $N), i: int) -> T {
    return v.data[i];
}
dense_vector_set :: inline (v: *DenseVector($T, $N), i: int, value: T) {
    v.data[i] = value;
}
#scope_export




// ########  ######## ##    ##  ######  ########    ##     ## ########    ###    ########
// ##     ## ##       ###   ## ##    ## ##          ##     ## ##         ## ##   ##     ##
// ##     ## ##       ####  ## ##       ##          ##     ## ##        ##   ##  ##     ##
// ##     ## ######   ## ## ##  ######  ######      ######### ######   ##     ## ########
// ##     ## ##       ##  ####       ## ##          ##     ## ##       ######### ##
// ##     ## ##       ##   ### ##    ## ##          ##     ## ##       ##     ## ##
// ########  ######## ##    ##  ######  ########    ##     ## ######## ##     ## ##


DenseHeapVector :: struct(_T: Type) {
    using _v: VectorType(
        #this,
        _T,
        dense_heap_vector_dim,
        VECTOR_FLAGS.DENSE & VECTOR_FLAGS.HEAP,

        dense_heap_vector_get,
        dense_heap_vector_set,

        dense_heap_vector_n_raw,
        dense_heap_vector_get, // same as above for dense vectors
        dense_heap_vector_set, // same as above for dense vectors
        dense_heap_vector_raw_index
    );

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

#scope_file
dense_heap_vector_dim :: inline (v: DenseHeapVector($T)) -> int {
    return v.dims;
}
dense_heap_vector_n_raw :: inline (v: DenseHeapVector($T)) -> int {
    return v.dims;
}
dense_heap_vector_raw_index :: inline (v: DenseHeapVector($T), i: int) -> int {
    return i;
}
dense_heap_vector_get :: inline (v: DenseHeapVector($T), i: int) -> T {
    return v.data[i];
}
dense_heap_vector_set :: inline (v: *DenseHeapVector($T), i: int, value: T) {
    v.data[i] = value;
}
#scope_export