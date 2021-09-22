
// ######## ##          ###     ######    ######
// ##       ##         ## ##   ##    ##  ##    ##
// ##       ##        ##   ##  ##        ##
// ######   ##       ##     ## ##   ####  ######
// ##       ##       ######### ##    ##        ##
// ##       ##       ##     ## ##    ##  ##    ##
// ##       ######## ##     ##  ######    ######

VECTOR_FLAGS :: enum_flags u32 {
    DENSE :: 0;
    DENSE_HEAP :: 1;
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
    N: int,
    flags: VECTOR_FLAGS,

    init: (v: *vector_type, $b: bool),

    get: (v: vector_type, int) -> data_type,
    set: (v: *vector_type, int, data_type) -> (),

    n_raw: (v: vector_type) -> int,
    get_raw: (v: vector_type, int) -> data_type,
    set_raw: (v: *vector_type, int, data_type) -> (),
    from_raw_index: (v: vector_type, int) -> int
){}

init :: inline (v: *$V/VectorType, $initialize_variables: bool = true) {
    V.init(v, initialize_variables);
}


get :: inline (v: $V/VectorType, i: int) -> V.data_type {
    return V.get(v, i);
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
    while i < V.N {
        `it_index := i;
        `it := get(<<v, i);
        #insert body;
        i += 1;
    }
}


raw :: inline (v: *$V/VectorType) -> VectorRawForExpansion(V) {
    res : VectorRawForExpansion(V) = ---;
    res.vector = v;
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


// ########  ######## ##    ##  ######  ########
// ##     ## ##       ###   ## ##    ## ##
// ##     ## ##       ####  ## ##       ##
// ##     ## ######   ## ## ##  ######  ######
// ##     ## ##       ##  ####       ## ##
// ##     ## ##       ##   ### ##    ## ##
// ########  ######## ##    ##  ######  ########

DenseVector :: struct(_T: Type, _N: int) {
    using _v: VectorType(
        DenseVector(_T, _N),
        _T,
        _N,
        .DENSE,

        dense_vector_init,

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
        #run assert(is_minor_number(T, U));
        assert(values.count <= N);
    }
    res : DenseVector(T,N);
    for values {
        set(*res, it_index, ncast(T, it));
    }
    return res;
}

#scope_file
dense_vector_init :: inline (v: *DenseVector($T, $N), $initialize_variables: bool = true) {
    #if initialize_variables {
        for 0..N-1 {
            v.data[it] = zero(T);
        }
    }
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


DenseHeapVector :: struct(_T: Type, _N: int) {
    using _v: VectorType(
        DenseHeapVector(_T, _N),
        _T,
        _N,
        .DENSE,

        dense_heap_vector_init,

        dense_heap_vector_get,
        dense_heap_vector_set,

        dense_heap_vector_n_raw,
        dense_heap_vector_get, // same as above for dense vectors
        dense_heap_vector_set, // same as above for dense vectors
        dense_heap_vector_raw_index
    );

    data: []_T;
}

free :: inline (v: $V/DenseHeapVector) {
    array_free(v.data);
}

dhvec :: inline($N: int, values: ..$T) -> DenseHeapVector(T,N) {
    #if CHECKS {
        assert(values.count <= N);
    }
    res := make(DenseHeapVector(T,N), false);
    vec_init(*res, ..values);
    return res;
}

dhvec :: ($T: Type, $N: int, values: ..$U) -> DenseHeapVector(T,N) {
    #if CHECKS {
        #run assert(is_minor_number(T, U));
        assert(values.count <= N);
    }
    res := make(DenseHeapVector(T,N));
    for values {
        set(*res, it_index, ncast(T, it));
    }
    return res;
}

#scope_file
dense_heap_vector_init :: inline (v: *DenseHeapVector($T, $N), $initialize_variables: bool = true) {
    v.data = NewArray(N, T, false);
    #if initialize_variables {
        for 0..N-1 {
            v.data[it] = zero(T);
        }
    }
}
dense_heap_vector_n_raw :: inline (v: DenseHeapVector($T, $N)) -> int {
    return N;
}
dense_heap_vector_raw_index :: inline (v: DenseHeapVector($T, $N), i: int) -> int {
    return i;
}
dense_heap_vector_get :: inline (v: DenseHeapVector($T, $N), i: int) -> T {
    return v.data[i];
}
dense_heap_vector_set :: inline (v: *DenseHeapVector($T, $N), i: int, value: T) {
    v.data[i] = value;
}
#scope_export




// ##     ## ######## ##       ########  ######## ########   ######
// ##     ## ##       ##       ##     ## ##       ##     ## ##    ##
// ##     ## ##       ##       ##     ## ##       ##     ## ##
// ######### ######   ##       ########  ######   ########   ######
// ##     ## ##       ##       ##        ##       ##   ##         ##
// ##     ## ##       ##       ##        ##       ##    ##  ##    ##
// ##     ## ######## ######## ##        ######## ##     ##  ######


VectorHelper :: struct(A: Type, B: Type, N: int) {
    T :: #run vector_helper(A,B,N);
}
#scope_file
vector_helper :: ($A: Type, $B: Type, $N: int) -> Type {
    if A.flags == .DENSE_HEAP || B.flags == .DENSE_HEAP {
        return DenseHeapVector(UpCast(A.data_type, B.data_type).T, N);
    }
    return DenseVector(UpCast(A.data_type, B.data_type).T, N);
}
#scope_export

VectorCastHelper :: struct(V: Type, t: Type) {
    T :: #run vector_cast_helper(V, t);
}
#scope_file
vector_cast_helper :: ($V: Type, $t: Type) -> Type {
    if V.flags == .DENSE_HEAP {
        return DenseHeapVector(UpCast(V.data_type, t).T, V.N);
    }
    return DenseVector(UpCast(V.data_type,t).T, V.N); // TODO other types
}
#scope_export


VectorSizeHelper :: struct(V: Type, N: int) {
    T :: #run vector_size_helper(V, N);
}
#scope_file
vector_size_helper :: ($V: Type, $N: int) -> Type {
    if V.flags == .DENSE_HEAP {
        return DenseHeapVector(V.data_type, N);
    }
    return DenseVector(V.data_type, N);
}
#scope_export


VectorMatrixHelper :: struct(A: Type, B: Type) {
    T :: #run vector_matrix_helper(A,B);
}
#scope_file
vector_matrix_helper :: ($A: Type, $B: Type) -> Type {
    if A.flags == .DENSE_HEAP || B.flags == .DENSE_HEAP {
        return DenseHeapMatrix(UpCast(A.data_type, B.data_type).T, A.N, B.N);
    }
    return DenseMatrix(UpCast(A.data_type, B.data_type).T, A.N, B.N);
}
#scope_export

VectorVectorHelper :: struct(A: Type, B: Type) {
    T :: #run vector_vector_helper(A, B);
}
#scope_file
vector_vector_helper :: ($A: Type, $B: Type) -> Type {
    if A.flags == .DENSE_HEAP || B.flags == .DENSE_HEAP {
        return DenseHeapVector(UpCast(A.data_type, B.data_type).T, A.N);
    }
    return DenseVector(UpCast(A.data_type, B.data_type).T, A.N); // TODO other types
}
#scope_export