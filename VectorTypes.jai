
// ######## ##          ###     ######    ######
// ##       ##         ## ##   ##    ##  ##    ##
// ##       ##        ##   ##  ##        ##
// ######   ##       ##     ## ##   ####  ######
// ##       ##       ######### ##    ##        ##
// ##       ##       ##     ## ##    ##  ##    ##
// ##       ######## ##     ##  ######    ######

VECTOR_FLAGS :: enum_flags u32 {
    DENSE :: 0;
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
    N_raw: int,
    flags: VECTOR_FLAGS,
    get: (v: vector_type, int) -> data_type,
    set: (v: *vector_type, int, data_type) -> (),
    get_raw: (v: vector_type, int) -> data_type,
    set_raw: (v: *vector_type, int, data_type) -> (),
    from_raw_index: (v: vector_type, int) -> int
){}
get :: inline (v: $V/VectorType, i: int) -> V.data_type {
    return v.get(v, i);
}
set :: inline (v: *$V/VectorType, i: int, value: V.data_type) {
    v.set(v, i, value);
}
get_raw :: inline (v: $V/VectorType, i: int) -> V.data_type {
    return v.get_raw(v, i);
}
set_raw :: inline (v: *$V/VectorType, i: int, value: V.data_type) {
    v.set_raw(v, i, value);
}
for_expansion :: (v: *$V/VectorType, body: Code, flags: For_Flags) #expand {
    i := 0;
    while i < V.N {
        `it_index := i;
        `it := get(<<v, i);
        #insert body;
        i += 1;
    }
}

from_raw_index :: inline (v: $V/VectorType, i: int) -> int {
    return v.from_raw_index(v, i);
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
for_expansion :: (v: *VectorRawForExpansion($V), body: Code, flags: For_Flags) #expand {
    i := 0;
    while i < V.N_raw {
        `it_index := from_raw_index(<<v.vector, i);
        `it := get_raw(<<v.vector, i);
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


DenseVector :: struct(_T: Type, _N: int) {
    using _v: VectorType(
        DenseVector(_T, _N),
        _T,
        _N,
        _N,
        .DENSE,
        dense_vector_get,
        dense_vector_set,
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
    return DenseVector(UpCast(A.data_type, B.data_type).T, N);
}
#scope_export

VectorCastHelper :: struct(V: Type, t: Type) {
    T :: #run vector_cast_helper(V, t);
}
#scope_file
vector_cast_helper :: ($V: Type, $t: Type) -> Type {
    return DenseVector(UpCast(V.data_type,t).T, V.N); // TODO other types
}
#scope_export


VectorSizeHelper :: struct(V: Type, N: int) {
    T :: #run vector_size_helper(V, N);
}
#scope_file
vector_size_helper :: ($V: Type, $N: int) -> Type {
    return DenseVector(V.data_type, N);
}
#scope_export


VectorMatrixHelper :: struct(A: Type, B: Type) {
    T :: #run vector_matrix_helper(A,B);
}
#scope_file
vector_matrix_helper :: ($A: Type, $B: Type) -> Type {
    return DenseMatrix(UpCast(A.data_type, B.data_type), A.N, B.N);
}
#scope_export

VectorVectorHelper :: struct(A: Type, B: Type) {
    T :: #run vector_vector_helper(A, B);
}
#scope_file
vector_vector_helper :: ($A: Type, $B: Type) -> Type {
    return DenseVector(UpCast(A.data_type, B.data_type).T, A.N); // TODO other types
}
#scope_export