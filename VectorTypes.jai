
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
    flags: VECTOR_FLAGS,
    get: (m: vector_type, int) -> data_type,
    set: (m: *vector_type, int, data_type) -> ()
){}

get :: inline (v: $V/VectorType, i: int) -> V.data_type {
    return v.get(m, i);
}

set :: inline (v: *$V/VectorType, i: int, value: V.data_type) {
    v.set(v, i, value);
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
        .DENSE,
        dense_vector_get,
        dense_vector_set
    );

    data: [_N]_T;
}

#scope_file
dense_vector_get :: inline (v: DenseVector($T, $N), i: int) -> T {
    return v.data[i];
}
dense_vector_set :: inline (v: *DenseVector($T, $N), i: int, value: T) {
    v.data[i] = value;
}
#scope_export