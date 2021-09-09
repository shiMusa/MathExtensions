
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