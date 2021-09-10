#import "Basic";
// #import "System";


// ##     ##    ###    ######## ########  #### ##     ##
// ###   ###   ## ##      ##    ##     ##  ##   ##   ##
// #### ####  ##   ##     ##    ##     ##  ##    ## ##
// ## ### ## ##     ##    ##    ########   ##     ###
// ##     ## #########    ##    ##   ##    ##    ## ##
// ##     ## ##     ##    ##    ##    ##   ##   ##   ##
// ##     ## ##     ##    ##    ##     ## #### ##     ##


FLAGS :: enum_flags u32 {
    NONE             :: 0;
    LOWER_TRIANGULAR :: (1 << 0);
    UPPER_TRIANGULAR :: (1 << 1);
    SYMMETRIC        :: (1 << 2);
}


// This struct will not exist during runtime in memory since there are no fields.
MatrixTrait :: struct(
    matrix_type: Type,
    data_type: Type,
    N: int,
    M: int,
    N_raw: int,
    flags: FLAGS,
    get: (m: matrix_type, int, int) -> data_type,
    set: (m: *matrix_type, int, int, data_type) -> (),
    get_lin: (m: matrix_type, int) -> data_type,
    set_lin: (m: *matrix_type, int, data_type) -> (),
    from_lin_index: (m: matrix_type, int) -> (int, int),
    to_lin_index: (m: matrix_type, int, int) -> int,
    get_raw: (m: matrix_type, int) -> data_type,
    set_raw: (m: *matrix_type, int, data_type) -> (),
    from_raw_index: (m: matrix_type, int) -> int
){}
from_lin_index :: inline (m: $M/MatrixTrait, i: int) -> int, int {
    r, c := m.from_lin_index(m, i);
    return r, c;
}
to_lin_index :: inline (m: $M/MatrixTrait, i: int, j: int) -> int {
    return m.to_lin_index(m, i,j);
}
get :: inline (m: $M/MatrixTrait, i: int, j: int) -> M.data_type {
    return m.get(m, i, j);
}
set :: inline (m: *$M/MatrixTrait, i: int, j: int, value: M.data_type) {
    m.set(m, i, j, value);
}
get_lin :: inline (m: $M/MatrixTrait, i: int) -> M.data_type {
    return m.get_lin(m, i);
}
set_lin :: inline (m: *$M/MatrixTrait, i: int, value: M.data_type) {
    m.set_lin(m, i, value);
}
for_expansion :: (m: *$M/MatrixTrait, body: Code, flags: For_Flags) #expand {
    i := 0;
    while i < M.N*M.M {
        `it_index := i;
        `it := get_lin(<<m, i);
        `it_row, `it_column := from_lin_index(<<m, i);
        #insert body;
        i += 1;
    }
}

from_raw_index :: inline (m: $M/MatrixTrait, i: int) -> int {
    return m.from_raw_index(m, i);
}
get_raw :: inline (m: $M/MatrixTrait, i: int) -> M.data_type {
    return m.get_raw(m, i);
}
set_raw :: inline (m: *$M/MatrixTrait, i: int, value: M.data_type) {
    m.set_raw(m, i, value);
}
raw :: inline (m: *$M/MatrixTrait) -> MatrixRawForExpansion(M) {
    res : MatrixRawForExpansion(M) = ---;
    res.matrix = m;
    return res;
}
#scope_file
MatrixRawForExpansion :: struct(T: Type) {
    matrix: *T;
}
#scope_export
for_expansion :: (m: *MatrixRawForExpansion($M), body: Code, flags: For_Flags) #expand {
    i := 0;
    while i < M.N_raw {
        `it_index := i;
        `it := get_raw(<<m.matrix, i);
        `it_row, `it_column := from_lin_index(<<m.matrix, from_raw_index(<<m.matrix, i));
        #insert body;
        i += 1;
    }
}



// The only data that actually exists in memory during runtime is "data", the ONELY field of the struct.
FullMatrix :: struct(_T: Type, _N: int, _M: int) {
    using _m: MatrixTrait(
        FullMatrix(_T, _N, _M),
        _T,
        _N,
        _M,
        _N * _M,
        .NONE,
        full_matrix_get,
        full_matrix_set,
        full_matrix_get_lin,
        full_matrix_set_lin,
        full_matrix_from_lin_index,
        full_matrix_to_lin_index,
        full_matrix_get_lin, // raw versions same as normal
        full_matrix_set_lin, // raw versions same as normal
        full_matrix_raw
    ); 
    data: [_N*_M]_T;
}
#scope_file
full_matrix_get :: inline (m: FullMatrix($T, $N, $M), i: int, j: int) -> T {
    return m.data[i*M + j];
}
full_matrix_set :: inline (m: *FullMatrix($T, $N, $M), i: int, j: int, value: T) {
    m.data[i*M + j] = value;
}
full_matrix_from_lin_index :: inline (m: FullMatrix($T, $N, $M), i: int) -> int, int {
    c := i % M;
    r := (i-c)/M;
    return r, c;
}
full_matrix_to_lin_index :: inline (m: FullMatrix($T, $N, $M), i: int, j: int) -> int {
    return i*M+j;
}
full_matrix_get_lin :: inline (m: FullMatrix($T, $N, $M), i: int) -> T {
    return m.data[i];
}
full_matrix_set_lin :: inline (m: *FullMatrix($T, $N, $M), i: int, value: T) {
    m.data[i] = value;
}
full_matrix_raw :: inline (m: FullMatrix($T, $N, $M), i: int) -> int {
    return i;
}
#scope_export




// The only data that actually exists in memory during runtime is "data", the ONELY field of the struct.
UpperTriangularQuadMatrix :: struct(_T: Type, _N: int) {
    using _m: MatrixTrait(
        UpperTriangularQuadMatrix(_T, _N),
        _T,
        _N,
        _N,
        _N*(_N+1)/2,
        .UPPER_TRIANGULAR,
        upper_triangular_quad_matrix_get,
        upper_triangular_quad_matrix_set,
        upper_triangular_quad_matrix_get_lin,
        upper_triangular_quad_matrix_set_lin,
        upper_triangular_quad_matrix_from_lin_index,
        upper_triangular_quad_matrix_to_lin_index,
        upper_triangular_quad_matrix_get_raw,
        upper_triangular_quad_matrix_set_raw,
        upper_triangular_quad_matrix_raw
    ); 
    data: [_N*(_N+1)/2]_T;
}
upper_triangular_quad_matrix_get :: inline (m: UpperTriangularQuadMatrix($T, $N), i: int, j: int) -> T {
    return m.data[i]; // ! wrong, TODO correct math here
}
upper_triangular_quad_matrix_set :: inline (m: *UpperTriangularQuadMatrix($T, $N), i: int, j: int, value: T) {
    m.data[i] = value; // ! wrong, TODO correct math here
}
upper_triangular_quad_matrix_from_lin_index :: (m: UpperTriangularQuadMatrix($T, $N), i: int) -> int, int {
    return i%N, 0; // ! wrong, TODO correct math here
}
upper_triangular_quad_matrix_to_lin_index :: inline (m: UpperTriangularQuadMatrix($T, $N), i: int, j: int) -> int {
    return i; // ! wrong, TODO correct math here
}
upper_triangular_quad_matrix_get_lin :: inline (m: UpperTriangularQuadMatrix($T, $N), i: int) -> T {
    return m.data[i%N];
}
upper_triangular_quad_matrix_set_lin :: inline (m: *UpperTriangularQuadMatrix($T, $N), i: int, value: T) {
    m.data[i%N] = value;
}
// just for testing
upper_triangular_quad_matrix_get_raw :: upper_triangular_quad_matrix_get_lin;
upper_triangular_quad_matrix_set_raw ::upper_triangular_quad_matrix_set_lin;
upper_triangular_quad_matrix_raw :: inline (m: UpperTriangularQuadMatrix($T, $N), i: int) -> int {
    return i;
}


do_sth_to_matrix :: (x: *$M/MatrixTrait) {
    set(x, 1, 1, cast(M.data_type)42);

    y : M;
    print("y = %\n", y);
    print("y.flags = %\n", y.flags);
}






_MatrixHelper :: struct(A: Type, B: Type) {
    T :: #run _matrix_helper(A,B);
}
_matrix_helper :: ($a: Type, $b: Type) -> Type {
    print("a.flags = %\n", a.flags);
    print("b.flags = %\n", b.flags);
    if a.flags == .UPPER_TRIANGULAR && b.flags == .UPPER_TRIANGULAR {
        print("_matrix_helper: returning UpperTriangularQuadMatrix\n");
        return UpperTriangularQuadMatrix(a.data_type, a.N);
    }
    print("_matrix_helper: returning FullMatrix\n");
    return FullMatrix(a.data_type, a.N, a.M);
}

create_some_matrix :: (x: $X/MatrixTrait, y: $Y/MatrixTrait) -> _MatrixHelper(X, Y).T {
    res : _MatrixHelper(X, Y).T;
    print("create_some_matrix: res = %\n", res);
    return res;
}



test_matrix :: () {

    m : FullMatrix(int, 4, 4);
    print("m = %\n", m);

    do_sth_to_matrix(*m);
    print("m = %\n", m);

    utm : UpperTriangularQuadMatrix(int, 4);
    print("utm = %\n", utm);
    do_sth_to_matrix(*utm);
    print("utm = %\n", utm);


    res  := create_some_matrix(m, utm);
    res2 := create_some_matrix(utm, utm);


    for m {
        print("m: it_index %, it %, it_row %, it_column %\n", it_index, it, it_row, it_column);
    }
    for raw(*m) {
        print("raw m: it_index %, it %, it_row %, it_column %\n", it_index, it, it_row, it_column);
    }

    for utm {
        print("utm: it_index %, it %, it_row %, it_column %\n", it_index, it, it_row, it_column);
    }
    for raw(*utm) {
        print("raw utm: it_index %, it %, it_row %, it_column %\n", it_index, it, it_row, it_column);
    }
}


// ##     ##    ###    #### ##    ##
// ###   ###   ## ##    ##  ###   ##
// #### ####  ##   ##   ##  ####  ##
// ## ### ## ##     ##  ##  ## ## ##
// ##     ## #########  ##  ##  ####
// ##     ## ##     ##  ##  ##   ###
// ##     ## ##     ## #### ##    ##

main :: () {
    test_matrix();
}