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
    flags: FLAGS,
    get: (m: matrix_type, int, int) -> data_type,
    set: (m: *matrix_type, int, int, data_type) -> ()
){}
get :: inline (m: $M/MatrixTrait, i: int, j: int) -> M.data_type {
    return m.get(m, i, j);
}
set :: inline (m: *$M/MatrixTrait, i: int, j: int, value: M.data_type) {
    m.set(m, i, j, value);
}



// The only data that actually exists in memory during runtime is "data", the ONELY field of the struct.
FullMatrix :: struct(_T: Type, _N: int, _M: int) {
    using _m: MatrixTrait(
        FullMatrix(_T, _N, _M),
        _T,
        _N,
        _M,
        .NONE,
        full_matrix_get,
        full_matrix_set
    ); 
    data: [_N*_M]_T;
}
full_matrix_get :: inline (m: FullMatrix($T, $N, $M), i: int, j: int) -> T {
    return m.data[i*N + j];
}
full_matrix_set :: inline (m: *FullMatrix($T, $N, $M), i: int, j: int, value: T) {
    m.data[i*N + j] = value;
}



// The only data that actually exists in memory during runtime is "data", the ONELY field of the struct.
UpperTriangularQuadMatrix :: struct(_T: Type, _N: int) {
    using _m: MatrixTrait(
        UpperTriangularQuadMatrix(_T, _N),
        _T,
        _N,
        _N,
        .UPPER_TRIANGULAR,
        upper_triangular_quad_matrix_get,
        upper_triangular_quad_matrix_set
    ); 
    data: [_N*(_N+1)/2]_T;
}
upper_triangular_quad_matrix_get :: inline (m: UpperTriangularQuadMatrix($T, $N), i: int, j: int) -> T {
    return m.data[i]; // ! wrong, TODO correct math here
}
upper_triangular_quad_matrix_set :: inline (m: *UpperTriangularQuadMatrix($T, $N), i: int, j: int, value: T) {
    m.data[i] = value; // ! wrong, TODO correct math here
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
    print("%\n", m);

    do_sth_to_matrix(*m);
    print("%\n", m);

    utm : UpperTriangularQuadMatrix(int, 4);
    print("%\n", utm);
    do_sth_to_matrix(*utm);
    print("%\n", utm);


    res  := create_some_matrix(m, utm);
    res2 := create_some_matrix(utm, utm);
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