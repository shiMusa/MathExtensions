
/*

 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##

*/

pstr :: inline (m: $V/MatrixType) -> string {
    return str(m, print_options.offset, print_options.shift, "│", print_options.indent_charactoer_every);
}

str :: (
    m: $M/MatrixType, 
    offset: int = 0, 
    indent: int = 4, 
    indent_character: string = " ", 
    indent_charactoer_every: int = 4
) -> string {
    
    shift :: (ind: int = 0, additional_indent: int = 0) #expand {
        for 0..offset-1 {
            append(*builder, " ");
        }
        for 0..ind-1+additional_indent {
            if it % indent_charactoer_every == 0 {
                append(*builder, indent_character);
            } else {
                append(*builder, " ");
            }
        }
    }
    
    builder: String_Builder;
    defer free_buffers(*builder);
    init_string_builder(*builder);

    append(*builder, "\n");

    column_width := NewArray(cols(m), int);
    defer array_free(column_width);

    for m {
        column_width[it_column] = max(
            column_width[it_column], 
            str(it).count
        );
    }

    matrix_width := 0;
    for column_width {
        matrix_width += it + 2;
    }
    matrix_width -= 2;


    for i: 0..rows(m)-1 {
        shift(indent+4);

        if i == 0 {
            append(*builder, "⎛ ");
        } else if i == rows(m)-1 {
            append(*builder, "⎝ ");
        } else {
            append(*builder, "⎜ ");
        }

        for j: 0..cols(m)-1 {
            padding := column_width[j] - str(get(m,i,j)).count;
            for 0..padding-1 {
                append(*builder, " ");
            }
            append(*builder, str(get(m,i,j)) );
            if j != cols(m)-1 {
                append(*builder, ", ");
            }
        }

        if i == 0 {
            append(*builder, " ⎞");
        } else if i == rows(m)-1 {
            append(*builder, " ⎠");
        } else {
            append(*builder, " ⎟");
        }

        append(*builder, "\n");
    }

    shift(indent);
    append(*builder, "└╴");
    return builder_to_string(*builder);
}


//  ######   #######  ########  ##    ##
// ##    ## ##     ## ##     ##  ##  ##
// ##       ##     ## ##     ##   ####
// ##       ##     ## ########     ##
// ##       ##     ## ##           ##
// ##    ## ##     ## ##           ##
//  ######   #######  ##           ##

copy :: inline (to: *$M/MatrixType, from: M) {
    for raw(from) {
        set_lin(to, it_index, it);
    }
}

copy :: (m: DenseMatrix($T, $N, $M)) -> type_of(m) {
    res : type_of(m) = ---;
    inline copy(*res, m);
    return res;
}
copy :: (m: DenseHeapMatrix) -> type_of(m) {
    res := dhmat(m.data_type,m.rows,m.cols, false); // ! alloc
    inline copy(*res, m);
    return res;
}

//  ######     ###     ######  ########         ########  #######
// ##    ##   ## ##   ##    ##    ##               ##    ##     ##
// ##        ##   ##  ##          ##               ##    ##     ##
// ##       ##     ##  ######     ##               ##    ##     ##
// ##       #########       ##    ##               ##    ##     ##
// ##    ## ##     ## ##    ##    ##               ##    ##     ##
//  ######  ##     ##  ######     ##    #######    ##     #######

cast_to :: inline (to: *$M/MatrixType, from: $K/MatrixType) {
    for raw(from) {
        set_lin(to, it_index, ncast(M.data_type, it));
    }
}

// #### ##    ## #### ########
//  ##  ###   ##  ##     ##
//  ##  ####  ##  ##     ##
//  ##  ## ## ##  ##     ##
//  ##  ##  ####  ##     ##
//  ##  ##   ###  ##     ##
// #### ##    ## ####    ##

mat_init :: (m: *$M/MatrixType, vals: ..$T) {
    inline mat_init_default(m, ..vals);
}

#scope_file
mat_init_default :: (m: *$M/MatrixType, vals: ..$T) {
    #if CHECKS {
        #run assert(is_minor_number_type(M.data_type, T));
        assert(vals.count <= n(m), "too many values (%) for matrix of % entries\n", vals.count, n(m));
    }
    for vals {
        set_lin(m, it_index, ncast(M.data_type, it));
    }
}
#scope_export




mat_init_all :: (m: *$M/MatrixType, value: $T) {
    inline mat_init_all_default(m, value);
}

#scope_file
mat_init_all_default :: (m: *$M/MatrixType, value: $T) {
    v := ncast(M.data_type, value);
    for raw(m) {
        set_raw(m, it_raw_index, v);
    }
}
#scope_export



mat_init_unit :: (m: *$M/MatrixType, $init_zeroes: bool = true) {
    inline mat_init_unit_default(m, init_zeroes);
}

#scope_file
mat_init_unit_default :: (m: *$M/MatrixType, $init_zeroes: bool = true) {
    #if init_zeroes {
        z := zero(M.data_type);
        o := one(M.data_type);
        for m {
            if it_row == it_column {
                set_lin(m, it_index, o);
            } else {
                set_lin(m, it_index, z);
            }
        }
    } else {
        o := one(M.data_type);
        for 0..min(rows(m), cols(m))-1 {
            set(m, it, it, o);
        }
    }
}
#scope_export



mat_init_diag :: (m: *$M/MatrixType, values: ..$T, $init_zeroes: bool = true) {
    inline mat_init_diag_default(m,..values,init_zeroes);
}

#scope_file
mat_init_diag :: (m: *$M/MatrixType, values: ..$T, $init_zeroes: bool = true) {
    #if CHECKS {
        assert(values.count == min(rows(m), cols(m)), "not enough values (% values) for diagonal (% values)\n", values.count, min(rows(m), cols(m)));
    }
    #if init_zeroes {
        z := zero(M.data_type);
        for m {
            if it_column == it_row {
                set_lin(m, it_index, ncast(M.data_type, values[it_row]));
            } else {
                set_lin(m, it_index, z);
            }
        }
    } else {
        for 0..min(rows(m), cols(m))-1 {
            set(m, it, it, ncast(M.data_type, values[it]));
        }
    }
}
#scope_export




mat_init_hadamard_2 :: (m: *$M/MatrixType) {
    inline mat_init_hadamard_2_default(m);
}

#scope_file
mat_init_hadamard_2_default :: inline (m: *$M/MatrixType) {
    #if CHECKS {
        assert(rows(m) == 2 && cols(m) == 2);
    }
    mat_init(m, 1, 1, 1, -1);
}
#scope_export





mat_init_hadamard :: (m: *$M/MatrixType, nestings: int = 0) {
    inline mat_init_hadamard_default(m, nestings);
}

#scope_file
mat_init_hadamard_default :: (m: *$M/MatrixType, nestings: int = 0) {
    #if CHECKS {
        assert(rows(m) == cols(m), "matrix not square (%,%)\n", rows(m), cols(m));
        assert(rows(m) == 1<<nestings, "matrix not of correct size: needed %, given %\n", 1<<nestings, rows(m));
        assert(nestings >= 0, "Nestings must be >=0 but given %\n", nestings);
        assert(rows(m)%2 == 0, "matrix dimensions not even! %\n", rows(m));
    }
    #if nestings == 0 {
        inline mat_init_hadamard_2(m);
    } else {
        r2 := rows(m)/2;
        c2 := cols(m)/2;

        tl := mat_heap_view(m, r2, c2,  0,  0);
        tr := mat_heap_view(m, r2, c2,  0, c2);
        bl := mat_heap_view(m, r2, c2, r2,  0);
        br := mat_heap_view(m, r2, c2, r2, c2);

        mat_init_hadamard(*tl, nestings-1);
        mat_init_hadamard(*tr, nestings-1);
        mat_init_hadamard(*bl, nestings-1);
        mat_init_hadamard(*br, nestings-1);
        neg(*br);
    }
}
#scope_export


mat_init_reflector :: (m: *$M/MatrixType, normal: $V/VectorType) {
    inline mat_init_reflector_default(m,normal);
}

#scope_file
mat_init_reflector_default :: (m: *$M/MatrixType, normal: $V/VectorType) {
    #if CHECKS {
        #run assert(M.data_type == V.data_type, "data_types mismatch %, %\n", M.data_type, V.data_type);
        assert(rows(m) == cols(m), "matrix not square, (%, %)\n", rows(m), cols(m));
        assert(rows(m) == dim(normal), "dimensions mismatch mat(%,%), vec(%)\n", rows(m), cols(m), dim(normal));
    }
    T :: M.data_type;
    dot := zero(T);
    for normal {
        dot += conjugate(it) * it;
    }
    dot = ncast(T, 2)/dot;

    mat_init_unit(m);
    for m {
        // TODO this should be a gaxpy, right?
        // ! added the `conjugate` since it should be |normal><normal|
        set(m, it_row, it_column, get(m, it_row, it_column) - tau * get(normal, it_row) * conjugate(get(normal, it_column)) );
    }
}
#scope_export






//  ######  ##     ## ########  ##     ##    ###    ######## ########  #### ##     ##
// ##    ## ##     ## ##     ## ###   ###   ## ##      ##    ##     ##  ##   ##   ##
// ##       ##     ## ##     ## #### ####  ##   ##     ##    ##     ##  ##    ## ##
//  ######  ##     ## ########  ## ### ## ##     ##    ##    ########   ##     ###
//       ## ##     ## ##     ## ##     ## #########    ##    ##   ##    ##    ## ##
// ##    ## ##     ## ##     ## ##     ## ##     ##    ##    ##    ##   ##   ##   ##
//  ######   #######  ########  ##     ## ##     ##    ##    ##     ## #### ##     ##

// inclusive first bound, exclusive last bound
mat_submatrix :: (
    sm: *$SM/MatrixType,
    m: $M/MatrixType, 
    row_from: int, 
    row_to: int, 
    column_from: int, 
    column_to: int
) {
    inline mat_submatrix_default(sm, m, row_from, row_to, column_from, column_to);
}

#scope_file
mat_submatrix_default :: (
    sm: *$SM/MatrixType,
    m: $M/MatrixType, 
    row_from: int, 
    row_to: int, 
    column_from: int, 
    column_to: int
) {
    #if CHECKS {
        #run assert(SM.data_type == M.data_type, "data_types do not match: %, %\n", SM.data_type, M.data_type);
        assert(rows(sm) == row_to - row_from, "submatrix rows-dimensions mismatch: %, %->%\n", rows(sm), row_from, row_to);
        assert(cols(sm) == column_to - column_from, "submatrix cols-dimensions mismatch: %, %->%\n", cols(sm), column_from, column_to);
        assert(row_to > row_from, "cannot create submatrix with rows %->%\n", row_to, row_from);
        assert(column_to > column_from, "cannot create submatrix with columns %->%\n", column_to, column_from);
        assert(row_to <= rows(m));
        assert(column_to <= cols(m));
        assert(row_from >= 0);
        assert(column_from >= 0);
    }

    for i: 0..row_to-row_from-1 {
        for j: 0..column_to-column_from-1 {
            set(sm, i, j, get(m, i+row_from, j+column_from));
        }
    }
}
#scope_export








/*
########   #######  ##      ##          ######   #######  ##       ##     ## ##     ## ##    ##
##     ## ##     ## ##  ##  ##         ##    ## ##     ## ##       ##     ## ###   ### ###   ##
##     ## ##     ## ##  ##  ##         ##       ##     ## ##       ##     ## #### #### ####  ##
########  ##     ## ##  ##  ##         ##       ##     ## ##       ##     ## ## ### ## ## ## ##
##   ##   ##     ## ##  ##  ##         ##       ##     ## ##       ##     ## ##     ## ##  ####
##    ##  ##     ## ##  ##  ##         ##    ## ##     ## ##       ##     ## ##     ## ##   ###
##     ##  #######   ###  ###           ######   #######  ########  #######  ##     ## ##    ##
*/


row :: (v: *$V/VectorType, m: $M/MatrixType, i: int) {
    inline row_default(v,m,i);
}

#scope_file
row_default :: (v: *$V/VectorType, m: $M/MatrixType, i: int) {
    #if CHECKS {
        #run assert(V.data_type == M.data_type, "data_types do not match: %, %\n", V.data_type, M.data_type);
        assert(dim(v) == cols(m), "dimensions mismatch v%, m%\n", dim(v), cols(m));
        assert(i >= 0 && i < rows(m));
    }
    for v {
        set(v, it_index, get(m,i,it_index));
    }
}
#scope_export

column :: (v: *$V/VectorType, m: $M/MatrixType, i: int) {
    inline column_default(v,m,i);
}

#scope_file
column_default :: (v: *$V/VectorType, m: $M/MatrixType, i: int) {
    #if CHECKS {
        #run assert(V.data_type == M.data_type, "data_types do not match: %, %\n", V.data_type, M.data_type);
        assert(dim(v) == rows(m), "dimensions mismatch v%, m%\n", dim(v), rows(m));
        assert(i >= 0 && i < cols(m));
    }
    for v {
        set(v, it_index, get(m,it_index,i));
    }
}
#scope_export




/*
 #######  ########  ######## ########     ###    ########  #######  ########   ######
##     ## ##     ## ##       ##     ##   ## ##      ##    ##     ## ##     ## ##    ##
##     ## ##     ## ##       ##     ##  ##   ##     ##    ##     ## ##     ## ##
##     ## ########  ######   ########  ##     ##    ##    ##     ## ########   ######
##     ## ##        ##       ##   ##   #########    ##    ##     ## ##   ##         ##
##     ## ##        ##       ##    ##  ##     ##    ##    ##     ## ##    ##  ##    ##
 #######  ##        ######## ##     ## ##     ##    ##     #######  ##     ##  ######
*/




operator == :: (a: $A/MatrixType, b: $B/MatrixType) -> bool {
    return inline equals(a,b);
}

equals :: (a: $A/MatrixType, b: $B/MatrixType) -> bool {
    #if A.flags == .DENSE && B.flags == .DENSE {
        return equals_dense_dense(a,b);
    } else {
        // * default implementation
        return equals_default(a,b);
    }
}

#scope_file
equals_default :: (a: $A/MatrixType, b: $B/MatrixType) -> bool {
    if rows(a) != rows(b) || cols(a) != cols(b) then return false;

    for a {
        if it != get_lin(b, it_index) then return false;
    }

    return true;
}
equals_dense_dense :: (a: DenseMatrix($AT, $AN, $AM), b: DenseMatrix($BT, $BN, $BM)) -> bool {
    #if AN != BN || AM != BM {
        return false;
    } else {
        for a {
            if it != get_lin(b, it_index) then return false;
        }
        return true;
    }
}
#scope_export

// ##     ##    ###    ########                     ##     ##    ###    ########
// ###   ###   ## ##      ##            ##   ##     ###   ###   ## ##      ##
// #### ####  ##   ##     ##             ## ##      #### ####  ##   ##     ##
// ## ### ## ##     ##    ##           #########    ## ### ## ##     ##    ##
// ##     ## #########    ##             ## ##      ##     ## #########    ##
// ##     ## ##     ##    ##       ###  ##   ##     ##     ## ##     ##    ##
// ##     ## ##     ##    ##       ###              ##     ## ##     ##    ##


mul_el :: (a: DenseMatrix($A, $N, $M), b: DenseMatrix($B, N, M)) -> DenseMatrix(UpCast(A,B).T, N, M) {
    res : DenseMatrix(UpCast(A,B).T, N, M);
    inline mul_el(*res,a,b);
    return res;
}
mul_el :: (a: DenseHeapMatrix($A), b: DenseHeapMatrix($B)) -> DenseHeapMatrix(UpCast(A,B).T) {
    res := dhmat(UpCast(A,B).T, rows(a), cols(a)); // ! alloc
    inline mul_el(*res,a,b);
    return res;
}


mul_el :: (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    inline mul_el_default(c,a,b);
}

mul_el :: (a: *$A/MatrixType, b: $B/MatrixType) {
    inline mul_el_default(a, <<a, b);
}

#scope_file
mul_el_default :: (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B.data_type));
        assert(rows(a) == rows(c) && cols(a) == cols(c));
        assert(rows(b) == rows(c) && cols(b) == cols(c));
    }
    for c {
        set_lin(c, it_index, it +  get_lin(a, it_index) * get_lin(b, it_index));
    }
}
#scope_export

// ##     ##    ###    ########              ##    ##     ##    ###    ########
// ###   ###   ## ##      ##                ##     ###   ###   ## ##      ##
// #### ####  ##   ##     ##               ##      #### ####  ##   ##     ##
// ## ### ## ##     ##    ##              ##       ## ### ## ##     ##    ##
// ##     ## #########    ##             ##        ##     ## #########    ##
// ##     ## ##     ##    ##       ###  ##         ##     ## ##     ##    ##
// ##     ## ##     ##    ##       ### ##          ##     ## ##     ##    ##


div_el :: (a: DenseMatrix($A, $N, $M), b: DenseMatrix($B, N, M)) -> DenseMatrix(UpCast(A,B).T, N, M) {
    res : DenseMatrix(UpCast(A,B).T, N, M);
    inline div_el(*res,a,b);
    return res;
}
div_el :: (a: DenseHeapMatrix($A), b: DenseHeapMatrix($B)) -> DenseHeapMatrix(UpCast(A,B).T) {
    res := dhmat(UpCast(A,B).T, rows(a), cols(a)); // ! alloc
    inline div_el(*res,a,b);
    return res;
}


div_el :: (a: $A/MatrixType, b: $B/MatrixType, c: *$C/MatrixType) {
    inline div_el_default(c,a,b);
}

div_el :: (a: *$A/MatrixType, b: $B/MatrixType) {
    inline div_el_default(a,<<a,b);
}


#scope_file
div_el_default :: (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B.data_type));
        assert(rows(a) == rows(c) && cols(a) == cols(c));
        assert(rows(b) == rows(c) && cols(b) == cols(c));
    }
    for c {
        set_lin(c, it_index, it + get_lin(a, it_index) / get_lin(b, it_index));
    }
}
#scope_export

// ##     ##    ###    ########                 ##     ##    ###    ########
// ###   ###   ## ##      ##        ##   ##     ###   ###   ## ##      ##
// #### ####  ##   ##     ##         ## ##      #### ####  ##   ##     ##
// ## ### ## ##     ##    ##       #########    ## ### ## ##     ##    ##
// ##     ## #########    ##         ## ##      ##     ## #########    ##
// ##     ## ##     ##    ##        ##   ##     ##     ## ##     ##    ##
// ##     ## ##     ##    ##                    ##     ## ##     ##    ##

operator * :: (a: DenseMatrix($A, $AN, $M), b: DenseMatrix($B, M, $BM)) -> DenseMatrix(UpCast(A,B).T, AN, BM) {
    res : DenseMatrix(UpCast(A,B).T, AN, BM);
    inline mul(*res,a,b);
    return res;
}
operator * :: (a: DenseHeapMatrix($A), b: DenseHeapMatrix($B)) -> DenseHeapMatrix(UpCast(A,B).T) {
    res := dhmat(UpCast(A,B).T, rows(a), cols(b)); // ! alloc
    inline mul(*res,a,b);
    return res;
}

// TODO implement operator overloads for other `MatrixType`s where it makes sense


// * Assuming res is initialized to zeroes.
// * User has to know, what type the result has to be!
mul :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if A.flags & MATRIX_FLAGS.TRIANGULAR && 
        B.flags & MATRIX_FLAGS.TRIANGULAR {
        inline mat_mul_mat_triangular(res, a, b);
    } else {
        inline mat_mul_mat_default(res,a,b);
    }
}

#scope_file
// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 2.11 Case Study: Matrix-Matrix Multiplication, page 118
// ? According to Table 2.12, the best performing loop order is j/k/i with
// ? the least number of cache misses.
// ? 
// ? Similar discussion in:
// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.5 ijk Matrix Multiplication, page 8

// ? saxpy ----------------------------------------------------------------------------
// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.6 Dot Product Matrix Multiplication, page 10
// ? Algorithm 1.1.7 Saxpy Matrix Multiplication, page 10
// ? C = C + A * B
mat_mul_mat_default :: (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type) && is_minor_number_type(C.data_type, B.data_type));
        assert(cols(a) == rows(b)  , "Matrix dimensions don't match: a.cols = %, b.rows = %", cols(a), rows(b)  );
        assert(rows(a) == rows(c), "Matrix dimensions don't match: a.rows = %, c.rows = %", rows(a), rows(c));
        assert(cols(b) == cols(c), "Matrix dimensions don't match: b.cols = %, c.cols = %", cols(b), cols(c));
    }

    for j: 0..cols(b)-1 {
        for k: 0..cols(a)-1 {
            tmp := get(b,k,j);
            for i: 0..rows(c)-1 {
                set(c,i,j, get(c,i,j) + get(a,i,k) * tmp);
            }
        }
    }
}

// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.2.1 Triangular Matrix Multiplication, page 16
mat_mul_mat_triangular :: (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type) && is_minor_number_type(C.data_type, B.data_type));
        assert(cols(a) == rows(b)  , "Matrix dimensions don't match: a.cols = %, b.rows = %", cols(a), rows(b)  );
        assert(rows(a) == rows(c), "Matrix dimensions don't match: a.rows = %, c.rows = %", rows(a), rows(c));
        assert(cols(b) == cols(c), "Matrix dimensions don't match: b.cols = %, c.cols = %", cols(b), cols(c));
    }

    #if C.flags & MATRIX_FLAGS.TRIANGULAR {
        if !a.transposed && !b.transposed {
            #if CHECKS {
                assert(!c.transposed);
            }
            for i: 0..rows(a)-1 {
                for j: i..cols(b)-1 {
                    for k: i..min(j, rows(b)-1) {
                        set(c,i,j, get(c,i,j, false) + get(a,i,k, false) * get(b,k,j, false), false);
                    }
                }
            }
        } else if a.transposed && b.transposed {
            #if CHECKS {
                assert(c.transposed);
            }
            // ? just use C = A * B  =>  Cᵀ = Bᵀ * Aᵀ, so switch a and b, 
            //? pretend they're not transposed, and disable all transposed related features (check_transposition = false).
            for i: 0..rows(b, false)-1 {
                for j: i..cols(a, false)-1 {
                    for k: i..min(j, rows(a, false)-1) {
                        set(c,i,j, get(c,i,j, false) + get(b,i,k, false) * get(a,k,j, false), false);
                    }
                }
            }
        }
    } else {
        if !a.transposed && !b.transposed {
            for i: 0..rows(a, false)-1 {
                for j: i..cols(b, false)-1 {
                    for k: i..min(j, rows(b, false)-1) {
                        set(c,i,j, get(c,i,j) + get(a,i,k, false) * get(b,k,j, false));
                    }
                }
            }
        } else if a.transposed && b.transposed {
            // ? just use C = A * B  =>  Cᵀ = Bᵀ * Aᵀ, so switch a and b, 
            //? pretend they're not transposed, and disable all transposed related features (check_transposition = false).
            for i: 0..rows(b, false)-1 {
                for j: i..cols(a, false)-1 {
                    for k: i..min(j, rows(a, false)-1) {
                        set(c,j,i, get(c,j,i) + get(b,i,k, false) * get(a,k,j, false));
                    }
                }
            }
        } else if a.transposed && !b.transposed {
            n := rows(b, false);
            for i: 0..rows(a)-1 {
                for j: 0..cols(b)-1 {
                    for k: 0..min(i,j,n-1) {
                        set(c,i,j, get(c,i,j) + get(a,k,i,false) * get(b,k,j,false));
                    }
                }
            }
        } else if !a.transposed && b.transposed {
            n := cols(a);
            for i: 0..min(rows(a), n)-1 {
                for j: 0..min(cols(b), n)-1 {
                    for k: max(i,j)..n-1 {
                        set(c,i,j, get(c,i,j) + get(a,i,k,false) * get(b,j,k,false));
                    }
                }
            }
        }
    }
}






test_mat_mul_mat_triangular :: () {
    println_push("mat_mul_mat_triangular", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        a := tmat(int, 5, 3, 1,2,3,4,5,6);
        b := tmat(int, 3, 4, 1,2,3,4,5,6,7,8,9);
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : TriangularMatrix(int, 5, 4);
        mat_mul_mat_triangular(*c, a, b);
        println("c = %", pstr(c));
        res := tmat(int, 5, 4, 1,12,39,45,20,64,73,48,54);
        println("res = %", pstr(res));
        assert(c == res);
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        a := tmat(int, 3, 5, 1,2,3,4,5,6,7,8,9,10,11,12);
        b := tmat(int, 5, 4, 1,2,3,4,5,6,7,8,9,10);
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : TriangularMatrix(int, 3, 4);
        mat_mul_mat_triangular(*c, a, b);
        println("c = %", pstr(c));
        res := tmat(int, 3, 4, 1,12,39,85,30,92,185,80,200);
        println("res = %", pstr(res));
        assert(c == res);
    }

    {
        println_push("test 3: transposed", color = .FG_YELLOW);
        a := tmat(int, 5, 3, 1,2,3,4,5,6);
        a.transposed = true;
        b := tmat(int, 3, 4, 1,2,3,4,5,6,7,8,9);
        b.transposed = true;
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : TriangularMatrix(int, 5, 4);
        c.transposed = true;
        mat_mul_mat_triangular(*c, b, a);
        println("c = %", pstr(c));
        res := tmat(int, 5, 4, 1,12,39,45,20,64,73,48,54);
        res.transposed = true;
        println("res = %", pstr(res));
        assert(c == res);
    }

    {
        println_push("test 4: transposed", color = .FG_YELLOW);
        a := tmat(int, 3, 5, 1,2,3,4,5,6,7,8,9,10,11,12);
        a.transposed = true;
        b := tmat(int, 5, 4, 1,2,3,4,5,6,7,8,9,10);
        b.transposed = true;
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : TriangularMatrix(int, 3, 4);
        c.transposed = true;
        mat_mul_mat_triangular(*c, b, a);
        println("c = %", pstr(c));
        res := tmat(int, 3, 4, 1,12,39,85,30,92,185,80,200);
        res.transposed = true;
        println("res = %", pstr(res));
        assert(c == res);
    }

    

    {
        println_push("test 5: dense result", color = .FG_YELLOW);
        a := tmat(int, 5, 3, 1,2,3,4,5,6);
        b := tmat(int, 3, 4, 1,2,3,4,5,6,7,8,9);
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : DenseMatrix(int, 5, 4);
        mat_mul_mat_triangular(*c, a, b);
        println("c = %", pstr(c));
        res := tmat(int, 5, 4, 1,12,39,45,20,64,73,48,54);
        println("res = %", pstr(res));
        assert(c == res);
    }

    {
        println_push("test 6: dense result", color = .FG_YELLOW);
        a := tmat(int, 3, 5, 1,2,3,4,5,6,7,8,9,10,11,12);
        b := tmat(int, 5, 4, 1,2,3,4,5,6,7,8,9,10);
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : DenseMatrix(int, 3, 4);
        mat_mul_mat_triangular(*c, a, b);
        println("c = %", pstr(c));
        res := tmat(int, 3, 4, 1,12,39,85,30,92,185,80,200);
        println("res = %", pstr(res));
        assert(c == res);
    }

    {
        println_push("test 7: transposed, dense result", color = .FG_YELLOW);
        a := tmat(int, 5, 3, 1,2,3,4,5,6);
        a.transposed = true;
        b := tmat(int, 3, 4, 1,2,3,4,5,6,7,8,9);
        b.transposed = true;
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : DenseMatrix(int, 4, 5);
        mat_mul_mat_triangular(*c, b, a);
        println("c = %", pstr(c));
        res := tmat(int, 5, 4, 1,12,39,45,20,64,73,48,54);
        res.transposed = true;
        println("res = %", pstr(res));
        assert(c == res);
    }

    {
        println_push("test 8: transposed, dense result", color = .FG_YELLOW);
        a := tmat(int, 3, 5, 1,2,3,4,5,6,7,8,9,10,11,12);
        a.transposed = true;
        b := tmat(int, 5, 4, 1,2,3,4,5,6,7,8,9,10);
        b.transposed = true;
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : DenseMatrix(int, 4, 3);
        mat_mul_mat_triangular(*c, b, a);
        println("c = %", pstr(c));
        res := tmat(int, 3, 4, 1,12,39,85,30,92,185,80,200);
        res.transposed = true;
        println("res = %", pstr(res));
        assert(c == res);
    }

    {
        println_push("test 9: opposite diagonals", color = .FG_YELLOW);
        a := tmat(int, 3, 5, 
            1, 2, 3, 4, 5,
               6, 7, 8, 9,
                 10,11,12
        );
        a.transposed = true;
        b := tmat(int, 3, 4, 
            1,2,3,4,
              5,6,7,
                8,9
        );
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : DenseMatrix(int, 5, 4);
        mat_mul_mat_triangular(*c, a, b);
        println("c = %", pstr(c));
        res := dmat(int, 5, 4, 
            1,  2,  3,  4,
            2, 34, 42, 50,
            3, 41,131,151,
            4, 48,148,171,
            5, 55,165,191
        );
        println("res = %", pstr(res));
        assert(c == res);
    }
    
    {
        println_push("test 10: opposite diagonals, tranposed", color = .FG_YELLOW);
        a := tmat(int, 3, 4, 
            1, 2, 3, 4,
               5, 6, 7,
                  8, 9
        );
        b := tmat(int, 5, 4,
            1, 2, 4, 7,
               3, 5, 8,
                  6, 9,
                    10
        );
        b.transposed = true;
        println("a = %", pstr(a));
        println("b = %", pstr(b));
        c : DenseMatrix(int, 3, 5);
        mat_mul_mat_triangular(*c, a, b);
        println("c = %", pstr(c));
        res := dmat(int, 3, 5,
            45, 53, 54, 40, 0,
            83,101, 99, 70, 0,
            95,112,129, 90, 0
        );
        println("res = %", pstr(res));
        assert(c == res);
    }
} 

test_mat_mul_mat :: () {
    println_push("mat_mul_mat_...", color = .FG_WHITE);
    test_mat_mul_mat_triangular();
}
#scope_export

// ##     ##    ###    ########                 ##     ## ########  ######
// ###   ###   ## ##      ##        ##   ##     ##     ## ##       ##    ##
// #### ####  ##   ##     ##         ## ##      ##     ## ##       ##
// ## ### ## ##     ##    ##       #########    ##     ## ######   ##
// ##     ## #########    ##         ## ##       ##   ##  ##       ##
// ##     ## ##     ##    ##        ##   ##       ## ##   ##       ##    ##
// ##     ## ##     ##    ##                       ###    ########  ######



operator * :: (a: DenseMatrix($A, $N, $M), b: DenseVector($B, M)) -> DenseVector(UpCast(A,B).T, N) {
    res : DenseVector(UpCast(A,B).T, N);
    inline mul(*res,a,b);
    return res;
}
operator * :: (a: DenseHeapMatrix($A), b: DenseHeapVector($B)) -> DenseHeapVector(UpCast(A,B).T) {
    res : DenseVector(UpCast(A,B).T, rows(a)); // ! alloc
    inline mul(*res,a,b);
    return res;
}


// * Assuming res is initialized to zeroes.
// * User has to know, what type the result has to be!
mul :: (res: *$C/VectorType, a: $A/MatrixType, b: $B/VectorType) {
    inline mat_mul_vec_default(res,a,b);
}

#scope_file
mat_mul_vec_default :: (res: *$C/VectorType, a: $A/MatrixType, b: $B/VectorType) {
    #if CHECKS {
        assert(cols(a) == dim(b) , "Dimensions mismatch: a.cols = %, b.dim   = %", cols(a), dim(b) );
        assert(rows(a) == dim(res), "Dimensions mismatch: a.rows = %, res.dim = %", rows(a), dim(res));
    }

    for a {
        set(res, it_row, get(res, it_row) + it * get(b, it_column));
    }
}
#scope_export

// ##     ## ########  ######                  ##     ##    ###    ########
// ##     ## ##       ##    ##     ##   ##     ###   ###   ## ##      ##
// ##     ## ##       ##            ## ##      #### ####  ##   ##     ##
// ##     ## ######   ##          #########    ## ### ## ##     ##    ##
//  ##   ##  ##       ##            ## ##      ##     ## #########    ##
//   ## ##   ##       ##    ##     ##   ##     ##     ## ##     ##    ##
//    ###    ########  ######                  ##     ## ##     ##    ##


operator * :: (a: DenseVector($A, $N), b: DenseMatrix($B, N, $M)) -> DenseVector(UpCast(A,B).T, M) {
    res : DenseVector(UpCast(A,B).T, M);
    inline mul(*res,a,b);
    return res;
}
operator * :: (a: DenseHeapVector($A), b: DenseHeapMatrix($B)) -> DenseHeapVector(UpCast(A,B).T) {
    res := dhvec(UpCast(A,B).T, cols(a)); // ! alloc
    inline mul(*res,a,b);
    return res;
}


// * Assuming res is initialized to zeroes.
// * User has to know, what type the result has to be!
mul :: (res: *$C/VectorType, a: $A/VectorType, b: $B/MatrixType) {
    inline vec_mul_mat_default(res,a,b);
}

#scope_file
vec_mul_mat_default :: (res: *$C/VectorType, a: $A/VectorType, b: $B/MatrixType) {
    #if CHECKS {
        assert(dim(a)  == rows(b), "Dimensions mismatch: a.dim  = %, b.rows  = %", dim(a) , rows(b));
        assert(cols(b) == dim(c) , "Dimensions mismatch: b.cols = %, res.dim = %", cols(b), dim(c) );
    }

    for j: 0..rows(b)-1 {
        for i: 0..dim(a)-1 {
            set(res,j, get(res,j) + get(a,i) * get(b,i,j));
        }
    }
}
#scope_export


// ##     ##    ###    ########                  ######   ######     ###    ##
// ###   ###   ## ##      ##        ##   ##     ##    ## ##    ##   ## ##   ##
// #### ####  ##   ##     ##         ## ##      ##       ##        ##   ##  ##
// ## ### ## ##     ##    ##       #########     ######  ##       ##     ## ##
// ##     ## #########    ##         ## ##            ## ##       ######### ##
// ##     ## ##     ##    ##        ##   ##     ##    ## ##    ## ##     ## ##
// ##     ## ##     ##    ##                     ######   ######  ##     ## ########


operator * :: (a: DenseMatrix($T, $N, $M), b: $B) -> DenseMatrix(UpCast(T,B).T, N, M) #modify {
    if !is_number(B) {B = null; return;}
} #symmetric {
    res : DenseMatrix(UpCast(T,B).T, N, M) = ---;
    inline mul(*res,a, b);
    return res;
}
operator * :: (a: DenseHeapMatrix($T), b: $B) -> DenseHeapMatrix(UpCast(T,B).T) #modify {
    if !is_number(B) {B = null; return;}
} #symmetric {
    res := dhmat(UpCast(T,B).T, rows(a), cols(a), false); // ! alloc
    inline mul(*res,a, b);
    return res;
}


mul :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B) #modify {
    if !is_number(B) {B = null; return;}
} {
    #if CHECKS {
        assert(rows(a) == rows(c), "Dimensions mismatch: a.rows  = %, res.rows  = %", rows(a), rows(c));
        assert(cols(a) == cols(c), "Dimensions mismatch: a.cols  = %, res.cols  = %", cols(a), cols(c));
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B));
    }

    for raw(a) {
        set_raw(res, it_raw_index, it * b);
    }
}


mul :: (a: *$A/MatrixType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    #if CHECKS {
        #run assert(is_minor_number_type(A.data_type, B));
    }

    for raw(a) {
        set_raw(a, it_raw_index, it * b);
    }
}



// // ##     ##    ###    ########          ##    ##     ##    ###    ########
// // ###   ###   ## ##      ##            ##     ###   ###   ## ##      ##
// // #### ####  ##   ##     ##           ##      #### ####  ##   ##     ##
// // ## ### ## ##     ##    ##          ##       ## ### ## ##     ##    ##
// // ##     ## #########    ##         ##        ##     ## #########    ##
// // ##     ## ##     ##    ##        ##         ##     ## ##     ##    ##
// // ##     ## ##     ##    ##       ##          ##     ## ##     ##    ##

// // mat_div_mat :: (a: $A/MatrixType, b: $B/MatrixType, res: *$C/MatrixType) {
// // }
// operator / :: inline (a: $A/MatrixType, b: $B/MatrixType) -> MatrixMultiplicationHelper(A,B,A.M,B.M).T { 
//     return inverse(b) * a;
// }

// operator / :: (a: $A/VectorType, b: $B/MatrixType) -> MatrixVectorMultiplicationHelper(B, A, A.N).T { 
//     return inverse(b) * a;
// }


// ##     ##    ###    ########          ##     ######   ######     ###    ##
// ###   ###   ## ##      ##            ##     ##    ## ##    ##   ## ##   ##
// #### ####  ##   ##     ##           ##      ##       ##        ##   ##  ##
// ## ### ## ##     ##    ##          ##        ######  ##       ##     ## ##
// ##     ## #########    ##         ##              ## ##       ######### ##
// ##     ## ##     ##    ##        ##         ##    ## ##    ## ##     ## ##
// ##     ## ##     ##    ##       ##           ######   ######  ##     ## ########

operator / :: (a: DenseMatrix($T, $N, $M), b: $B) -> DenseMatrix(UpCast(T,B).T, N, M) #modify {
    if !is_number(B) {B = null; return;}
} {
    res : DenseMatrix(UpCast(T,B).T, N, M) = ---;
    inline div(*res,a, b);
    return res;
}
operator / :: (a: DenseHeapMatrix($T), b: $B) -> DenseHeapMatrix(UpCast(T,B).T) #modify {
    if !is_number(B) {B = null; return;}
} {
    res := dhmat(UpCast(T,B).T, rows(a), cols(a), false); // ! alloc
    inline div(*res,a, b);
    return res;
}



div :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B) #modify {
    if !is_number(B) { B = null; return; }
} {
    #if CHECKS {
        assert(rows(a) == rows(c), "Dimensions mismatch: a.rows  = %, res.rows  = %", rows(a), rows(c));
        assert(cols(a) == cols(c), "Dimensions mismatch: a.cols  = %, res.cols  = %", cols(a), cols(c));
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B));
    }

    for raw(a) {
        set_raw(res, it_raw_index, it / b);
    }
}

div :: (a: *$A/MatrixType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    #if CHECKS {
        #run assert(is_minor_number_type(A.data_type, B));
    }

    for raw(a) {
        set_raw(a, it_raw_index, it / b);
    }
}


// ##     ##    ###    ########              ##     ##    ###    ########
// ###   ###   ## ##      ##         ##      ###   ###   ## ##      ##
// #### ####  ##   ##     ##         ##      #### ####  ##   ##     ##
// ## ### ## ##     ##    ##       ######    ## ### ## ##     ##    ##
// ##     ## #########    ##         ##      ##     ## #########    ##
// ##     ## ##     ##    ##         ##      ##     ## ##     ##    ##
// ##     ## ##     ##    ##                 ##     ## ##     ##    ##


operator + :: (a: DenseMatrix($A, $N, $M), b: DenseMatrix($B, N, M)) -> DenseMatrix(UpCast(A,B).T, N, M) {
    res : DenseMatrix(UpCast(A,B).T, N, M) = ---;
    inline add(*res,a,b);
    return res;
}
operator + :: (a: DenseHeapMatrix($A), b: DenseHeapMatrix($B)) -> DenseHeapMatrix(UpCast(A,B).T) {
    res := dhmat(UpCast(A,B).T, rows(a), cols(a), false); // ! alloc
    inline add(*res,a,b);
    return res;
}


add :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    inline mat_add_mat_default(res,a,b);
}

#scope_file
mat_add_mat_default :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        assert(rows(a) == rows(b)  , "Matrix dimensions don't match: a.rows = %, b.rows = %", rows(a), rows(b)  );
        assert(rows(a) == rows(res), "Matrix dimensions don't match: a.rows = %, b.rows = %", rows(a), rows(res));
        assert(cols(a) == cols(b)  , "Matrix dimensions don't match: a.cols = %, c.cols = %", cols(a), cols(b)  );
        assert(cols(a) == cols(res), "Matrix dimensions don't match: a.cols = %, c.cols = %", cols(a), cols(res));
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B.data_type));
    }

    for a {
        set_lin(res, it_index, it + get_lin(b, it_index));
    }
}
#scope_export

add :: (a: *$A/MatrixType, b: $B/MatrixType) {
    inline add_default(a,b);
}

#scope_file
add_default :: (a: *$A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(is_minor_number_type(A.data_type, B.data_type), "cannot cast % to %.", B.data_type, A.data_type);
        assert(rows(a)== rows(b) && cols(a) == cols(b), "matrix dimensions do not match A(% x %), B(% x %)", rows(a), cols(a), rows(b), cols(b));
    }

    for a {
        set_lin(a, it_index, it + get_lin(b, it_index));
    }
}
#scope_export


// ##     ##    ###    ########               ##     ##    ###    ########
// ###   ###   ## ##      ##                  ###   ###   ## ##      ##
// #### ####  ##   ##     ##                  #### ####  ##   ##     ##
// ## ### ## ##     ##    ##       #######    ## ### ## ##     ##    ##
// ##     ## #########    ##                  ##     ## #########    ##
// ##     ## ##     ##    ##                  ##     ## ##     ##    ##
// ##     ## ##     ##    ##                  ##     ## ##     ##    ##


operator - :: (a: DenseMatrix($A, $N, $M), b: DenseMatrix($B, N, M)) -> DenseMatrix(UpCast(A,B).T, N, M) {
    res : DenseMatrix(UpCast(A,B).T, N, M) = ---;
    matt_sub_mat(*res,a,b);
    return res;
}
operator + :: (a: DenseHeapMatrix($A), b: DenseHeapMatrix($B)) -> DenseHeapMatrix(UpCast(A,B).T) {
    res := dhmat(UpCast(A,B).T, rows(a), cols(a), false); // ! alloc
    inline sub(*res,a,b);
    return res;
}

sub :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    inline mat_sub_mat_default(res,a,b);
}

#scope_file
mat_sub_mat_default :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        assert(rows(a) == rows(b)  , "Matrix dimensions don't match: a.rows = %, b.rows = %", rows(a), rows(b)  );
        assert(rows(a) == rows(res), "Matrix dimensions don't match: a.rows = %, b.rows = %", rows(a), rows(res));
        assert(cols(a) == cols(b)  , "Matrix dimensions don't match: a.cols = %, c.cols = %", cols(a), cols(b)  );
        assert(cols(a) == cols(res), "Matrix dimensions don't match: a.cols = %, c.cols = %", cols(a), cols(res));
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B.data_type));
    }

    for a {
        set_lin(res, it_index, it - get_lin(b, it_index));
    }
}
#scope_export

sub :: (a: *$A/MatrixType, b: $B/MatrixType) {
    inline sub_default(a,b);
}

#scope_file
sub_default :: (a: *$A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(is_minor_number_type(A.data_type, B.data_type), "cannot cast % to %.", B.data_type, A.data_type);
        assert(rows(a)== rows(b) && cols(a) == cols(b), "matrix dimensions do not match A(% x %), B(% x %)", rows(a), cols(a), rows(b), cols(b));
    }

    for a {
        set_lin(a, it_index, it - get_lin(b, it_index));
    }
}
#scope_export


// ##    ## ########  ######
// ###   ## ##       ##    ##
// ####  ## ##       ##
// ## ## ## ######   ##   ####
// ##  #### ##       ##    ##
// ##   ### ##       ##    ##
// ##    ## ########  ######


operator - :: (a: DenseMatrix($T, $N, $M)) -> DenseMatrix(T,N,M) {
    res := copy(a);
    inline neg(*res);
    return res;
}
operator - :: (a: DenseHeapMatrix($T)) -> DenseHeapMatrix(T) {
    res := copy(a); // ! alloc
    inline neg(*res);
    return res;
}

// ? by iterating through the raw data, there is no double-negation issue!
neg :: (a: *$A/MatrixType) {
    for raw(a) {
        set_raw(a, it_raw_index, -it);
    }
}


//              ###    ##     ## ########  ##    ##
//  ##   ##    ## ##    ##   ##  ##     ##  ##  ##
//   ## ##    ##   ##    ## ##   ##     ##   ####
// ######### ##     ##    ###    ########     ##
//   ## ##   #########   ## ##   ##           ##
//  ##   ##  ##     ##  ##   ##  ##           ##
//           ##     ## ##     ## ##           ##


// ? see `mul`
saxpy :: inline (c: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    inline mul(c,a,b);
}


// ########  ######## ########  ##     ## ##     ## ######## ########
// ##     ## ##       ##     ## ###   ### ##     ##    ##    ##
// ##     ## ##       ##     ## #### #### ##     ##    ##    ##
// ########  ######   ########  ## ### ## ##     ##    ##    ######
// ##        ##       ##   ##   ##     ## ##     ##    ##    ##
// ##        ##       ##    ##  ##     ## ##     ##    ##    ##
// ##        ######## ##     ## ##     ##  #######     ##    ########


// there are faster ways of doing this.
// check for compatible matrix types: e.g. this is not possible with triangular matrices
permute_rows :: (m: *$M/MatrixType, order: ..int) {
    inline permute_rows_default(m, ..order);
}

#scope_file
permute_rows_default :: (m: *$M/MatrixType, order: ..int) {
    #if CHECKS {
        // #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(order.count == rows(m));
    }

    for 0..rows(m)-1 {
        if order[it] <= it then continue;
        inline swap_rows(m, it, order[it]);
    }
}
#scope_export

permute_cols :: (m: *$M/MatrixType, order: ..int) {
    inline permute_cols_default(m, ..order);
}

#scope_file
permute_cols_default :: (m: *$M/MatrixType, order: ..int) {
    #if CHECKS {
        // #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(order.count == cols(m));
    }

    for 0..cols(m)-1 {
        if order[it] <= it then continue;
        inline swap_cols(m, it, order[it]);
    }
}
#scope_export




// there are faster ways of doing this.
// check for compatible matrix types: e.g. this is not possible with triangular matrices
permute_rows :: (m: *$M/MatrixType, order: $O/VectorType) {
    inline permute_rows_default(m, order);
}

#scope_file
permute_rows_default :: (m: *$M/MatrixType, order: $O/VectorType) {
    #if CHECKS {
        // #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(dim(order) == rows(m));
    }

    for 0..rows(m)-1 {
        if get(order, it) <= it then continue;
        inline swap_rows(m, it, get(order,it));
    }
}
#scope_export

permute_cols :: (m: *$M/MatrixType, order: $O/VectorType) {
    inline permute_cols_default(m, order);
}

#scope_file
permute_cols_default :: (m: *$M/MatrixType, order: $O/VectorType) {
    #if CHECKS {
        // #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(dim(order) == cols(m));
    }

    for 0..cols(m)-1 {
        if get(order,it) <= it then continue;
        swap_cols(m, it, get(order,it));
    }
}
#scope_export


//  ######  ##      ##    ###    ########
// ##    ## ##  ##  ##   ## ##   ##     ##
// ##       ##  ##  ##  ##   ##  ##     ##
//  ######  ##  ##  ## ##     ## ########
//       ## ##  ##  ## ######### ##
// ##    ## ##  ##  ## ##     ## ##
//  ######   ###  ###  ##     ## ##


swap_lin :: (m: *$M/MatrixType, i: int, j: int) {
    inline swap_lin_default(m,i,j);
}

#scope_file
swap_lin_default :: (m: *$M/MatrixType, i: int, j: int) {
    tmp := get_lin(m, i);
    set_lin(m, i, get_lin(m, j));
    set_lin(m, j, tmp);
}
#scope_export



swap :: (m: *$M/MatrixType, from_r: int, from_c: int, to_r: int, to_c: int) {
    inline swap_default(m, from_r, from_c, to_r, to_c);
}

#scope_file
swap_default :: (m: *$M/MatrixType, from_r: int, from_c: int, to_r: int, to_c: int) {
    tmp := get(m, from_r, from_c);
    set(m, from_r, from_c, get(m, to_r, to_c));
    set(m, to_r, to_c, tmp);
}
#scope_export



swap_cols :: (m: *$M/MatrixType, i: int, j: int) {
    inline swap_cols_default(m,i,j);
}

#scope_file
swap_cols_default :: (m: *$M/MatrixType, i: int, j: int) {
    #if CHECKS {
        // #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense)\n");
        assert(i >= 0 && j >= 0);
        assert(i < cols(m) && j < cols(m));
    }

    if i == j then return;

    for r: 0..rows(m)-1 {
        tmp := get(m, r, i);
        set(m, r, i, get(m, r, j));
        set(m, r, j, tmp);
    }
}
#scope_export

swap_rows :: (m: *$M/MatrixType, i: int, j: int) {
    inline swap_rows_default(m,i,j);
}

#scope_file
swap_rows_default :: (m: *$M/MatrixType, i: int, j: int) {
    #if CHECKS {
        // #run assert(M.flags == .DENSE || M.flags == .DENSE_HEAP, "not a compatible matrix type (not dense) : %\n", M.flags);
        assert(i >= 0 && j >= 0);
        assert(i < rows(m) && j < rows(m));
    }

    if i == j then return;

    for c: 0..cols(m)-1 {
        tmp := get(m, i, c);
        set(m, i, c, get(m, j, c));
        set(m, j, c, tmp);
    }
}
#scope_export

// ######## ########     ###    ##    ##  ######  ########   #######   ######  ########
//    ##    ##     ##   ## ##   ###   ## ##    ## ##     ## ##     ## ##    ## ##
//    ##    ##     ##  ##   ##  ####  ## ##       ##     ## ##     ## ##       ##
//    ##    ########  ##     ## ## ## ##  ######  ########  ##     ##  ######  ######
//    ##    ##   ##   ######### ##  ####       ## ##        ##     ##       ## ##
//    ##    ##    ##  ##     ## ##   ### ##    ## ##        ##     ## ##    ## ##
//    ##    ##     ## ##     ## ##    ##  ######  ##         #######   ######  ########


transpose :: (m: DenseMatrix($T, $N, $M)) -> DenseMatrix(T, N, M) {
    res := copy(m);
    inline transpose(*res);
    return res;
}
transpose :: (m: DenseHeapMatrix($T)) -> DenseHeapMatrix(T) {
    res := copy(m); // ! alloc
    inline transpose(*res);
    return res;
}


transpose :: (m: *$M/MatrixType) {
    #if M.flags == .DENSE {
        inline transpose_dense(m);
    } else if M.flags == MATRIX_FLAGS.DENSE & MATRIX_FLAGS.HEAP {
        inline transpose_dense_heap(m);
    }
}

#scope_file
// ? only square matrices here
// ? maybe add an additional runtime variable `transpose: bool` to DenseMatrix?
transpose_dense :: (m: *DenseMatrix($T, $N, N)) {
    for i: 0..N-1 {
        for j: i+1..N-1 {
            inline swap(m, i, j, j, i);
        }
    }
}

// ? https://www.wikiwand.com/en/In-place_matrix_transposition
transpose_dense_heap :: (m: *DenseHeapMatrix) {
    if rows(m) == cols(m) {
        for i: 0..rows(m)-1 {
            for j: i+1..cols(m)-1 {
                inline swap(m, i, j, j, i);
            }
        }
    } else {
        //     0 1 2 3
        //     4 5 6 7
        //     8 9 a b
        // ->  
        //     0 4 8     0 4 8 1
        //     1 5 9  ~  5 9 2 6
        //     2 6 a     a 3 7 b
        //     3 7 b 

        r := rows(m);
        c := cols(m);
        for m {
            pi := ifx it_index == r*c-1 then r*c-1 else (r*it_index) % (r*c-1);
            if pi <= it_index then continue; // ! this is my tweak. I hope it works?!
            inline swap_lin(m, pi, it_index);
        }
        m.rows = c;
        m.cols = r;
    }

}
#scope_export



//  ######   #######  ##    ##       ## ##     ##  ######      ###    ######## ########
// ##    ## ##     ## ###   ##       ## ##     ## ##    ##    ## ##      ##    ##
// ##       ##     ## ####  ##       ## ##     ## ##         ##   ##     ##    ##
// ##       ##     ## ## ## ##       ## ##     ## ##   #### ##     ##    ##    ######
// ##       ##     ## ##  #### ##    ## ##     ## ##    ##  #########    ##    ##
// ##    ## ##     ## ##   ### ##    ## ##     ## ##    ##  ##     ##    ##    ##
//  ######   #######  ##    ##  ######   #######   ######   ##     ##    ##    ########

conjugate :: (m: DenseMatrix($T, $N, $M)) -> DenseMatrix(T, N, M) {
    res := copy(m);
    inline conjugate(*res);
    return res;
}
conjugate :: (m: DenseHeapMatrix($T)) -> DenseHeapMatrix(T) {
    res := copy(m); // ! alloc
    inline conjugate(*res);
    return res;
}

conjugate :: (m: *$M/MatrixType) {
    #if #run !is_non_real_number(M.data_type) {
        return;
    } else {
        for raw(m) {
            set_raw(m, it_raw_index, conjugate(it));
        }
    }
}


// ########     ###     ######    ######   ######## ########
// ##     ##   ## ##   ##    ##  ##    ##  ##       ##     ##
// ##     ##  ##   ##  ##        ##        ##       ##     ##
// ##     ## ##     ## ##   #### ##   #### ######   ########
// ##     ## ######### ##    ##  ##    ##  ##       ##   ##
// ##     ## ##     ## ##    ##  ##    ##  ##       ##    ##
// ########  ##     ##  ######    ######   ######## ##     ##

dagger :: conjugate_transpose;
hermitian_conjugate :: conjugate_transpose;
hermitian_transpose :: conjugate_transpose;

conjugate_transpose :: (m: DenseMatrix($T, $N, $M)) -> DenseMatrix(T, N, M) {
    res := copy(m);
    inline conjugate_transpose(*res);
    return res;
}
conjugate_transpose :: (m: DenseHeapMatrix($T)) -> DenseHeapMatrix(T) {
    res := copy(m); // ! alloc
    inline conjugate_transpose(*res);
    return res;
}

conjugate_transpose :: (m: *$M/MatrixType) {
    inline transpose(m);
    inline conjugate(m);
}



// ######## ######## ##    ##  ######   #######  ########
//    ##    ##       ###   ## ##    ## ##     ## ##     ##
//    ##    ##       ####  ## ##       ##     ## ##     ##
//    ##    ######   ## ## ##  ######  ##     ## ########
//    ##    ##       ##  ####       ## ##     ## ##   ##
//    ##    ##       ##   ### ##    ## ##     ## ##    ##
//    ##    ######## ##    ##  ######   #######  ##     ##


tensor :: (a: DenseMatrix($A, $AN, $AM), b: DenseMatrix($B, $BN, $BM)) -> DenseMatrix(UpCast(A,B).T, AN*BN, AM*BM) {
    res : DenseMatrix(UpCast(A,B).T, AN*BN, AM*BM) = ---;
    inline tensor(a,b, *res);
    return res;
}
tensor :: (a: DenseHeapMatrix($A), b: DenseHeapMatrix($B)) -> DenseHeapMatrix(UpCast(A,B).T) {
    res := dhmat(UpCast(A,B).T, rows(a)*rows(b), cols(a)*cols(b), false);
    inline tensor(a,b, *res);
    return res;
}

// https://www.wikiwand.com/en/Tensor_product
tensor :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    inline tensor_default(res,a,b);
}

#scope_file
tensor_default :: (res: *$C/MatrixType, a: $A/MatrixType, b: $B/MatrixType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type) && is_minor_number_type(C.data_type, B.data_type));
        assert(rows(res) == rows(a)*rows(b) && cols(res) == cols(a)*cols(b), "dimensions mismatch: needed (%,%), but res(%,%)", rows(a)*rows(b), cols(a)*cols(b), rows(res), cols(res));
    }

    for r: 0..rows(res)-1 {
        ar := r/rows(b);
        br := r % rows(b);
        
        for c: 0..cols(res)-1 {
            ac := c/cols(b);
            bc := c % cols(b);

            set(res, r, c, get(a, ar, ac) * get(b, br, bc));
        }
    }
}
#scope_export



// ##    ##  #######  ########  ##     ##
// ###   ## ##     ## ##     ## ###   ###
// ####  ## ##     ## ##     ## #### ####
// ## ## ## ##     ## ########  ## ### ##
// ##  #### ##     ## ##   ##   ##     ##
// ##   ### ##     ## ##    ##  ##     ##
// ##    ##  #######  ##     ## ##     ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Lemma 3.5.4 Matrix Norms, page 176

// this needs the largest eigenvalue
// norm_2 :: inline (v: Vector) -> float64 {
//     return norm(v, 2.0);
// }


norm_1 :: inline (m: $M/MatrixType) -> float64 {
    return norm_1_default(m);
}

#scope_file
norm_1_default :: (m: $M/MatrixType) -> float64 {
    res : float64 = 0.0;

    for j: 0..cols(m)-1 {
        sum : float64 = 0.0;
        for i: 0..rows(m)-1 {
            sum += abs(get(m,i,j));
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}
#scope_export




norm_inf :: inline (m: $M/MatrixType) -> float64 {
    return norm_inf_default(m);
}

#scope_file
norm_inf_default :: (m: $M/MatrixType) -> float64 {
    res : float64 = 0.0;

    for i: 0..rows(m)-1 {
        sum : float64 = 0.0;
        for j: 0..cols(m)-1 {
            sum += abs(get(m,i,j));
        }
        res = ifx res < sum then sum else res;
    }

    return res;
}
#scope_export



norm_frobenius :: inline (m: $M/MatrixType) -> float64 {
    return norm_frobenius_default(m);
}

#scope_file
norm_frobenius_default :: (m: $M/MatrixType) -> float64 {
    sum : float64 = 0.0;

    for m {
        sum += abs_sq(it);
    }

    return sqrt(sum);
}
#scope_export


norm_max :: inline (m: $M/MatrixType) -> float64 {
    return norm_max_default(m);
}

#scope_file
norm_max_default :: (m: $M/MatrixType) -> float64 {
    res : float64 = 0.0;
    for m {
        a := abs(it);
        if a > res then res = a;
    }
    return res;
}
#scope_export





////////////////////////////////////////////////////////////////////////////////////////
/*
      ######## ########  ######  ########  ######
           ##    ##       ##    ##    ##    ##    ##
           ##    ##       ##          ##    ##
           ##    ######    ######     ##     ######
           ##    ##             ##    ##          ##
           ##    ##       ##    ##    ##    ##    ##
           ##    ########  ######     ##     ######
*/
////////////////////////////////////////////////////////////////////////////////////////

#scope_module

test_matrix :: () {

    println_push("Matrix", color = .FG_WHITE);

    {
        println_push("DenseMatrix", color = .FG_GREEN);
        {
            m : DenseMatrix(float64, 3, 3);
            println("%", pstr(m));
        }
        {
            a :: 3;
            m : DenseMatrix(float64, a, a);
            println("%", pstr(m));
        }



        {
            m : DenseMatrix(float64, 3, 3);
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }
        {
            a :: 3;
            m : DenseMatrix(float64, a, a);
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }


        {
            m : DenseMatrix(C64, 3, 3);
            mat_init_all(*m, cpx(42.0));
            set(*m,1,2,R);
            assert( get(m,1,2) == cpx(1.0, 0.0) );
            println("m = %", pstr(m));
        }

        {
            m : DenseMatrix(C64, 4, 4) = ---;
            mat_init_unit(*m);
            println("m = %", pstr(m));
            assert(get(m,0,1) == C_ZERO);
            assert(get(m,1,1) == R);
            assert(get(m,2,1) == C_ZERO);
            assert(get(m,3,3) == R);
        }

        {
            m : DenseMatrix(C64, 3, 3);
            set(*m,0,2,R);
            transpose(*m);
            assert(get(m,2,0) == R);
        }

        {
            a : DenseMatrix(float64, 2, 2);
            b : DenseMatrix(float64, 2, 2);
            mat_init_unit(*a);
            mat_init_unit(*b);
            println("a = %", pstr(a));
            println("b = %", pstr(b));
            add(*a, b);
            println("a = %", pstr(a));
            mul(*b, 2.0);
            println("2 b = %", pstr(b));
            assert(a == b);

            add(*a, b);
            mul(*b, 2);
            assert(a == b);
        }

        {
            m : DenseMatrix(C64, 2, 2);
            mat_init(*m, R, 2*R, 3*R, 4*R);
            println("m = %", pstr(m));
            res : DenseMatrix(C64, 2, 2);
            mul(*res, m, m);
            println("res = %", pstr(res));
        }
        {
            m : DenseMatrix(C64, 2, 2);
            mat_init(*m, R, 2.0*R, 3.0*R, 4.0*R);
            println("m = %", pstr(m));
            res : DenseMatrix(C64, 2, 2);
            mul(*res, m, m);
            println("res = %", pstr(res));
        }
        // {
        //     m := mat(1, 3, 1.0, 2.0, 3.0);
        //     print("m = %\n", m);
        //     mt := transpose(m);
        //     m = m * mt;
        //     print("m = %\n", m);
        // }
        {
            m : DenseMatrix(C64, 2, 2);
            mat_init(*m, 
                C_ZERO, -I, 
                I, C_ZERO
            );
            k : DenseMatrix(C64, 2, 2);
            copy(*k, m);
            dagger(*k);
            assert( k == m );
        }


        {
            m : DenseMatrix(C64, 2, 2);
            mat_init(*m, R, I, 2*R-I, 2*R + 3*I);
            println("m = %", pstr(m));
            t1 : DenseMatrix(C64, 4, 4);
            tensor(*t1, m, m);
            println("t1 = %", pstr(t1));
            t2 : DenseMatrix(C64, 4, 4);
            mat_init(*t2,
                R            ,             I,             I,   -R           ,
                2.0 -       I, 2.0 + 3.0 * I, 1.0 + 2.0 * I, -3.0 +  2.0 * I,
                2.0 -       I, 1.0 + 2.0 * I, 2.0 + 3.0 * I, -3.0 +  2.0 * I,
                3.0 - 4.0 * I, 7.0 + 4.0 * I, 7.0 + 4.0 * I, -5.0 + 12.0 * I
            );
            println("t2 = %", pstr(t2));
            assert(t1 == t2);
        }

        {
            M : DenseMatrix(C64, 3, 3);
            mat_init(*M,
                C_ZERO, R, 2*R,
                3*R, 4*R, 5*R,
                6*R, 7*R, 8*R
            );

            S := mat_view(*M, 2, 2, 1, 1);
            
            res : DenseMatrix(C64, 2, 2);
            mat_init(*res, 4*R, 5*R, 7*R, 8*R);
            
            assert(S == res);
        }

        // {
        //     print("%\n", MATRIX_FLAGS.UNIT_LOWER_TRIANGULAR);
        //     m := make(DenseMatrix(float64, 4, 4)) = ---;
        //     m.flags = .UNIT_LOWER_TRIANGULAR;
        //     assert( is(m, .TRAPEZOIDAL) );
        // }
    }



    {
        println_push("DenseHeapMatrix", color = .FG_GREEN);
        {
            println_push("test 1", color = .FG_YELLOW);
            m : DenseHeapMatrix(float64);
            defer free(m);
            println("%", pstr(m));
        }
        {
            println_push("test 2", color = .FG_YELLOW);
            a := 3;
            m := dhmat(float64, a, a);
            defer free(m);
            println("%", pstr(m));
        }



        {
            println_push("test 3", color = .FG_YELLOW);
            m := dhmat(float64, 3, 3);
            defer free(m);
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }
        {
            println_push("test 4", color = .FG_YELLOW);
            a := 3;
            m := dhmat(float64, a, a);
            defer free(m);
            mat_init_all(*m, 42.0);
            println("%", pstr(m));
        }


        {
            println_push("test 5", color = .FG_YELLOW);
            m := dhmat(C64, 3, 3, false);
            defer free(m);
            mat_init_all(*m, cpx(42.0));
            set(*m,1,2,R);
            assert( get(m,1,2) == cpx(1.0, 0.0) );
            println("m = %", pstr(m));
        }

        {
            println_push("test 6", color = .FG_YELLOW);
            m := dhmat(C64, 4, 4);
            defer free(m);
            mat_init_unit(*m);
            println("m = %", pstr(m));
            assert(get(m,0,1) == C_ZERO);
            assert(get(m,1,1) == R);
            assert(get(m,2,1) == C_ZERO);
            assert(get(m,3,3) == R);
        }

        {
            println_push("test 7", color = .FG_YELLOW);
            m := dhmat(C64, 3, 3);
            defer free(m);
            set(*m,0,2,R);
            println("transposing: ");
            transpose(*m);
            assert(get(m,2,0) == R);
        }

        {
            println_push("test 8", color = .FG_YELLOW);
            a := dhmat(float64, 2, 2);
            defer free(a);
            b := dhmat(float64, 2, 2);
            defer free(b);
            mat_init_unit(*a);
            mat_init_unit(*b);
            println("a = %", pstr(a));
            println("b = %", pstr(b));
            add(*a, b);
            println("a = %", pstr(a));
            mul(*b, 2.0);
            println("2 b = %", pstr(b));
            assert(a == b);

            add(*a, b);
            mul(*b, 2.0);
            assert(a == b);
        }

        {
            println_push("test 9", color = .FG_YELLOW);
            m := dhmat(C64, 2, 2);
            defer free(m);
            mat_init(*m, R, 2*R, 3*R, 4*R);
            println("m = %", pstr(m));
            res := dhmat(C64, 2, 2);
            defer free(res);
            mul(*res, m, m);
            println("res = %", pstr(res));
        }
        {
            println_push("test 10", color = .FG_YELLOW);
            m := dhmat(C64, 2, 2);
            defer free(m);
            mat_init(*m, R, 2.0*R, 3.0*R, 4.0*R);
            println("m = %", pstr(m));
            res := dhmat(C64, 2, 2);
            defer free(res);
            mul(*res, m, m);
            println("res = %", pstr(res));
        }
        {
            println_push("test 11", color = .FG_YELLOW);
            m := dhmat(C64, 2, 2);
            defer free(m);
            mat_init(*m, 
                C_ZERO, -I, 
                I, C_ZERO
            );
            b := dhmat(C64, 2, 2);
            defer free(b);
            copy(*b, m);

            dagger(*b);
            assert( b == m );
        }


        {
            println_push("test 12", color = .FG_YELLOW);
            m := dhmat(C64, 2, 2);
            defer free(m);
            mat_init(*m, R, I, 2*R-I, 2*R + 3*I);
            println("m = %", pstr(m));
            t1 := dhmat(C64, 4, 4);
            tensor(*t1, m, m);
            println("t1 = %", pstr(t1));
            t2 := dhmat(C64, 4, 4);
            mat_init(*t2,
                R            ,             I,             I,   -R           ,
                2.0 -       I, 2.0 + 3.0 * I, 1.0 + 2.0 * I, -3.0 +  2.0 * I,
                2.0 -       I, 1.0 + 2.0 * I, 2.0 + 3.0 * I, -3.0 +  2.0 * I,
                3.0 - 4.0 * I, 7.0 + 4.0 * I, 7.0 + 4.0 * I, -5.0 + 12.0 * I
            );
            println("t2 = %", pstr(t2));
            assert(t1 == t2);
        }

        {
            println_push("test 13", color = .FG_YELLOW);
            M := dhmat(C64, 3, 3);
            defer free(M);
            mat_init(*M,
                C_ZERO, R, 2*R,
                3*R, 4*R, 5*R,
                6*R, 7*R, 8*R
            );

            S := mat_heap_view(*M, 2, 2, 1, 1);
            println("S = %", pstr(S));
            
            res := dhmat(C64, 2, 2);
            defer free(res);
            mat_init(*res, 4*R, 5*R, 7*R, 8*R);
            
            assert(S == res);
        }
    }

    test_mat_mul_mat();
}