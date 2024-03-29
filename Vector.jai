






//TODO ########  #######  ########   #######
//TODO    ##    ##     ## ##     ## ##     ##
//TODO    ##    ##     ## ##     ## ##     ##
//TODO    ##    ##     ## ##     ## ##     ##
//TODO    ##    ##     ## ##     ## ##     ##
//TODO    ##    ##     ## ##     ## ##     ##
//TODO    ##     #######  ########   #######


// TODO loop unrolling for stack-allocated vectors of size <= 4









/*
 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##
*/


pstr :: inline (v: $V/VectorType) -> string {
    return str(v, print_options.offset, print_options.shift, "│", print_options.indent_character_every);
}

str :: (
    v: $V/VectorType, 
    offset: int = 0, 
    indent: int = 4, 
    indent_character: string = " ", 
    indent_character_every: int = 4
) -> string {
    builder: String_Builder;
    defer free_buffers(*builder);

    init_string_builder(*builder);
    append(*builder, "\n");

    shift :: (ind: int = 0) #expand {
        for 0..offset-1 {
            append(*builder, " ");
        }
        for 0..ind-1 {
            if it % indent_character_every == 0 {
                append(*builder, indent_character);
            } else {
                append(*builder, " ");
            }
        }
    }

    shift(indent+4);
    append(*builder, "(");
    for v {
        append(*builder, str(it));
        if it_index != dim(v)-1 {
            append(*builder, ", ");
        }
    }
    append(*builder, ")\n");

    shift(indent);
    append(*builder, "└╴");

    return builder_to_string(*builder);
}

//  ######     ###     ######  ########         ########  #######
// ##    ##   ## ##   ##    ##    ##               ##    ##     ##
// ##        ##   ##  ##          ##               ##    ##     ##
// ##       ##     ##  ######     ##               ##    ##     ##
// ##       #########       ##    ##               ##    ##     ##
// ##    ## ##     ## ##    ##    ##               ##    ##     ##
//  ######  ##     ##  ######     ##    #######    ##     #######

// * Assuming `to` is initialized to zeroes.
// * User has to know, what type the result has to be!
cast_to :: (to: *$T/VectorType, from: $F/VectorType) {
    for raw(*from) {
        set(to, it_index, ncast(T.data_type, it));
    }
}

//  ######   #######  ########  ##    ##
// ##    ## ##     ## ##     ##  ##  ##
// ##       ##     ## ##     ##   ####
// ##       ##     ## ########     ##
// ##       ##     ## ##           ##
// ##    ## ##     ## ##           ##
//  ######   #######  ##           ##

// * Assuming `to` is initialized to zeroes.
// * User has to know, what type the result has to be!
copy :: (to: *$V/VectorType, from: V) {
    for raw(*from) {
        set(to, it_index, it);
    }
}

copy :: (v: DenseVector) -> type_of(v) {
    res : type_of(v) = ---;
    inline copy(*res, v);
    return res;
}
copy :: (v: DenseHeapVector) -> type_of(v) {
    res := dhvec(v.data_type, dim(v), false); // ! alloc
    inline copy(*res, v);
    return v;
}

// ########  #######  ##     ##    ###    ##        ######
// ##       ##     ## ##     ##   ## ##   ##       ##    ##
// ##       ##     ## ##     ##  ##   ##  ##       ##
// ######   ##     ## ##     ## ##     ## ##        ######
// ##       ##  ## ## ##     ## ######### ##             ##
// ##       ##    ##  ##     ## ##     ## ##       ##    ##
// ########  ##### ##  #######  ##     ## ########  ######

operator == :: (a: $V/VectorType, b: $B/VectorType) -> bool {
    return inline equals(a,b);
}

equals :: inline (a: $A/VectorType, b: $B/VectorType) -> bool {
    #if A.flags == .DENSE && B.flags == .DENSE {
        return inline equals_dense_dense(a,b);
    } else {
        // * default implementation
        return inline equals_default(a,b);
    }
}

#scope_file
equals_default :: (a: $A/VectorType, b: $B/VectorType) -> bool {
    if dim(a) != dim(b) then return false;

    for a {
        if it != get(b, it_index) then return false;
    }

    return true;
}
equals_dense_dense :: (a: DenseVector, b: DenseVector) -> bool {
    #if a.N != b.N {
        return false;
    } else {
        for a {
            if it != get(b, it_index) then return false;
        }
        return true;
    }
}
#scope_export

// ##     ## ########  ######                      ##     ## ########  ######
// ##     ## ##       ##    ##         ##   ##     ##     ## ##       ##    ##
// ##     ## ##       ##                ## ##      ##     ## ##       ##
// ##     ## ######   ##              #########    ##     ## ######   ##
//  ##   ##  ##       ##                ## ##       ##   ##  ##       ##
//   ## ##   ##       ##    ##    ###  ##   ##       ## ##   ##       ##    ##
//    ###    ########  ######     ###                 ###    ########  ######


mul_el :: (a: DenseVector, b: DenseVector) -> DenseVector(UpCast(a.data_type,b.data_type).T, a.N) {
    #assert(a.N == b.N);
    res : DenseVector(UpCast(a.data_type,b.data_type).T, N);
    inline mul_el(*res,a,b);
    return res;
}
mul_el :: (a: DenseHeapVector($A), b: DenseHeapVector($B)) -> DenseHeapVector(UpCast(a.data_type,b.data_type).T) {
    res := dhvec(UpCast(a.data_type,b.data_type).T, dim(a)); // ! alloc
    inline mul_el(*res,a,b);
    return res;
}


mul_el :: (c: *$C/VectorType, a: $A/VectorType, b: $B/VectorType) {
    inline mul_el_default(c,a,b);
}

mul_el :: (a: *$A/VectorType, b: $B/VectorType) {
    inline mul_el_default(a,<<a, b);
}

#scope_file
mul_el_default :: (c: *$C/VectorType, a: $A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B.data_type));
        assert(dim(a) == dim(c));
        assert(dim(b) == dim(c));
    }
    for c {
        set(c, it_index, it + get(a, it_index) * get(b, it_index));
    }
}
#scope_export


// ##     ## ########  ######               ##    ##     ## ########  ######
// ##     ## ##       ##    ##             ##     ##     ## ##       ##    ##
// ##     ## ##       ##                  ##      ##     ## ##       ##
// ##     ## ######   ##                 ##       ##     ## ######   ##
//  ##   ##  ##       ##                ##         ##   ##  ##       ##
//   ## ##   ##       ##    ##    ###  ##           ## ##   ##       ##    ##
//    ###    ########  ######     ### ##             ###    ########  ######

div_el :: (a: DenseVector, b: DenseVector) -> DenseVector(UpCast(a.data_type,b.data_type).T, a.N) {
    #assert(a.N == b.N);
    res : DenseVector(UpCast(a.data_type,b.data_type).T, a.N);
    inline div_el(*res,a,b);
    return res;
}
div_el :: (a: DenseHeapVector, b: DenseHeapVector) -> DenseHeapVector(UpCast(a.data_type,b.data_type).T) {
    #if CHECKS {
        assert(dim(a) == dim(b));
    }
    res := dhvec(UpCast(a.data_type,b.data_type).T, dim(a)); // ! alloc
    inline div_el(*res,a,b);
    return res;
}


div_el :: (c: *$C/VectorType, a: $A/VectorType, b: $B/VectorType) {
    inline div_el_default(c,a,b);
}

div_el :: (a: *$A/VectorType, b: $B/VectorType) {
    inline div_el_default(a,<<a, b);
}

#scope_file
div_el_default :: (c: *$C/VectorType, a: $A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #assert(is_minor_number_type(C.data_type, A.data_type));
        #assert(is_minor_number_type(C.data_type, B.data_type));
        assert(dim(a) == dim(c));
        assert(dim(b) == dim(c));
    }
    for c {
        set(c, it_index, it + get(a, it_index) / get(b, it_index));
    }
}
#scope_export


// ########   #######  ########
// ##     ## ##     ##    ##
// ##     ## ##     ##    ##
// ##     ## ##     ##    ##
// ##     ## ##     ##    ##
// ##     ## ##     ##    ##
// ########   #######     ##

// ? by default, the left vector will be conjugated, 
// ? like `⟨ψ|φ⟩ = ∑ₙ conjugate(ψₙ) φₙ` in quantum mechanics
operator * :: (a: $A/VectorType, b: $B/VectorType) -> UpCast(A.data_type,B.data_type).T {
    return inline dot(a,b);
}

dot :: inline (a: $A/VectorType, b: $B/VectorType, $conjugate: bool = true) -> UpCast(A.data_type,B.data_type).T {
    #if A.flags == .DENSE && B.flags == .DENSE {
        return inline dot_dense_dense(a,b,conjugate);
    } else {
        // * default implementation
        return inline dot_default(a,b,conjugate);
    }
}
#scope_file
dot_default :: (a: $A/VectorType, b: $B/VectorType, $conjugate: bool = true) -> UpCast(A.data_type,b.data_type).T {
    #if CHECKS {
        assert(dim(a) == dim(b), "vectors not of same dimensions (%, %)\n", dim(a), dim(b));
    }
    res := zero(UpCast(A.data_type,b.data_type).T);
    for a {
        #if conjugate {
            res += conjugate(it) * get(b, it_index);
        } else {
            res += it * get(b, it_index);
        }
    }
    return res;
}
dot_dense_dense :: (a: DenseVector, b: DenseVector, $conjugate: bool = true) -> UpCast(a.data_type,b.data_type).T {
    #assert(a.N == b.N);
    res := zero(UpCast(a.data_type,b.data_type).T);
    for a {
        #if conjugate {
            res += conjugate(it) * get(b, it_index);
        } else {
            res += it * get(b, it_index);
        }
    }
    return res;
}
#scope_export



//  ######   ######     ###    ##          ##     ## ##     ## ##
// ##    ## ##    ##   ## ##   ##          ###   ### ##     ## ##
// ##       ##        ##   ##  ##          #### #### ##     ## ##
//  ######  ##       ##     ## ##          ## ### ## ##     ## ##
//       ## ##       ######### ##          ##     ## ##     ## ##
// ##    ## ##    ## ##     ## ##          ##     ## ##     ## ##
//  ######   ######  ##     ## ########    ##     ##  #######  ########


operator * :: (a: DenseVector, b: $B) -> DenseVector(UpCast(a.data_type,B).T, a.N) #modify {
    if !is_number(B) {B = null; return;}    
} #symmetric {
    res : DenseVector(UpCast(a.data_type,B).T,a.N) = ---;
    inline copy(*res, a);
    inline mul(*res, b);
    return res;
}
operator * :: (a: DenseHeapVector, b: $B) -> DenseHeapVector(UpCast(a.data_type,B).T) #modify {
    if !is_number(B) {B = null; return;}    
} #symmetric {
    res := dhvec(UpCast(a.data_type,B).T, dim(a), false); // ! alloc
    inline copy(*res, a);
    inline mul(*res, b);
    return res;
}


mul :: (a: *$A/VectorType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    inline mul_default(a,b);
}

#scope_file
mul_default :: (a: *$A/VectorType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    #if CHECKS {
        #assert(is_minor_number_type(A.data_type,B));
    }
    for raw(a) {
        set_raw(a, it_raw_index, it * b);
    }
}
#scope_export






operator / :: inline (a: DenseVector, b: $B) -> DenseVector(UpCast(a.data_type,B).T, a.N) #modify {
    if !is_number(B) {B = null; return;}    
} {
    res : DenseVector(UpCast(a.data_type,B).T,a.N) = ---;
    inline copy(*res, a);
    inline div(*res, b);
    return res;
}
operator / :: inline (a: DenseHeapVector, b: $B) -> DenseHeapVector(UpCast(a.data_type,B).T) #modify {
    if !is_number(B) {B = null; return;}    
} {
    res := dhvec(UpCast(a.data_type,B).T, dim(a), false); // ! alloc
    inline copy(*res, a);
    inline div(*res, b);
    return res;
}


div :: (a: *$A/VectorType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    inline div_default(a,b);
}

#scope_file
div_default :: (a: *$A/VectorType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    ib := inv(b);
    inline mul(a, ib);
}
#scope_export

//  #######  ##     ## ######## ######## ########
// ##     ## ##     ##    ##    ##       ##     ##
// ##     ## ##     ##    ##    ##       ##     ##
// ##     ## ##     ##    ##    ######   ########
// ##     ## ##     ##    ##    ##       ##   ##
// ##     ## ##     ##    ##    ##       ##    ##
//  #######   #######     ##    ######## ##     ##


operator / :: (a: DenseVector, b: DenseVector) -> DenseMatrix(UpCast(a.data_type,b.data_type).T, a.N, b.N) {
    res : DenseMatrix(UpCast(a.data_type,b.data_type).T, a.N, b.N);
    inline outer_product(*res,a,b);
    return res;
}
operator / :: (a: DenseHeapVector, b: DenseHeapVector) -> DenseHeapMatrix(UpCast(a.data_type,b.data_type).T) {
    res := dhmat(UpCast(a.data_type,b.data_type).T, dim(a), dim(b));
    inline outer_product(*res,a,b);
    return res;
}

// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? 1.1.9 The Outer Product Update, p.7
// ? calculates M = M + a*bᵀ 
outer_product :: (M: *$C/MatrixType, a: $A/VectorType, b: $B/VectorType) {
    inline outer_product_default(a,b,M);
}

#scope_file
outer_product_default :: inline (M: *$C/MatrixType, a: $A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type));
        #run assert(is_minor_number_type(C.data_type, B.data_type));
        assert(rows(M) == dim(a));
        assert(cols(M) == dim(b));
    }
    for M {
        set(*M, it_row, it_column, it + get(a, it_row) * conjugate(get(b, it_column)));
    }
}
#scope_export

//    ###    ########  ########      ######  ##     ## ########
//   ## ##   ##     ## ##     ##    ##    ## ##     ## ##     ##
//  ##   ##  ##     ## ##     ##    ##       ##     ## ##     ##
// ##     ## ##     ## ##     ##     ######  ##     ## ########
// ######### ##     ## ##     ##          ## ##     ## ##     ##
// ##     ## ##     ## ##     ##    ##    ## ##     ## ##     ##
// ##     ## ########  ########      ######   #######  ########


operator + :: (a: DenseVector, b: DenseVector) -> DenseVector(UpCast(a.data_type,b.data_type).T, a.N) {
    #assert(a.N == b.N);
    res : DenseVector(UpCast(a.data_type,b.data_type).T, a.N) = ---;
    inline cast_to(*res, a);
    inline add(*res, b);
    return res;
}
operator + :: (a: DenseHeapVector, b: DenseHeapVector) -> DenseHeapVector(UpCast(a.data_type,b.data_type).T) {
    res := dhvec(UpCast(a.data_type,b.data_type).T, a.N, false);
    inline cast_to(*res, a);
    inline add(*res, b);
    return res;
}


add :: (a: *$A/VectorType, b: $B/VectorType) {
    inline add_default(a,b);
}

#scope_file
add_default :: (a: *$A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #assert(is_minor_number_type(A.data_type, B.data_type));
        assert(dim(a) == dim(b), "vectors not of same dimensions (%, %)\n", dim(a), dim(b));
    }

    for raw(*b) {
        set(a, it_index, get(a, it_index) + it);
    }
}
#scope_export




operator - :: (a: DenseVector, b: DenseVector) -> DenseVector(UpCast(a.data_type,b.data_type).T, a.N) {
    #assert(a.N == b.N);
    res : DenseVector(UpCast(a.data_type,b.data_type).T, a.N) = ---;
    inline cast_to(*res, a);
    inline sub(*res, b);
    return res;
}
operator - :: (a: DenseHeapVector, b: DenseHeapVector) -> DenseHeapVector(UpCast(a.data_type,b.data_type).T) {
    res := dhvec(UpCast(a.data_type,b.data_type).T, a.N, false);
    inline cast_to(*res, a);
    inline sub(*res, b);
    return res;
}


sub :: (a: *$A/VectorType, b: $B/VectorType) {
    inline sub_default(a,b);
}

#scope_file
sub_default :: (a: *$A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #assert(is_minor_number_type(A.data_type, B.data_type));
        assert(dim(a) == dim(b), "vectors not of same dimensions (%, %)\n", dim(a), dim(b));
    }

    for raw(*b) {
        set(a, it_index, get(a, it_index) - it);
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


operator - :: (a: DenseVector) -> type_of(a) {
    res := copy(a);
    inline neg(*res);
    return res;
}
operator - :: (a: DenseHeapVector) -> type_of(a) {
    res := copy(a);
    inline neg(*res);
    return res;
}

neg :: (a: *$A/VectorType) {
    inline neg_default(a);
}

#scope_file
neg_default :: (a: *$A/VectorType) {
    for raw(a) {
        set_raw(a, it_raw_index, -it);
    }
}
#scope_export


//  ######  ########   #######   ######   ######
// ##    ## ##     ## ##     ## ##    ## ##    ##
// ##       ##     ## ##     ## ##       ##
// ##       ########  ##     ##  ######   ######
// ##       ##   ##   ##     ##       ##       ##
// ##    ## ##    ##  ##     ## ##    ## ##    ##
//  ######  ##     ##  #######   ######   ######


cross :: (a: DenseVector, b: DenseVector) -> DenseVector(UpCast(a.data_type,b.data_type).T, 3) {
    #assert(a.N == 3 && b.N == 3);
    res : DenseVector(UpCast(a.data_type,b.data_type).T, 3);
    inline cross(*res,a,b);
    return res;
}
cross :: (a: DenseHeapVector, b: DenseHeapVector) -> DenseHeapVector(UpCast(a.data_type,b.data_type).T) {
    res := dhvec(UpCast(a.data_type,b.data_type).T, 3, false);
    inline cross(*res,a,b);
    return res;
}

// ! is this correct for complex vectors, or am I missing a conjugate somewhere?
cross :: (res: *$C/VectorType, a: $A/VectorType, b: $B/VectorType) {
    inline cross_default(res,a,b);
}

#scope_file
cross_default :: (res: *$C/VectorType, a: $A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type));
        #run assert(is_minor_number_type(C.data_type, B.data_type));
        assert(dim(a) == 3, "vector a not of same dimension 3 (%)\n", dim(a));
        assert(dim(b) == 3, "vector b not of same dimension 3 (%)\n", dim(b));
        assert(dim(res) == 3, "vector res not of same dimension 3 (%)\n", dim(res));
    }
    set(*res, 0, get(a,1)*get(b,2) - get(a,2)*get(b,1));
    set(*res, 1, get(a,2)*get(b,0) - get(a,0)*get(b,2));
    set(*res, 2, get(a,0)*get(b,1) - get(a,1)*get(b,0));
    return res;
}
#scope_export

// ########  ######## ######## ##       ########  ######  ########
// ##     ## ##       ##       ##       ##       ##    ##    ##
// ##     ## ##       ##       ##       ##       ##          ##
// ########  ######   ######   ##       ######   ##          ##
// ##   ##   ##       ##       ##       ##       ##          ##
// ##    ##  ##       ##       ##       ##       ##    ##    ##
// ##     ## ######## ##       ######## ########  ######     ##

// * default implementation
reflect :: (res: *$C/VectorType, vec: $A/VectorType, normal: $B/VectorType) {
    inline reflect_default(res, vec, normal);
}

#scope_file
reflect_default :: (res: *$C/VectorType, vec: $A/VectorType, normal: $B/VectorType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type));
        #run assert(is_minor_number_type(C.data_type, B.data_type));
        assert(dim(a) == dim(b) && dim(a) == dim(res), "vectors not of same dimensions (%, %, %)\n", dim(a), dim(b), dim(res));
    }

    tau := 2 / norm_2(normal, true); // ? norm_2(normal, squared=true) = ⟨normal|normal⟩
    alpha := tau * (normal * vec);

    add(res, vec);
    saxpy(res, normal, -alpha);
}
#scope_export



/*
             ###    ##     ## ########  ##    ##
 ##   ##    ## ##    ##   ##  ##     ##  ##  ##
  ## ##    ##   ##    ## ##   ##     ##   ####
######### ##     ##    ###    ########     ##
  ## ##   #########   ## ##   ##           ##
 ##   ##  ##     ##  ##   ##  ##           ##
          ##     ## ##     ## ##           ##
*/

// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.2 Saxpy, page 4
// ? y = a * x + y
saxpy :: (y: *$A/VectorType, x: $B/VectorType, a: $C) #modify {
    if !is_number(C) {C = null; return;}    
} {
    inline saxpy_default(y, x, a);
}

#scope_file
saxpy_default :: (y: *$A/VectorType, x: $B/VectorType, a: $C) #modify {
    if !is_number(C) {C = null; return;}    
} {
    #if CHECKS {
        assert(dim(y) == dim(x), "vectors not of same dimensions (%, %)\n", dim(y), dim(x));
        #assert(is_minor_number_type(A.data_type, B) && is_minor_number_type(A.data_type, C));
    }

    for y {
        set(y, it_index, it + a * get(x, it_index));
    }
}
#scope_export



// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.3 Row-Oriented Gaxpy, page 5
// ? y = A * x + y
gaxpy :: (y: *$Y/VectorType, A: $M/MatrixType, x: $X/VectorType) {
    inline gaxpy_default(y, A, x);
}

#scope_file
gaxpy_default :: (y: *$Y/VectorType, A: $M/MatrixType, x: $X/VectorType) {
    #if CHECKS {
        assert(dim(y) == rows(A), "dimensions do not match (%, %)\n", dim(y), rows(A));
        assert(dim(x) == cols(A), "dimensions do not match (%, %)\n", dim(x), cols(A));
        #assert(is_minor_number_type(Y.data_type, M.data_type) && is_minor_number_type(Y.data_type, X.data_type));
    }

    for A {
        set(y, it_row, get(y, it_row) + it * get(x, it_column));
    }
}
#scope_export


// ########  ######## ########  ##     ## ##     ## ######## ########
// ##     ## ##       ##     ## ###   ### ##     ##    ##    ##
// ##     ## ##       ##     ## #### #### ##     ##    ##    ##
// ########  ######   ########  ## ### ## ##     ##    ##    ######
// ##        ##       ##   ##   ##     ## ##     ##    ##    ##
// ##        ##       ##    ##  ##     ## ##     ##    ##    ##
// ##        ######## ##     ## ##     ##  #######     ##    ########

// permute :: (v: $A/VectorType, order: ..int) -> A {
//     #if CHECKS {
//         assert(A.N == order.count);
//     }
//     res := make(A, false);
//     for res {
//         set(*res, it_index, get(v, order[it_index]));
//     }
//     return res;
// }

// TODO: check, if this in-place permutation is correct!
permute :: (v: *$V/VectorType, order: ..int) {
    inline permute_default(v, ..order);
}


#scope_file
permute_default :: (v: *$V/VectorType, order: ..int) {
    #if CHECKS {
        assert(dim(v) == order.count);
    }

    for v {
        if order[it_index] <= it_index then continue;
        swap(v, it_index, order[it_index]);
    }

    // example:
    // order 1 5 4 6 2 3
    //       1 2 3 4 5 6
    //    -> 1 5 3 4 2 6
    //    -> 1 5 4 3 2 6
    //    -> 1 5 4 6 2 3 done

    // * original code
    // for 0..N-1 {
    //     if order[it] <= it then continue;

    //     v.data[it], v.data[order[it]] = swap(v.data[it], v.data[order[it]]);
    // }
}
#scope_export


permute :: (v: *$V/VectorType, order: $O/VectorType) {
    inline permute_default(v,order);
}

#scope_file
permute_default :: (v: *$V/VectorType, order: $O/VectorType) {
    #if CHECKS {
        assert(dim(v) == dim(order));
    }

    for v {

        // if get(order, it_index) > it_index {
        //     swap(v, it_index, get(order, it_index));
        // }

        if get(order,it_index) <= it_index then continue;
        swap(v, it_index, get(order,it_index));
    }

    // example:
    // order 1 5 4 6 2 3
    //       1 2 3 4 5 6
    //    -> 1 5 3 4 2 6
    //    -> 1 5 4 3 2 6
    //    -> 1 5 4 6 2 3 done

    // * original code
    // for 0..N-1 {
    //     if order[it] <= it then continue;

    //     v.data[it], v.data[order[it]] = swap(v.data[it], v.data[order[it]]);
    // }
}
#scope_export


//  ######  ##      ##    ###    ########
// ##    ## ##  ##  ##   ## ##   ##     ##
// ##       ##  ##  ##  ##   ##  ##     ##
//  ######  ##  ##  ## ##     ## ########
//       ## ##  ##  ## ######### ##
// ##    ## ##  ##  ## ##     ## ##
//  ######   ###  ###  ##     ## ##

swap :: inline (v: *$V/VectorType, i: int, j: int) {
    inline swap_default(v,i,j);
}

#scope_file
swap_default :: (v: *$V/VectorType, i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && i < dim(v));
        assert(j >= 0 && j < dim(v));
    }
    tmp := get(v, i);
    set(v, i, get(v, j));
    set(v, j, tmp);
}
#scope_export




//  ######   #######  ##    ##       ## ##     ##  ######      ###    ######## ########
// ##    ## ##     ## ###   ##       ## ##     ## ##    ##    ## ##      ##    ##
// ##       ##     ## ####  ##       ## ##     ## ##         ##   ##     ##    ##
// ##       ##     ## ## ## ##       ## ##     ## ##   #### ##     ##    ##    ######
// ##       ##     ## ##  #### ##    ## ##     ## ##    ##  #########    ##    ##
// ##    ## ##     ## ##   ### ##    ## ##     ## ##    ##  ##     ##    ##    ##
//  ######   #######  ##    ##  ######   #######   ######   ##     ##    ##    ########

conjugate :: (v: DenseVector) -> type_of(v) {
    res := copy(v);
    inline conjugate(*res);
    return res;
}
conjugate :: (v: DenseHeapVector) -> type_of(v) {
    res := copy(v);
    inline conjugate(*res);
    return res;
}

conjugate :: (v: *$V/VectorType) {
    inline conjugate_default(v);
}

#scope_file
conjugate_default :: (v: *$V/VectorType) {
    // * no-op for real numbers
    #if is_non_real_number(V.data_type) {
        for raw(v) {
            set_raw(v, it_raw_index, conjugate(it));
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
// ? Definition 3.5.2 Norms, page 171
// TODO return other type than float64?
norm :: (v: $V/VectorType, $$l: float64 = 2.0) -> float64 {
    #if is_constant(l) {
        #if l == Math.FLOAT64_INFINITY {
            return inline norm_inf(v);
        } else #if l == 1.0 {
            return inline norm_1(v);
        } else #if l == 2.0 {
            return inline norm_2(v);
        } else {
            return inline norm_l(v, l);
        }
    } else {
        if l == Math.FLOAT64_INFINITY {
            return inline norm_inf(v);
        } else if l == 1.0 {
            return inline norm_1(v);
        } else if l == 2.0 {
            return inline norm_2(v);
        }
        return inline norm_l(v, l);
    }
}

norm_l :: (v: $V/VectorType, l: float64) -> float64 {
    return inline norm_l_default(v, l);
}
#scope_file
norm_l_default :: (v: $V/VectorType, l: float64) -> float64 {
    res : float64 = 0.0;
    for v {
        res += pow(abs(it), l, float64);
    }
    return pow(res, 1.0/l, float64);
}
#scope_export


norm_2 :: (v: $V/VectorType, $squared: bool = false) -> float64 {
    return inline norm_2_default(v, squared);
}
#scope_file
norm_2_default :: (v: $V/VectorType, $squared: bool = false) -> float64 {
    res : float64 = 0.0;
    for v {
        res += abs_sq(it);
    }
    #if squared {
        return res;
    } else {
        return sqrt(res, float64);
    }
}
#scope_export

norm_1 :: (v: $V/VectorType) -> float64 {
    return inline norm_1_default(v);
}
#scope_file
norm_1_default :: (v: $V/VectorType) -> float64 {
    res : float64 = 0.0;
    for v {
        res += abs(it);
    }
    return res;
}
#scope_export

norm_inf :: (v: $V/VectorType) -> float64 {
    return inline norm_inf_default(v);
}
#scope_file
norm_inf_default :: (v: $V/VectorType) -> float64 {
    m := abs(get(v,0));
    for 1..dim(v)-1 {
        tmp := abs(get(v,it));
        m = ifx m < tmp then tmp else m;
    }
    return m;
}
#scope_export



//    ###    ##    ##  ######   ##       ########
//   ## ##   ###   ## ##    ##  ##       ##
//  ##   ##  ####  ## ##        ##       ##
// ##     ## ## ## ## ##   #### ##       ######
// ######### ##  #### ##    ##  ##       ##
// ##     ## ##   ### ##    ##  ##       ##
// ##     ## ##    ##  ######   ######## ########


angle :: (a: $A/VectorType, b: $B/VectorType) -> float64 {
    return inline angle_default(a,b);
}
#scope_file
angle_default :: (a: $A/VectorType, b: $B/VectorType) -> float64 {
    #if CHECKS {
        #run assert(dim(a) == dim(b), "vectors not of same dimensions (%, %)\n", dim(a), dim(b));
    }
    return acos( (a*b) / (norm_2(a) * norm_2(b)) );
}
#scope_export










// ######## ########  ######  ########  ######
//    ##    ##       ##    ##    ##    ##    ##
//    ##    ##       ##          ##    ##
//    ##    ######    ######     ##     ######
//    ##    ##             ##    ##          ##
//    ##    ##       ##    ##    ##    ##    ##
//    ##    ########  ######     ##     ######

#scope_module
test_vector :: () {

    println_push("Vector", color = .FG_WHITE);

    {
        println_push("DenseVector", color = .FG_GREEN);    
        {
            v := dvec(float64, 5, 1, 2, 3, 4, 5);
            println("v = %", pstr(v));
            print("%, type_of(v) = %\n", v, type_of(v));

            add(*v, v);
            println("v = %", pstr(v));
            add(*v, v);
            println("v = %", pstr(v));
        }
    }

    {
        println_push("DenseHeapVector", color = .FG_GREEN);    
        {
            v := dhvec(float64, 5, 1,2,3,4,5);
            defer free(v);
            println("v = %", pstr(v));
            add(*v,v);
            println("v = %", pstr(v));
            add(*v,v);
            println("v = %", pstr(v));
        }
    }
}

