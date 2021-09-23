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
    return str(v, print_options.offset, print_options.shift + 4, "│");
}

str :: (v: $V/VectorType, offset: int = 0, indent: int = 4, indent_character: string = " ", indent_every: int = 4) -> string {
    builder: String_Builder;
    defer free_buffers(*builder);

    init_string_builder(*builder);
    print_to_builder(*builder, "\n");

    shift :: (additional_indent: int = 0) #expand {
        for 0..offset-1 {
            append(*builder, " ");
        }
        for 0..indent-1+additional_indent {
            if it % indent_every == 0 {
                append(*builder, indent_character);
            } else {
                append(*builder, " ");
            }
        }
    }

    shift();

    append(*builder, "(");
    for v {
        append(*builder, str(it));
        if it_index != v.N-1 {
            append(*builder, ", ");
        }
    }
    append(*builder, ")\n");

    shift(-4);
    append(*builder, "└╴");

    return builder_to_string(*builder);
}



/*
#### ##    ## #### ######## ####    ###    ##       #### ######## ########
 ##  ###   ##  ##     ##     ##    ## ##   ##        ##       ##  ##
 ##  ####  ##  ##     ##     ##   ##   ##  ##        ##      ##   ##
 ##  ## ## ##  ##     ##     ##  ##     ## ##        ##     ##    ######
 ##  ##  ####  ##     ##     ##  ######### ##        ##    ##     ##
 ##  ##   ###  ##     ##     ##  ##     ## ##        ##   ##      ##
#### ##    ## ####    ##    #### ##     ## ######## #### ######## ########
*/


vec_cast :: ($type: Type, v: $V/VectorType) -> VectorCastHelper(V, type).T {
    res := make(VectorCastHelper(V, type).T); // initialize all values to zero, important for sparse?

    for raw(*v) {
        set_raw(*res, it_raw_index, ncast(type, it));
    }

    return res;
}

copy :: (v: $V/VectorType) -> V {
    res := make(V, false);

    for v {
        set_lin(*res, it_lin_index, it);
    }

    return res;
}

copy :: (to: *$V/VectorType, from: V) {
    for raw(*from) {
        set_raw(to, it_raw_index, it);
    }
}


vec_init :: inline (v: *$V/VectorType, values: ..$T) {
    #if CHECKS {
        #assert(is_minor_number_type(V.data_type, T));
        assert(values.count <= V.N, "too many values (%) for vector of % entries\n", values.count, V.N);
    }
    for values {
        set(v, it_index, ncast(V.data_type, it));
    }
}

vec_init_all :: inline (v: *$V/VectorType, value: $T) {
    #if CHECKS {
        #assert(is_minor_number_type(V.data_type, T));
    }
    for raw(v) {
        set_raw(v, it_raw_index, ncast(V.data_type, value));
    }
}




// inclusive left, exclusive right bounds
vec_subvector :: (v: $V/VectorType, $from: int, $to: int) -> VectorSizeHelper(V, to-from).T {
    #if CHECKS {
        #run assert(to > from, "bounds not possible: to(%) <= from(%)\n", to, from);
        #run assert(to <= V.N, "bound not possible: to(%) >= N(%)\n", to, V.N);
        #run assert(from >= 0, "bound not possible: from(%) < 0\n", from);
    }
    res := make(VectorSizeHelper(V, to-from).T, false);

    for res {
        set(*res, it_index, get(v, it_index + from));
    }

    return res;
}


/*
 #######  ########  ######## ########     ###    ########  #######  ########   ######
##     ## ##     ## ##       ##     ##   ## ##      ##    ##     ## ##     ## ##    ##
##     ## ##     ## ##       ##     ##  ##   ##     ##    ##     ## ##     ## ##
##     ## ########  ######   ########  ##     ##    ##    ##     ## ########   ######
##     ## ##        ##       ##   ##   #########    ##    ##     ## ##   ##         ##
##     ## ##        ##       ##    ##  ##     ##    ##    ##     ## ##    ##  ##    ##
 #######  ##        ######## ##     ## ##     ##    ##     #######  ##     ##  ######
*/



operator == :: (a: $V/VectorType, b: $B/VectorType) -> bool {
    #if V.N != B.N {
        return false;
    } else {

        for a {
            if it != get(b, it_index) then return false;
        }

        return true;
    }
}




operator * :: (a: $A/VectorType, b: $B/VectorType) -> UpCast(A.data_type,B.data_type).T {
    #if CHECKS {
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
    }
    res := zero(UpCast(A.data_type,B.data_type).T);
    for a {
        res += it * get(b, it_index);
    }
    return res;
}

operator * :: (a: $A/VectorType, b: $B) -> VectorCastHelper(A,B).T #modify {
    if !is_number(B) {B = null; return;}    
} #symmetric {
    res := make(VectorCastHelper(A,B).T, false);
    for a {
        set(*res, it_index, it * b);
    }
    return res;
}



mul :: inline (a: *$A/VectorType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    #if CHECKS {
        #assert(is_minor_number_type(A.data_type,B));
    }
    for raw(a) {
        set_raw(a, it_raw_index, it * b);
    }
}










outer_product :: (a: $A/VectorType, b: $B/VectorType) -> VectorMatrixHelper(A,B).T {
    res := make(VectorMatrixHelper(A,B).T, false);

    for res {
        set(*res, it_row, it_column, get(a, it_row) * get(b, it_column));
    }

    return res;
}

operator / :: inline (a: $A/VectorType, b: $B/VectorType) -> VectorMatrixHelper(A,B).T {
    return outer_product(a, b);
}

operator / :: (a: $A/VectorType, b: $B) -> VectorCastHelper(A,B).T #modify {
    if !is_number(B) {B = null; return;}    
} {
    res := make(VectorCastHelper(A,B).T, false);

    for a {
        set(*res, it_index, it / b);
    }

    return res;
}

div :: inline (a: *$A/VectorType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    #if CHECKS {
        #assert(is_minor_number_type(A,B));
    }

    for raw(*a) {
        set_raw(*a, it_raw_index, it / b);
    }
}






operator + :: (a: $A/VectorType, b: $B/VectorType) -> VectorVectorHelper(A,B).T {
    #if CHECKS {
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
    }
    res := make(VectorVectorHelper(A,B).T, false);

    for a {
        set(*res, it_index, it + get(b, it_index));
    }

    return res;
}

add :: (a: *$A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #assert(is_minor_number_type(A.data_type, B.data_type));
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
    }

    for raw(*b) {
        set(a, it_index, get(<<a, it_index) + it);
    }
}




operator - :: (a: $A/VectorType, b: $B/VectorType) -> VectorVectorHelper(A,B).T {
    #if CHECKS {
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
    }
    res := make(VectorVectorHelper(A,B).T, false);

    for a {
        set(*res, it_index, it - get(b, it_index));
    }

    return res;
}

sub :: (a: *$A/VectorType, b: $B/VectorType) {
    #if CHECKS {
        #assert(is_minor_number_type(A.data_type, B.data_type));
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
    }

    for raw(*b) {
        set(a, it_index, get(<<a, it_index) - it);
    }
}




operator - :: (a: $A/VectorType) -> A {
    res := make(A, false);

    for raw(*a) {
        set_raw(*res, it_raw_index, -it);
    }

    return res;
}

neg :: (a: *$A/VectorType) {
    for raw(a) {
        set_raw(a, it_raw_index, -it);
    }
}









cross :: (a: $A/VectorType, b: $B/VectorType) -> VectorVectorHelper(A,B).T {
    #if CHECKS {
        #run assert(A.N == 3, "vector not of same dimension 3 (%)\n", A.N);
        #run assert(B.N == 3, "vector not of same dimension 3 (%)\n", B.N);
    }
    res := make(VectorVectorHelper(A,B).T, false);
    set(*res, 0, get(a,1)*get(b,2) - get(a,2)*get(b,1));
    set(*res, 1, get(a,2)*get(b,0) - get(a,0)*get(b,2));
    set(*res, 2, get(a,0)*get(b,1) - get(a,1)*get(b,0));
    return res;
}


reflect :: (vec: $A/VectorType, normal: $B/VectorType) -> VectorVectorHelper(A,B).T {
    #if CHECKS {
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
    }
    tau := 2.0 / (conjugate(normal) * normal);
    alpha := tau * conjugate(normal) * vec;
    return vec - normal * alpha;
}










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
saxpy :: (a: *$A/VectorType, b: $B/VectorType, factor: $C) #modify {
    if !is_number(C) {C = null; return;}    
} {
    #if CHECKS {
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
        #assert(is_minor_number_type(A.data_type, B) && is_minor_number_type(A.data_type, C));
    }

    for a {
        set(a, it_index, it + factor.value * get(b, it_index));
    }
}



// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.3 Row-Oriented Gaxpy, page 5
gaxpy :: (y: *$Y/VectorType, A: $M/MatrixType, x: $X/VectorType) {
    #if CHECKS {
        #run assert(Y.N == M.N, "dimensions do not match (%, %)\n", Y.N, M.N);
        #run assert(X.N == M.M, "dimensions do not match (%, %)\n", X.N, M.M);
        #assert(is_minor_number_type(Y.data_type, M.data_type) && is_minor_number_type(Y.data_type, X.data_type));
    }

    for A {
        set(y, it_row, get(<<y, it_row) + it * get(x, it_column));
    }
}











/*
######## ########     ###    ##    ##  ######  ########
   ##    ##     ##   ## ##   ###   ## ##    ## ##
   ##    ##     ##  ##   ##  ####  ## ##       ##
   ##    ########  ##     ## ## ## ##  ######  ######
   ##    ##   ##   ######### ##  ####       ## ##
   ##    ##    ##  ##     ## ##   ### ##    ## ##
   ##    ##     ## ##     ## ##    ##  ######  ##
*/

permute :: (v: $A/VectorType, order: ..int) -> A {
    #if CHECKS {
        assert(A.N == order.count);
    }
    res := make(A, false);

    for res {
        set(*res, it_index, get(v, order[it_index]));
    }

    return res;
}

// TODO: check, if this in-place permutation is correct!
permute :: (v: *$V/VectorType, order: ..int) {
    #if CHECKS {
        assert(V.N == order.count);
    }

    for v {
        if order[it_index] <= it_index then continue;

        tmp := it;
        set(v, it_index, get(v, order[it_index]));
        set(v, order[it_index], tmp);
    }

    // * original code
    // for 0..N-1 {
    //     if order[it] <= it then continue;

    //     v.data[it], v.data[order[it]] = swap(v.data[it], v.data[order[it]]);
    // }
}

swap :: (v: *$V/VectorType, i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && i < V.N);
        assert(j >= 0 && j < V.N);
    }
    tmp := get(v, i);
    set(v, i, get(v, j));
    set(v, j, tmp);
}



conjugate :: (v: $V/VectorType) -> V {
    res := make(V, false);
    for raw(*v) {
        set_raw(*res, it_raw_index, conjugate(it));
    }
    return res;
}
conjugate :: (v: *$V/VectorType) {
    // * no-op for real numbers
    // TODO might update later when including quaternions, octonions
    #if V.data_type == Complex64 {
        for raw(v) {
            set_raw(v, it_raw_index, conjugate(it));
        }
    }
}




/*
        ##        ######   ######     ###    ##          ###    ########
         ##      ##    ## ##    ##   ## ##   ##         ## ##   ##     ##
          ##     ##       ##        ##   ##  ##        ##   ##  ##     ##
#######    ##     ######  ##       ##     ## ##       ##     ## ########
          ##           ## ##       ######### ##       ######### ##   ##
         ##      ##    ## ##    ## ##     ## ##       ##     ## ##    ##
        ##        ######   ######  ##     ## ######## ##     ## ##     ##
*/

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Definition 3.5.2 Norms, page 171
norm :: (v: $V/VectorType, $$n: float64) -> float64 {
    #if is_constant(n) {
        #if n == Math.FLOAT64_INFINITY {
            m := abs(get(v,0));
            for 1..V.N-1 {
                tmp := abs(get(v,it));
                m = ifx m < tmp then tmp else m;
            }
            return m;
        } else #if n == 1.0 {
            res : float64 = 0.0;
            for 0..V.N-1 {
                res := abs(get(v,it));
            }
            return res;
        } else #if n == 2.0 {
            return sqrt( (conjugate(v) * v).real, float64 );
        } else {
            return _generic_norm(v, n);
        }
    } else {
        return _generic_norm(v, n);
    }
}

#scope_file
_generic_norm :: (v: $V/VectorType, n: float64) -> float64 {
    res : float64 = 0.0;

    for v {
        res += pow(abs(it), n, float64);
    }

    return pow(res, 1.0/n, float64);
}
#scope_export


norm_2 :: inline (v: $V/VectorType) -> float64 {
    return norm(v, 2.0);
}

norm_1 :: inline (v: $V/VectorType) -> float64 {
    return norm(v, 1.0);
}

norm_inf :: inline (v: $V/VectorType) -> float64 {
    return norm(v, Math.FLOAT64_INFINITY);
}


angle :: (a: $A/VectorType, b: $B/VectorType) -> float64 {
    #if CHECKS {
        #run assert(A.N == B.N, "vectors not of same dimensions (%, %)\n", A.N, B.N);
    }
    return acos( abs(conjugate(a) * b) / (norm_2(a) * norm_2(b)) );
}













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
test_vector :: () {

    println_push("Vector", color = .FG_WHITE);

    {
        println_push("DenseVector", color = .FG_GREEN);    
        {
            v := make(DenseVector(float64, 5));
            println("v = %", pstr(v));
            vec_init(*v, 1.0, 2.0, 3.0, 4.0, 5.0);
            println("v = %", pstr(v));
            // print("%, type_of(v) = %, Vector(T) = Vector(%)\n", v, T.vector_type, T.data_type);
            a := v + v;
            println("a = %", pstr(v));
            a = a + v;
            println("a = %", pstr(v));
        }
    }

    {
        println_push("DenseHeapVector", color = .FG_GREEN);    
        {
            v := make(DenseHeapVector(float64, 5));
            defer free(v);
            println("v = %", pstr(v));
            vec_init(*v, 1.0, 2.0, 3.0, 4.0, 5.0);
            println("v = %", pstr(v));
            // print("%, type_of(v) = %, Vector(T) = Vector(%)\n", v, T.vector_type, T.data_type);
            a := v + v;
            println("a = %", pstr(v));
            a = a + v;
            println("a = %", pstr(v));
        }
    }
}

