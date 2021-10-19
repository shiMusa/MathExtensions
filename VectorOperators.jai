

// ########  #######  ##     ##    ###    ##        ######
// ##       ##     ## ##     ##   ## ##   ##       ##    ##
// ##       ##     ## ##     ##  ##   ##  ##       ##
// ######   ##     ## ##     ## ##     ## ##        ######
// ##       ##  ## ## ##     ## ######### ##             ##
// ##       ##    ##  ##     ## ##     ## ##       ##    ##
// ########  ##### ##  #######  ##     ## ########  ######

operator == :: inline (a: $V/VectorType, b: $B/VectorType) -> bool {
    return equals(a,b);
}

equals :: inline (a: $A/VectorType, b: $B/VectorType) -> bool {
    #if A.flags == .DENSE && B.flags == .DENSE {
        return equals_dense_dense(a,b);
    } else {
        // * default implementation
        return equals_default(a,b);
    }
}

#scope_file
equals_default :: (a: $A/VectorType, b: $B/VectorType) -> bool {
    #if CHECKS {
        assert(dim(a) == dim(b), "vectors not of same dimension (%,%)\n", dim(a), dim(b));
    }

    for a {
        if it != get(b, it_index) then return false;
    }

    return true;
}
equals_dense_dense :: (a: DenseVector($T, $N), b: DenseVector($U, $M)) -> bool {
    #if N != M {
        return false;
    } else {
        for a {
            if it != get(b, it_index) then return false;
        }
        return true;
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

operator * :: inline (a: $A/VectorType, b: $B/VectorType) -> UpCast(A.data_type,B.data_type).T {
    return dot(a,b);
}
dot :: inline (a: $A/VectorType, b: $B/VectorType) -> UpCast(A.data_type,B.data_type).T {
    #if A.flags == .DENSE && B.flags == .DENSE {
        return dot_dense_dense(a,b);
    } else {
        // * default implementation
        return dot_default(a,b);
    }
}
#scope_file
dot_default :: (a: $A/VectorType, b: $B/VectorType) -> UpCast(A.data_type,b.data_type).T {
    #if CHECKS {
        assert(dim(a) == dim(b), "vectors not of same dimensions (%, %)\n", dim(a), dim(b));
    }
    res := zero(UpCast(A.data_type,b.data_type).T);
    for a {
        res += it * get(b, it_index);
    }
    return res;
}
dot_dense_dense :: (a: DenseVector($T, $N), b: DenseVector($U, $M)) -> UpCast(T,U).T {
    #if CHECKS {
        #run assert(N == M, "vectors not of same dimensions (%, %)\n", N, M);
    }
    res := zero(UpCast(T,U).T);
    for a {
        res += it * get(b, it_index);
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


div :: (a: *$A/VectorType, b: $B) #modify {
    if !is_number(B) {B = null; return;}    
} {
    inline div(a,b);
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
outer_product :: (a: $A/VectorType, b: $B/VectorType, res: *$C/MatrixType) {
    inline outer_product_default(a,b,res);
}

#scope_file
outer_product_default :: (a: $A/VectorType, b: $B/VectorType, res: *$C/MatrixType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type));
        #run assert(is_minor_number_type(C.data_type, B.data_type));
    }
    for res {
        set(*res, it_row, it_column, get(a, it_row) * get(b, it_column));
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
cross :: (a: $A/VectorType, b: $B/VectorType, res: *$C/VectorType) {
    inline cross_default(a,b,res);
}

#scope_file
cross_default :: (a: $A/VectorType, b: $B/VectorType, res: *$C/VectorType) {
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
reflect :: (vec: $A/VectorType, normal: $B/VectorType, res: *$C/VectorType) {
    inline reflect_default(vec, normal, res);
}

#scope_file
reflect_default :: (vec: $A/VectorType, normal: $B/VectorType, res: *$C/VectorType) {
    #if CHECKS {
        #run assert(is_minor_number_type(C.data_type, A.data_type));
        #run assert(is_minor_number_type(C.data_type, B.data_type));
        #run assert(dim(a) == dim(b) && dim(a) == dim(res), "vectors not of same dimensions (%, %, %)\n", dim(a), dim(b), dim(res));
    }
    cn = conjugate(normal);
    #if B.flags & .HEAP {
        defer free(cn);
    }

    tau := 2 / (cn * normal);
    alpha := tau * (cn * vec);

    copy(res, vec);
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
saxpy :: (a: *$A/VectorType, b: $B/VectorType, factor: $C) #modify {
    if !is_number(C) {C = null; return;}    
} {
    inline saxpy_default(a, b, factor);
}

#scope_file
saxpy_default :: (a: *$A/VectorType, b: $B/VectorType, factor: $C) #modify {
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
#scope_export



// ? Matrix Computations, 4th edition
// ? The Johns Hopkins University Press
// ? Algorithm 1.1.3 Row-Oriented Gaxpy, page 5
gaxpy :: (y: *$Y/VectorType, A: $M/MatrixType, x: $X/VectorType) {
    inline gaxpy_default(y, A, x);
}

#scope_file
gaxpy_default :: (y: *$Y/VectorType, A: $M/MatrixType, x: $X/VectorType) {
    #if CHECKS {
        #run assert(Y.N == M.N, "dimensions do not match (%, %)\n", Y.N, M.N);
        #run assert(X.N == M.M, "dimensions do not match (%, %)\n", X.N, M.M);
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
        assert(V.N == order.count);
    }

    for v {
        if order[it_index] <= it_index then continue;
        inline swap(v, it_index, order[it_index]);
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

swap :: (v: *$V/VectorType, i: int, j: int) {
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
            return norm_inf(v);
        } else #if l == 1.0 {
            return norm_1(v);
        } else #if l == 2.0 {
            return norm_2(v);
        } else {
            return norm_l(v, l);
        }
    } else {
        if l == Math.FLOAT64_INFINITY {
            return inline norm_inf(v);
        } else if l == 1.0 {
            return inline norm_1(v);
        } else if l == 2.0 {
            return inline norm_2(v);
        }
        return norm_l(v, l);
    }
}

norm_l :: inline (v: $V/VectorType, l: float64) -> float64 {
    return norm_l_default(v, l);
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


norm_2 :: inline (v: $V/VectorType) -> float64 {
    return norm_2_default(v);
}
#scope_file
norm_2_default :: (v: $V/VectorType) -> float64 {
    res : float64 = 0.0;
    for v {
        res += abs_sq(it);
    }
    return sqrt(res, float64);
}
#scope_export

norm_1 :: inline (v: $V/VectorType) -> float64 {
    return norm_1_default(v);
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

norm_inf :: inline (v: $V/VectorType) -> float64 {
    return norm_inf_default(v);
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
    dot : float64 = 0.0;
    for a {
        dot += conjugate(it) * get(b, it_index);
    }
    return acos( abs(dot) / (norm_2(a) * norm_2(b)) );
}
#scope_export