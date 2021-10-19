

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

equals :: (a: $A/VectorType, b: $B/VectorType) -> bool {
    #if A.flags == .DENSE && B.flags == .DENSE {
        return inline equals_dense_dense(a,b);
    } else {
        // * default implementation
        return inline equals_default(a,b);
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

operator * :: (a: $A/VectorType, b: $B) -> VectorCastHelper(A,B).T #modify {
    if !is_number(B) {B = null; return;}    
} #symmetric {
    res := make(VectorCastHelper(A,B).T, false);
    for a {
        set(*res, it_index, it * b);
    }
    return res;
}

mul :: (a: DenseVector($T, $N), b: $B) -> DenseVector(UpCast(T,B).T, N) #symmetric {
    res := make(DenseVector(UpCast(T,B).T, N));
    mul(*res, b);
    return res;
}
mul :: (a: DenseHeapVector($T), b: $B) -> DenseHeapVector(UpCast(T,B).T) #symmetric {
    res := make(DenseHeapVector(UpCast(T,B).T));
    mul(*res, b);
    return res;
}
mul :: (a: VectorView(DenseVector($T, $N), $M), b: $B) -> DenseVector(UpCast(T,B).T, M) #symmetric {
    res : DenseVector(UpCast(T,B).T), M);
    init(*res);
    mul(*res, b);
    return res;
}
mul :: (a: VectorView(DenseHeapVector($T), $M), b: $B) -> DenseVector(UpCast(T,B).T, M) #symmetric {
    res : DenseVector(UpCast(T,B).T), M);
    init(*res);
    mul(*res, b);
    return res;
}
mul :: (a: HeapVectorView(DenseVector($T, $N)), b: $B) -> DenseHeapVector(UpCast(T,B).T) #symmetric {
    res : DenseHeapVector(UpCast(T,B).T));
    res.dim = a.dim;
    init(*res);
    mul(*res, b);
    return res;
}
mul :: (a: HeapVectorView(DenseHeapVector($T)), b: $B) -> DenseHeapVector(UpCast(T,B).T) #symmetric {
    res : DenseHeapVector(UpCast(T,B).T));
    res.dim = a.dim;
    init(*res);
    mul(*res, b);
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