



/*
 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##
*/


str :: (v: Vector($T), indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    print_to_builder(*builder, "\n%", indent);
    for i: 0..v.dim-1 {
        #if T == Complex64 {
            print_to_builder(*builder, "%, ", str(v[i]) );
        } else {
            print_to_builder(*builder, "%, ", v[i] );
        }
    }
    return builder_to_string(*builder);
}




/*
 ######  ######## ########  ##     ##  ######  ########
##    ##    ##    ##     ## ##     ## ##    ##    ##
##          ##    ##     ## ##     ## ##          ##
 ######     ##    ########  ##     ## ##          ##
      ##    ##    ##   ##   ##     ## ##          ##
##    ##    ##    ##    ##  ##     ## ##    ##    ##
 ######     ##    ##     ##  #######   ######     ##
*/

Vector :: struct (T: Type) {
    dim: int;
    data: []T;
}

operator [] :: inline (v: Vector($T), i: int) -> T {
    #if CHECKS {
        assert(i < v.dim, "index out of dimensions: i = %, dim = %\n", i, v.dim);
    }
    return v.data[i];
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


vec :: ($type: Type, values: ..$T) -> Vector(type) {
    res : Vector(type) = ---;
    res.dim = values.count;
    res.data = NewArray(values.count, type);
    for values {
        res.data[it_index] = ncast(type, it);
    }
    return res;
}


copy :: (v: Vector($T)) -> Vector(T) {
    res : Vector(T) = ---;
    res.dim = v.dim;
    res.data = array_copy(v.data);
    return res;
}

vec :: ($type: Type, dim: int) -> Vector(type) {
    res: Vector(type) = ---;
    res.dim = dim;
    res.data = NewArray(dim, type, false);
    return res;
}

vec_zeros :: ($type: Type, dim: int) -> Vector(type) {
    res: Vector(type) = ---;
    res.dim = dim;
    res.data = NewArray(dim, type, false);
    for 0..dim-1 {
        res.data[it] = zero(type);
    }
    return res;
}

vec_ones :: ($type: Type, dim: int) -> Vector(type) {
    res: Vector(type);
    res.dim = dim;
    res.data = NewArray(dim, type, false);
    for 0..dim-1 {
        res.data[it] = one(type);
    }
    return res;
}

vec_basis :: inline ($type: Type, dim: int, el: int) -> Vector(type) {
    res := vec_zeros(type, dim);
    res.data[el] = one(type);
    return res;
}

vec :: ($type: Type, v: Vector) -> Vector(type) {
    res := vec(type, v.dim);
    for 0..res.dim-1 {
        res.data[it] = ncast(type, v[it]);
    }
    return res;
}

// inclusive left, exclusive right bounds
vec_subvector :: (v: Vector($T), from: int, to: int) -> Vector(T) {
    #if CHECKS {
        assert(to > from);
    }
    res := vec(T, to-from);
    for 0..to-from-1 {
        res.data[it] = v[it + from];
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



operator == :: (a: Vector($A), b: Vector($B)) -> bool {
    if a.dim != b.dim {
        return false;
    }
    for a.data {
        if it != b.data[it_index] {
            return false;
        }
    }
    return true;
}




operator * :: (a: Vector($A), b: Vector($B)) -> $C #modify {
    C = dominant_type(A, B);
}{
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := zero(C);
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}

operator * :: (a: Vector($A), b: $B) -> $C #modify {
    if !is_scalar(B) { C = null; return; }
    C = dominant_type(A, B);
} #symmetric {
    res := vec(C, a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] * b;
    }
    return res;
}



mul :: inline (a: *Vector($A), val: $B) {
    #if CHECKS {
        assert(is_minor_type(A,B));
    }
    for 0..a.dim-1 {
        a.data[it] *= val;
    }
}









operator / :: inline (a: Vector($A), b: Vector($B)) -> Matrix($C) #modify {
    C = dominant_type(A, B);
}{
    return outer_product(a, b);
}

operator / :: (a: Vector($A), b: $B) -> Vector($C) #modify {
    if !is_scalar(B) { C = null; return; }
    C = dominant_type(A, B);
} {
    res := vec(C, a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] / b;
    }
    return res;
}

div :: inline (a: *Vector($A), val: $B) {
    #if CHECKS {
        assert(is_minor_type(A,B));
    }
    for 0..a.dim-1 {
        a.data[it] /= val;
    }
}

outer_product :: (a: Vector($A), b: Vector($B)) -> Vector($C) #modify {
    C = dominant_type(A, B);
} {
    res := mat(C, a.dim, b.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.dim-1 {
            res[i][j] = a[i] * b[j];
        }
    }
    return res;
}






operator + :: (a: Vector($A), b: Vector($B)) -> Vector($C) #modify {
    C = dominant_type(A, B);
} {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }
    res := vec(C, a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

add :: (a: *Vector($A), b: Vector($B)) {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
        assert(is_minor_type(A, B));
    }
    for 0..a.dim-1 {
        a.data[it] += b[it];
    }
}




operator - :: (a: Vector($A), b: Vector($B)) -> Vector($C) #modify {
    C = dominant_type(A, B);
} {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }
    res := vec(C, a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}

sub :: (a: *Vector($A), b: Vector($B)) {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
        assert(is_minor_type(A, B));
    }
    for 0..a.dim-1 {
        a.data[it] -= b[it];
    }
}




operator - :: (a: Vector($A)) -> Vector(A) {
    res := vec(C, a.dim);
    for 0..a.dim-1 {
        res.data[it] = -a[it];
    }
    return res;
}

neg :: (a: *Vector($A)) {
    for 0..a.dim-1 {
        a.data[it] = -a.data[it];
    }
}









cross :: (a: Vector($A), b: Vector($B)) -> Vector($C) #modify {
    C = dominant_type(A, B);
} {
    #if CHECKS {
        assert(a.dim == 3 && b.dim == 3);
    }
    return vec(C, 
        a[1]*b[2] - a[2]*b[1], 
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}


reflect :: (a: Vector($A), b: Vector($B)) -> Vector($C) #modify {
    C = dominant_type(A, B);
} {
    tau := 2.0 / (conjugate(normal) * normal);
    alpha := tau * conjugate(normal) * vec;
    return vec - normal * alpha;
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

permute :: (v: Vector($V), order: ..int) -> Vector(V) {
    #if CHECKS {
        assert(v.dim == order.count);
    }
    res := vec(V, v.dim);
    for 0..v.dim-1 {
        res.data[it] = v[order[it]];
    }

    return res;
}

permute :: (v: *Vector($V), order: ..int) {
    #if CHECKS {
        assert(v.dim == order.count);
    }

    for 0..v.dim-1 {
        if order[it] <= it then continue;

        v.data[it], v.data[order[it]] = swap(v.data[it], v.data[order[it]]);
    }
}

swap :: (v: *Vector($V), i: int, j: int) {
    v.data[i], v.data[j] = swap(v.data[i], v.data[j]);
}



conjugate :: (v: Vector($V)) -> Vector(V) {
    res := vec(V, v.dim);
    for v.data {
        res.data[it_index] = conjugate(it);
    }
    return res;
}
conjugate :: (v: *Vector($V)) {
    for v.data {
        v.data[it_index] = conjugate(it);
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

// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Definition 3.5.2 Norms, page 171
norm :: (v: Vector($T), $$n: float64) -> float64 {

    #if is_constant(n) {
        #if n == Math.FLOAT64_INFINITY {
            m := abs(v[0]);
            for 1..v.dim-1 {
                tmp := abs(v[it]);
                m = ifx m < tmp then tmp else m;
            }
            return m;
        } else #if n == 1.0 {
            res : float64 = 0.0;
            for 0..v.dim-1 {
                res := abs(v[it]);
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
_generic_norm :: inline (v: Vector($T), n: float64) -> float64 {
    res : float64 = 0.0;
    for 0..v.dim-1 {
        res += pow(abs(v[it]), n, float64);
    }
    return pow(res, 1.0/n, float64);
}
#scope_export


norm_2 :: inline (v: Vector($T)) -> float64 {
    return norm(v, 2.0);
}

norm_1 :: inline (v: Vector($T)) -> float64 {
    return norm(v, 1.0);
}

norm_inf :: inline (v: Vector($T)) -> float64 {
    return norm(v, Math.FLOAT64_INFINITY);
}


angle :: (a: Vector($A), b: Vector($B)) -> float64 {
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
    print("\n\nVector -------------------------------\n\n");
    {
        v := vec(float64, 1.0, 2.0, 3.0, 4.0, 5.0);
        print("%, type_of(v) = %, Vector(T) = Vector(%)\n", v, type_of(v), v.T);
        a := v + v;
        print("%, type_of(a) = %, Vector(T) = Vector(%)\n", a, type_of(a), a.T);
        a = a + v;
        print("%, type_of(a) = %, Vector(T) = Vector(%)\n", a, type_of(a), a.T);
    }
}

