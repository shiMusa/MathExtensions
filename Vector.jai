/*
 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##
*/


str :: (v: $V/VectorType, indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    print_to_builder(*builder, "\n%", indent);
    for i: 0..V.N-1 {
        #if T == Complex64 {
            print_to_builder(*builder, "%, ", str(get(v,i)) );
        } else {
            print_to_builder(*builder, "%, ", get(v,i) );
        }
    }
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


#scope_file
VectorCastHelper :: struct(V: Type, t: Type) {
    T :: vector_cast_helper(V, t);
}
vector_cast_helper :: ($V: Type, $t: Type) -> Type {
    return DenseVector(t, V.N); // TODO other types
}
#scope_export

vec_cast :: ($type: Type, v: $V/VectorType) -> VectorCastHelper(V, type).T {
    res : VectorCastHelper(V, type).T; // initialize all values to zero, important for sparse?

    for v {
        set_lin(*res, it_lin_index, ncast(type, it));
    }

    return res;
}

copy :: (v: $V/VectorType) -> V {
    res : V = ---;

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
        #assert(is_minor_number(V.data_type, T));
        assert(values.count <= V.N, "too many values (%) for vector of % entries\n", values.count, V.N);
    }
    for values {
        set(v, it_index, ncast(V.data_type, it));
    }
}

vec_init_all :: inline (v: *$V/VectorType, value: $T) {
    #if CHECKS {
        #assert(is_minor_number(V.data_type, T));
    }
    for raw(v) {
        set_raw(v, it_raw_index, ncast(V.data_type, value));
    }
}



#scope_file
VectorSubvectorHeper :: struct(V: Type, N: int) {
    T :: #run vector_subvector_helper(V, N);
}
vector_subvector_helper :: ($V: Type, $N: int) -> Type {
    return DenseVector(V.data_type, N);
}
#scope_export
// inclusive left, exclusive right bounds
vec_subvector :: (v: $V/VectorType, $from: int, $to: int) -> VectorSubvectorHeper(V, to-from).T {
    #if CHECKS {
        #run assert(to > from, "bounds not possible: to(%) <= from(%)\n", to, from);
        #run assert(to < N, "bound not possible: to(%) >= N(%)\n", to, N);
        #run assert(from >= 0, "bound not possible: from(%) < 0\n", from);
    }
    res : VectorSubvectorHeper(V, to-from).T = ---;

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



operator == :: (a: Vector($A, $N), b: Vector($B, $M)) -> bool {
    #if N != M {
        return false;
    } else {
        for a.data {
            if it != b[it_index] {
                return false;
            }
        }
        return true;
    }
}




operator * :: (a: Vector($A, $N), b: Vector($B, N)) -> UpCast(A,B).T {
    res := zero(UpCast(A,B).T);
    for 0..N-1 {
        res += a[it] * b[it];
    }
    return res;
}

operator * :: (a: Vector($A, $N), b: Scalar($B)) -> UpCast(A,B).T #symmetric {
    res : Vector(UpCast(A,B).T, N) = ---;
    for 0..N-1 {
        res.data[it] = a[it] * b;
    }
    return res;
}



mul :: inline (a: *Vector($A, $N), b: Scalar($B)) {
    #if CHECKS {
        #assert(is_minor_number(A,B));
    }
    for 0..N-1 {
        a.data[it] *= b.value;
    }
}









operator / :: inline (a: Vector($A, $N), b: Vector($B, $M)) -> Matrix(UpCast(A,B).T, N, M) {
    return outer_product(a, b);
}

operator / :: (a: Vector($A, $N), b: Scalar($B)) -> Vector(UpCast(A,B).T, N) {
    res : Vector(UpCast(A,B).T, N) = ---;
    for 0..N-1 {
        res.data[it] = a[it] / b.value;
    }
    return res;
}

div :: inline (a: *Vector($A, $N), b: Scalar($B)) {
    #if CHECKS {
        #assert(is_minor_number(A,B));
    }
    for 0..N-1 {
        a.data[it] /= b.value;
    }
}

outer_product :: (a: Vector($A, $N), b: Vector($B, $M)) -> Matrix(UpCast(A,B).T, N, M) {
    res : Matrix(UpCast(A,B).T, N, M) = ---;
    for i: 0..N-1 {
        for j: 0..M-1 {
            res[i][j] = a[i] * b[j];
        }
    }
    return res;
}






operator + :: (a: Vector($A, $N), b: Vector($B, N)) -> Vector(UpCast(A,B).T, N) {
    res : Vector(UpCast(A,B).T, N) = ---;
    for 0..N-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

add :: (a: *Vector($A, $N), b: Vector($B, N)) {
    #if CHECKS {
        #assert(is_minor_number(A, B));
    }
    for 0..N-1 {
        a.data[it] += b[it];
    }
}




operator - :: (a: Vector($A, $N), b: Vector($B, N)) -> Vector(UpCast(A,B).T, N) {
    res : Vector(UpCast(A,B).T, N) = ---;
    for 0..N-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}

sub :: (a: *Vector($A, $N), b: Vector($B, N)) {
    #if CHECKS {
        #assert(is_minor_number(A, B));
    }
    for 0..N-1 {
        a.data[it] -= b[it];
    }
}




operator - :: (a: Vector($A, $N)) -> Vector(A, N) {
    res : Vector(C, N) = ---;
    for 0..N-1 {
        res.data[it] = -a[it];
    }
    return res;
}

neg :: (a: *Vector($A, $N)) {
    for 0..N-1 {
        a.data[it] = -a.data[it];
    }
}









cross :: (a: Vector($A, $N), b: Vector($B, N)) -> Vector(UpCast(A,B).T, N) {
    #assert(N == 3);
    res : Vector(UpCast(A,B).T, N) = ---;
    res.data[0] = a[1]*b[2] - a[2]*b[1];
    res.data[1] = a[2]*b[0] - a[0]*b[2];
    res.data[2] = a[0]*b[1] - a[1]*b[0];
    return res;
}


reflect :: (a: Vector($A, $N), b: Vector($B, N)) -> Vector(UpCast(A,B).T, N) {
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

// Matrix Computations, 4th edition
// The Johns Hopkins University Press
// Algorithm 1.1.2 Saxpy, page 4
saxpy :: (a: *Vector($A, $N), b: Vector($B, N), factor: $C) {
    #if CHECKS {
        #assert(is_minor_number(A, B) && is_minor_number(A, C));
    }

    for 0..N-1 {
        a.data[it] += factor * b.data[it];
    }
}



// Matrix Computations, 4th edition
// The Johns Hopkins University Press
// Algorithm 1.1.3 Row-Oriented Gaxpy, page 5
gaxpy :: (y: *Vector($T, $N), A: Matrix($U, N, M), x: Vector($V, $M)) {
    #if CHECKS {
        #assert(is_minor_number(T, U) && is_minor_number(T, V));
    }

    for i: 0..N {
        for j: 0..M {
            y.data[i] += A[i][j] * x.data[j];
        }
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

permute :: (v: Vector($V, $N), order: ..int) -> Vector(V, N) {
    #if CHECKS {
        assert(N == order.count);
    }
    res : Vector(V, N) = ---;
    for 0..N-1 {
        res.data[it] = v[order[it]];
    }

    return res;
}

permute :: (v: *Vector($V, $N), order: ..int) {
    #if CHECKS {
        assert(N == order.count);
    }

    for 0..N-1 {
        if order[it] <= it then continue;

        v.data[it], v.data[order[it]] = swap(v.data[it], v.data[order[it]]);
    }
}

swap :: (v: *Vector($V, $N), i: int, j: int) {
    #if CHECKS {
        assert(i >= 0 && i < N);
        assert(j >= 0 && j < N);
    }
    v.data[i], v.data[j] = swap(v.data[i], v.data[j]);
}



conjugate :: (v: Vector($V, $N)) -> Vector(V, N) {
    res : Vector(V, N) = ---;
    for v.data {
        res.data[it_index] = conjugate(it);
    }
    return res;
}
conjugate :: (v: *Vector($V, $N)) {
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
norm :: (v: Vector($T, $N), $$n: float64) -> float64 {

    #if is_constant(n) {
        #if n == Math.FLOAT64_INFINITY {
            m := abs(v[0]);
            for 1..N-1 {
                tmp := abs(v[it]);
                m = ifx m < tmp then tmp else m;
            }
            return m;
        } else #if n == 1.0 {
            res : float64 = 0.0;
            for 0..N-1 {
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
_generic_norm :: inline (v: Vector($T, $N), n: float64) -> float64 {
    res : float64 = 0.0;
    for 0..N-1 {
        res += pow(abs(v[it]), n, float64);
    }
    return pow(res, 1.0/n, float64);
}
#scope_export


norm_2 :: inline (v: Vector($T, $N)) -> float64 {
    return norm(v, 2.0);
}

norm_1 :: inline (v: Vector($T, $N)) -> float64 {
    return norm(v, 1.0);
}

norm_inf :: inline (v: Vector($T, $N)) -> float64 {
    return norm(v, Math.FLOAT64_INFINITY);
}


angle :: (a: Vector($A, $N), b: Vector($B, N)) -> float64 {
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
        v := vec(5, float64, 1.0, 2.0, 3.0, 4.0, 5.0);
        print("%, type_of(v) = %, Vector(T) = Vector(%)\n", v, type_of(v), v.T);
        a := v + v;
        print("%, type_of(a) = %, Vector(T) = Vector(%)\n", a, type_of(a), a.T);
        a = a + v;
        print("%, type_of(a) = %, Vector(T) = Vector(%)\n", a, type_of(a), a.T);
    }
}

