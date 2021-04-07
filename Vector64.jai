




str :: (v: Vector64, indent: string = "    ") -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    print_to_builder(*builder, "\n%", indent);
    for i: 0..v.dim-1 {
        print_to_builder(*builder, "%, ", v[i] );
    }
    return builder_to_string(*builder);
}








Vector64 :: struct {
    dim: int;
    data: []float64;
};


operator [] :: inline (v: Vector64, i: int) -> float64 {
    return v.data[i];
}


copy :: inline (v: Vector64) -> Vector64 #must {
    res : Vector64 = ---;
    res.dim = v.dim;
    res.data = array_copy(v.data);
    return res;
}


vec64 :: (values: ..float64) -> Vector64 #must {
    res : Vector64 = ---;
    res.dim = values.count;
    res.data = array_copy(values);
    return res;
}

vec64 :: inline (dim: int) -> Vector64 #must {
    res : Vector64 = ---;
    res.dim = dim;
    res.data = NewArray(dim, float64, initialized=false);
    return res;
}

vec64_zeros :: inline (dim: int) -> Vector64 #must {
    res := vec64(dim);
    for 0..dim-1 {
        res.data[it] = 0.0;
    }
    return res;
}

vec64_ones :: inline (dim: int) -> Vector64 #must {
    res := vec64(dim);
    for 0..dim-1 {
        res.data[it] = 1.0;
    }
    return res;
}

vec64_basis :: inline (dim: int, el: int) -> Vector64 #must {
    res := vec64_zeros(dim);
    res.data[el] = 1.0;
    return res;
}








// inclusive bounds
vec64_subvector :: (v: Vector64, from: int, to: int) -> Vector64 #must {
    #if CHECKS {
        assert(to > from);
    }
    res := vec64(to-from);
    for 0..to-from-1 {
        res.data[it] = v[it + from];
    }
    return res;
}









operator + :: inline (a: Vector64, b: Vector64) -> Vector64 #must {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}
add :: (a: Vector64, b: Vector64) {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    for 0..a.dim-1 {
        a.data[it] += b[it];
    }
}

operator - :: inline (a: Vector64, b: Vector64) -> Vector64 #must {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}
sub :: (a: Vector64, b: Vector64) {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    for 0..a.dim-1 {
        a.data[it] -= b[it];
    }
}

operator - :: inline (a: Vector64) -> Vector64 #must {
    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = -a[it];
    }
    return res;
}
neg :: (a: Vector64) {
    for 0..a.dim-1 {
        a.data[it] = -a.data[it];
    }
}

operator * :: inline (a: Vector64, val: float64) -> Vector64 #must #symmetric {
    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] * val;
    }
    return res;
}
mul :: (a: Vector64, val: float64) {
    for 0..a.dim-1 {
        a.data[it] *= val;
    }
}

operator / :: inline (a: Vector64, val: float64) -> Vector64 #must {
    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] / val;
    }
    return res;
}
div :: (a: Vector64, val: float64) {
    for 0..a.dim-1 {
        a.data[it] /= val;
    }
}


operator * :: inline (a: Vector64, b: Vector64) -> float64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res : float64 = 0.0;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}







permute :: (v: Vector64, order: ..int) -> Vector64 #must {
    #if CHECKS {
        assert(v.dim == order.count);
    }
    res := vec64(v.dim);
    for 0..v.dim-1 {
        res.data[it] = v[order[it]];
    }

    return res;
}

permute :: (v: *Vector64, order: ..int) {
    #if CHECKS {
        assert(v.dim == order.count);
    }

    for 0..v.dim-1 {
        if order[it] <= it then continue;

        v.data[it], v.data[order[it]] = swap(v.data[it], v.data[order[it]]);
    }
}

swap :: (v: Vector64, i: int, j: int) {
    v.data[i], v.data[j] = swap(v.data[i], v.data[j]);
}




// for generic functions
conjugate :: inline(v: Vector64) -> Vector64 #must {
    return copy(v);
}
conjugate :: inline(v: *Vector64) {}









outer_product :: (a: Vector64, b: Vector64) -> Matrix64 #must {
    res := mat64(a.dim, b.dim);
    for i: 0..a.dim {
        for j: 0..b.dim {
            res[i][j] = a[i] * b[j];
        }
    }
    return res;
}


operator / :: inline (a: Vector64, b: Vector64) -> Matrix64 #must {
    return outer_product(a, b);
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.10 Apply Elementary Reflector, page 170
reflect :: (vec: Vector64, normal: Vector64) -> Vector64 #must {
    tau := 2.0 / norm_2(normal);
    alpha := tau * normal * vec;
    return vec - normal * alpha;
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Definition 3.5.2 Norms, page 171
norm :: (v: Vector64, $$n: float64) -> float64 {

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
            return Math.sqrt( v * v );
        }
    } else {
        res : float64 = 0.0;
        for 0..v.dim-1 {
            res += Math.pow(abs(v[it]), n);
        }
        return Math.pow(res, 1.0/n);
    }
}

norm_2 :: inline (v: Vector64) -> float64 {
    return norm(v, 2.0);
}

norm_1 :: inline (v: Vector64) -> float64 {
    return norm(v, 1.0);
}

norm_inf :: inline (v: Vector64) -> float64 {
    return norm(v, Math.FLOAT64_INFINITY);
}





cross :: (a: Vector64, b: Vector64) -> Vector64 #must {
    #if CHECKS {
        assert(a.dim == 3 && b.dim == 3);
    }
    return vec64( 
        a[1]*b[2] - a[2]*b[1], 
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}


angle :: (a: Vector64, b: Vector64) -> float64 {
    return acos( abs(a * b) / (norm_2(a) * norm_2(b)) );
}







test_real_vector :: () {
    print("\n\nVector -------------------------------\n\n");

    {
        v := vec64(1.0, 2.0, 3.0);
        print("%\n", v);
    }

    {
        v := vec64_basis(5, 3);
        w := copy(v);
        assert(w[3] == 1.0);
        print("%\n", w);
    }
}