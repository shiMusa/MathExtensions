


Vector64 :: struct {
    dim: int;
    data: []float64;
};


operator [] :: inline (v: Vector64, i: int) -> float64 {
    return v.data[i];
}


copy :: inline (v: Vector64) -> Vector64 {
    res : Vector64 = ---;
    res.dim = v.dim;
    res.data = array_copy(v.data);
    return res;
}


vec64 :: (values: ..float64) -> Vector64 {
    res : Vector64 = ---;
    res.dim = values.count;
    res.data = array_copy(values);
    return res;
}

vec64 :: inline (dim: int) -> Vector64 {
    res : Vector64 = ---;
    res.dim = dim;
    res.data = NewArray(dim, float64, initialized=false);
    return res;
}

vec64_zeros :: inline (dim: int) -> Vector64 {
    res := vec64(dim);
    for 0..dim-1 {
        res.data[it] = 0.0;
    }
    return res;
}

vec64_ones :: inline (dim: int) -> Vector64 {
    res := vec64(dim);
    for 0..dim-1 {
        res.data[it] = 1.0;
    }
    return res;
}

vec64_basis :: inline (dim: int, el: int) -> Vector64 {
    res := vec64_zeros(dim);
    res.data[el] = 1.0;
    return res;
}








operator + :: inline (a: Vector64, b: Vector64) -> Vector64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

operator - :: inline (a: Vector64, b: Vector64) -> Vector64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}

operator - :: inline (a: Vector64) -> Vector64 {
    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = -a[it];
    }
    return res;
}

operator * :: inline (a: Vector64, val: float64) -> Vector64 #symmetric {
    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] * val;
    }
    return res;
}

operator / :: inline (a: Vector64, val: float64) -> Vector64 {
    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] / val;
    }
    return res;
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




norm_2 :: (v: Vector64) -> float64 {
    return v * v;
}




outer_product :: (a: Vector64, b: Vector64) -> Matrix64 {
    res := mat64(a.dim, b.dim);
    for i: 0..a.dim {
        for j: 0..b.dim {
            res[i][j] = a[i] * b[j];
        }
    }
    return res;
}


operator / :: inline (a: Vector64, b: Vector64) -> Matrix64 {
    return outer_product(a, b);
}






test_vector :: () {
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