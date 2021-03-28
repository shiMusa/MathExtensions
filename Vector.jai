


Vector64 :: struct {
    dim: int;
    data: [..]float64;
};


operator [] :: inline (v: Vector64, i: int) -> float64 {
    return v.data[i];
}


set :: inline (v: Vector64, i: int, val: float64) {
    v.data[i] = val;
}



vec64 :: (values: ..float64) -> Vector64 {
    res : Vector64;
    res.dim = values.count;
    array_resize(*res.data, values.count, false);
    for values {
        res.data[it_index] = it;
    }
    return res;
}

vec64 :: inline (dim: int) -> Vector64 {
    res : Vector64;
    res.dim = dim;
    array_resize(*res.data, dim, false);
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

copy :: inline (v: Vector64) -> Vector64 {
    res := vec64(v.dim);
    array_copy(*res.data, v.data);
    return res;
}






operator + :: inline (a: Vector64, b: Vector64) -> Vector64 {
    res := vec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

operator - :: inline (a: Vector64, b: Vector64) -> Vector64 {
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
    res : float64 = 0.0;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
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