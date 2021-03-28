


Vector64 :: [..]float64;



vec64 :: (values: ..float64) -> Vector64 {
    res : Vector64;
    array_resize(*res, values.count, false);
    for values {
        res[it_index] = it;
    }
    return res;
}

vec64 :: inline (dim: int) -> Vector64 {
    res : Vector64;
    array_resize(*res, dim, false);
    return res;
}

vec64_zeros :: inline (dim: int) -> Vector64 {
    res := vec64(dim);
    for 0..dim-1 {
        res[it] = 0.0;
    }
    return res;
}

vec64_ones :: inline (dim: int) -> Vector64 {
    res := vec64(dim);
    for 0..dim-1 {
        res[it] = 1.0;
    }
    return res;
}

vec64_basis :: inline (dim: int, el: int) -> Vector64 {
    res := vec64_zeros(dim);
    res[el] = 1.0;
    return res;
}

copy :: inline (v: Vector64) -> Vector64 {
    res := vec64(v.count);
    array_copy(*res, v);
    return res;
}






operator + :: inline (a: Vector64, b: Vector64) -> Vector64 {
    res := vec64(a.count);
    for a.count-1 {
        res[it] = a[it] + b[it];
    }
    return res;
}

operator - :: inline (a: Vector64, b: Vector64) -> Vector64 {
    res := vec64(a.count);
    for a.count-1 {
        res[it] = a[it] - b[it];
    }
    return res;
}

operator - :: inline (a: Vector64) -> Vector64 {
    res := vec64(a.count);
    for a.count-1 {
        res[it] = -a[it];
    }
    return res;
}

operator * :: inline (a: Vector64, val: float64) -> Vector64 #symmetric {
    res := vec64(a.count);
    for a.count-1 {
        res[it] = a[it] * val;
    }
    return res;
}

operator / :: inline (a: Vector64, val: float64) -> Vector64 {
    res := vec64(a.count);
    for a.count-1 {
        res[it] = a[it] / val;
    }
    return res;
}


operator * :: inline (a: Vector64, b: Vector64) -> float64 {
    res := 0.0;
    for a.count-1 {
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