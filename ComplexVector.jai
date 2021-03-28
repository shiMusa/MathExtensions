


CVector64 :: struct {
    dim: int;
    data: [..]Complex64;
}

operator [] :: inline (v: CVector64, i: int) -> Complex64 {
    return v.data[i];
}

operator == :: inline (a: CVector64, b: CVector64) -> bool {
    r := true;
    for a.data {
        r = r && it == b.data[it_index];
    }
    return a.dim == b.dim && r;
}


cvec64 :: (values: ..Complex64) -> CVector64 {
    res : CVector64;
    res.dim = values.count;
    array_resize(*res.data, values.count, false);
    for values {
        res.data[it_index] = it;
    }
    return res;
}

// TODO: wait for bug fix in compiler
// cvec64 :: (values: ..float64) -> CVector64 {
//     res : CVector64;
//     res.dim = values.count;
//     array_resize(*res.data, values.count, false);
//     for values {
//         res.data[it_index] = c64(it, 0.0);
//     }
//     return res;
// }

cvec64 :: inline (dim: int) -> CVector64 {
    res : CVector64;
    res.dim = dim;
    array_resize(*res.data, dim, false);
    return res;
}

cvec64_zeros :: inline (dim: int) -> CVector64 {
    res := cvec64(dim);
    for 0..dim-1 {
        res.data[it] = c64(0.0, 0.0);
    }
    return res;
}

cvec64_ones :: inline (dim: int) -> CVector64 {
    res := cvec64(dim);
    for 0..dim-1 {
        res.data[it] = R;
    }
    return res;
}

cvec64_basis :: inline (dim: int, el: int) -> CVector64 {
    res := cvec64_zeros(dim);
    res.data[el] = R;
    return res;
}

copy :: inline (v: CVector64) -> CVector64 {
    res := cvec64(v.dim);
    array_copy(*res.data, v.data);
    return res;
}






operator + :: inline (a: CVector64, b: CVector64) -> CVector64 {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

operator + :: inline (a: CVector64, b: Vector64) -> CVector64 #symmetric {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

operator - :: inline (a: CVector64, b: CVector64) -> CVector64 {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}
operator - :: inline (a: CVector64, b: Vector64) -> CVector64 {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}
operator - :: inline (a: Vector64, b: CVector64) -> CVector64 {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}

operator - :: inline (a: CVector64) -> CVector64 {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = -a[it];
    }
    return res;
}

operator * :: inline (a: CVector64, val: Complex64) -> CVector64 #symmetric {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] * val;
    }
    return res;
}
operator * :: inline (a: CVector64, val: float64) -> CVector64 #symmetric {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] * val;
    }
    return res;
}

operator / :: inline (a: CVector64, val: Complex64) -> CVector64 {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] / val;
    }
    return res;
}
operator / :: inline (a: CVector64, val: float64) -> CVector64 {
    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] / val;
    }
    return res;
}


operator * :: inline (a: CVector64, b: CVector64) -> Complex64 {
    res := CZERO;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}
operator * :: inline (a: CVector64, b: Vector64) -> Complex64 {
    res := CZERO;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}
operator * :: inline (a: Vector64, b: CVector64) -> Complex64 {
    res := CZERO;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}



conjugate :: (v: CVector64) -> CVector64 {
    res := cvec64(v.dim);
    for v.data {
        res.data[it_index] = conjugate(it);
    }
    return res;
}










test_complex_vector :: () {
    print("\n\nComplex Vector -------------------------------\n\n");

    {
        v := cvec64(R, c64(2.0), c64(3.0));
        print("%\n", v);
    }
    {
        v := cvec64(R, 2.0*R, 3.0*R);
        print("%\n", v);
    }
    // {
    //     v := cvec64(R, 2.0, 3.0);
    //     print("%\n", v);
    // }

    {
        v := cvec64_basis(5, 3);
        w := copy(v);
        assert(w[3] == R);
        print("%\n", w);
    }
}