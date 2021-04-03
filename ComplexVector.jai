


CVector64 :: struct {
    dim: int;
    data: [..]Complex64;
}

operator [] :: inline (v: CVector64, i: int) -> Complex64 {
    return v.data[i];
}

operator == :: inline (a: CVector64, b: CVector64) -> bool {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

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


cvec64 :: (vec: Vector64) -> CVector64 {
    res := cvec64(vec.dim);
    for vec.data {
        res.data[it_index] = c64(it);
    }
    return res;
}


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


// inclusive bounds
cvec64_subvector :: (v: CVector64, from: int, to: int) -> CVector64 {
    #if CHECKS {
        assert(to > from);
    }
    res := cvec64(to-from);
    for 0..to-from-1 {
        res.data[it] = v[it + from];
    }
    return res;
}






operator + :: inline (a: CVector64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

operator + :: inline (a: CVector64, b: Vector64) -> CVector64 #symmetric {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] + b[it];
    }
    return res;
}

operator - :: inline (a: CVector64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}
operator - :: inline (a: CVector64, b: Vector64) -> CVector64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := cvec64(a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] - b[it];
    }
    return res;
}
operator - :: inline (a: Vector64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

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
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := CZERO;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}
operator * :: inline (a: CVector64, b: Vector64) -> Complex64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := CZERO;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}
operator * :: inline (a: Vector64, b: CVector64) -> Complex64 {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    res := CZERO;
    for 0..a.dim-1 {
        res += a[it] * b[it];
    }
    return res;
}




permute :: (v: CVector64, order: ..int) -> CVector64 {
    #if CHECKS {
        assert(v.dim == order.count);
    }
    res := cvec64(v.dim);
    for 0..v.dim-1 {
        res.data[it] = v[order[it]];
    }

    return res;
}

swap :: (v: CVector64, i: int, j: int) {
    v.data[i], v.data[j] = swap(v.data[i], v.data[j]);
}



conjugate :: (v: CVector64) -> CVector64 {
    res := cvec64(v.dim);
    for v.data {
        res.data[it_index] = conjugate(it);
    }
    return res;
}






outer_product :: (a: CVector64, b: CVector64) -> CMatrix64 {
    res := cmat64(a.dim, b.dim);
    for i: 0..a.dim {
        for j: 0..b.dim {
            res[i][j] = a[i] * b[j];
        }
    }
    return res;
}

outer_product :: (a: Vector64, b: CVector64) -> CMatrix64 {
    res := cmat64(a.dim, b.dim);
    for i: 0..a.dim {
        for j: 0..b.dim {
            res[i][j] = a[i] * b[j];
        }
    }
    return res;
}

outer_product :: (a: CVector64, b: Vector64) -> CMatrix64 {
    res := cmat64(a.dim, b.dim);
    for i: 0..a.dim {
        for j: 0..b.dim {
            res[i][j] = a[i] * b[j];
        }
    }
    return res;
}


operator / :: inline (a: CVector64, b: CVector64) -> CMatrix64 {
    return outer_product(a, b);
}
operator / :: inline (a: Vector64, b: CVector64) -> CMatrix64 {
    return outer_product(a, b);
}
operator / :: inline (a: CVector64, b: Vector64) -> CMatrix64 {
    return outer_product(a, b);
}


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.10 Apply Elementary Reflector, page 170
reflect :: (vec: CVector64, normal: CVector64) -> CVector64 {
    tau := 2.0 / (conjugate(normal) * normal);
    alpha := tau * conjugate(normal) * vec;
    return vec - normal * alpha;
}




// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Definition 3.5.2 Norms, page 171
norm :: (v: CVector64, n: float64) -> float64 {
    if n == Math.FLOAT64_INFINITY {
        m := abs(v[0]);
        for 1..v.dim-1 {
            tmp := abs(v[it]);
            m = ifx m < tmp then tmp else m;
        }
        return m;
    } else if n == 1.0 {
        res : float64 = 0.0;
        for 0..v.dim-1 {
            res := abs(v[it]);
        }
        return res;
    } else if n == 2.0 {
        return sqrt( (conjugate(v) * v).real, float64 );
    } else {
        res : float64 = 0.0;
        for 0..v.dim-1 {
            res += pow(abs(v[it]), n, float64);
        }
        return pow(res, 1.0/n, float64);
    }
}

norm_2 :: inline (v: CVector64) -> float64 {
    return norm(v, 2.0);
}

norm_1 :: inline (v: CVector64) -> float64 {
    return norm(v, 1.0);
}

norm_inf :: inline (v: CVector64) -> float64 {
    return norm(v, Math.FLOAT64_INFINITY);
}




cross :: (a: CVector64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(a.dim == 3 && b.dim == 3);
    }
    return cvec64( 
        a[1]*b[2] - a[2]*b[1], 
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}


angle :: (a: CVector64, b: CVector64) -> float64 {
    return acos( abs(conjugate(a) * b) / (norm_2(a) * norm_2(b)) );
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