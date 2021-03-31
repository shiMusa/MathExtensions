

VectorT :: struct (T: Type) {
    dim: int;
    data: [..]T;
}

operator [] :: inline (v: VectorT($T), i: int) -> T {
    return v.data[i];
}

operator == :: inline (a: VectorT, b: VectorT) -> bool {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    for a.data {
        if it != b.data[it_index] then return false;
    }
    return true;
}



vec :: ($type: Type, values: ..$T) -> VectorT(type) {
    res : VectorT(type);
    res.dim = values.count;
    array_resize(*res.data, values.count, false);
    for values {
        res.data[it_index] = ncast(type, it);
    }
    return res;
}

_vec64 :: inline (values: ..float64) -> VectorT(float64) {
    return vec(float64, ..values);
}
_cvec64 :: inline (values: ..Complex64) -> VectorT(Complex64) {
    return vec(Complex64, ..values);
}

vec :: ($type: Type, dim: int) -> VectorT(type) {
    res: VectorT(type);
    res.dim = dim;
    array_resize(*res.data, dim, false);
    return res;
}

vec_zeros :: ($type: Type, dim: int) -> VectorT(type) {
    res: VectorT(type);
    res.dim = dim;
    array_resize(*res.data, dim, false);
    for 0..dim-1 {
        res.data[it] = zero(type);
    }
    return res;
}

vec_ones :: ($type: Type, dim: int) -> VectorT(type) {
    res: VectorT(type);
    res.dim = dim;
    array_resize(*res.data, dim, false);
    for 0..dim-1 {
        res.data[it] = one(type);
    }
    return res;
}

vec_basis :: inline ($type: Type, dim: int, el: int) -> VectorT(type) {
    res := vec_zeros(type, dim);
    res.data[el] = one(type);
    return res;
}

copy :: inline(v: VectorT($T)) -> VectorT(T) {
    res := vec(T, v.dim);
    array_copy(*res.data, v.data);
    return res;
}

vec :: ($type: Type, v: VectorT) -> VectorT(type) {
    res := vec(type, v.dim);
    for 0..res.dim {
        res.data[it] = ncast(type, v[it]);
    }
    return res;
}





operator + :: inline (a: VectorT($A), b: VectorT($B)) -> VectorT(Super(A,B).T) {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    type :: #run super_type(A, B);
    res := vec(type, a.dim);
    for 0..a.dim-1 {
        res.data[it] = ncast(type, a[it]) + ncast(type, b[it]);
    }
    return res;
}

operator - :: inline (a: VectorT($A), b: VectorT($B)) -> VectorT(Super(A,B).T) {
    #if CHECKS {
        assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
    }

    type :: #run super_type(A, B);
    res := vec(type, a.dim);
    for 0..a.dim-1 {
        res.data[it] = ncast(type, a[it]) - ncast(type, b[it]);
    }
    return res;
}

operator - :: inline (a: VectorT($A)) -> VectorT(A) {
    res := vec(A, a.dim);
    for 0..a.dim-1 {
        res.data[it] = -a[it];
    }
    return res;
}

operator * :: inline (a: VectorT($A), val: $T) -> VectorT(Super(A,T).T) #symmetric {
    type :: #run super_type(A, T);

    res := vec(type, a.dim);
    for 0..a.dim-1 {
        res.data[it] = a[it] * val;
    }
    return res;
}

operator / :: inline (a: VectorT($A), val: $T) -> VectorT(Super(A,T).T) {
    return a * (1.0/val);
}


// causes compiler crash!

// operator * :: inline (a: VectorT($A), b: VectorT($B)) -> Super(A,B).T {
//     #if CHECKS {
//         assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
//     }
//     type :: #run super_type(A, B);

//     res := zero(type);
//     for 0..a.dim-1 {
//         res += a[it] * b[it];
//     }
//     return res;
// }





test_vectort :: () {
    print("\n\nVectorT -------------------------------\n\n");
    {
        v := _vec64(1.0, 2.0, 3.0);
        c := _cvec64(R, I, R+I);

        print("%\n", v);
        print("%\n", c);
    }

    {
        v := vec_basis(Complex64, 3, 1);
        print("%\n", v);
    }

    {
        v := _vec64(1.0, 2.0, 3.0);
        c := _cvec64(R, I, R+I);
        cpy := copy(c);
        assert( cpy == c );
        s := v + c;
    }

    {
        v := _vec64(1.0, 2.0, 3.0);
        r := v * I;
        s := v / I;
        assert( _cvec64(I, 2.0*I, 3.0*I) == r);
        assert( _cvec64(-I, -2.0*I, -3.0*I) == s);
    }
}