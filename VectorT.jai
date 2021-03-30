

VectorT :: struct (T: Type) {
    type: Type = T;
    dim: int;
    data: [..]T;
}

operator [] :: inline (v: VectorT, i: int) -> v.type {
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
    res.type = type;
    res.dim = values.count;
    array_resize(*res.data, values.count, false);
    for values {
        res.data[it_index] = ncast(type, it);
    }
    return res;
}


test_vectort :: () {
    print("\n\nVectorT -------------------------------\n\n");
    v := vec(float64, 1.0, 2.0, 3.0);
    c := vec(Complex64, R, I, R+I);

    print("%, %\n", v, v.type);
    print("%, %\n", c, c.type);
}