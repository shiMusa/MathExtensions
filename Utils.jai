





calculation :: (code: Code) #expand {
    _pool: Pool.Pool;
    Pool.set_allocators(*_pool);
    
    {
        push_allocator(Pool.pool_allocator, *_pool);
        #insert code;
    }

    print("The calculation pool contains % bytes.\n", _pool.memblock_size - _pool.bytes_left);
    print("Releasing the pool now.\n");
    Pool.release(*_pool);
}




ncast :: inline ($T: Type, num: $N) -> T {
    #if T == Complex64 {
        #if N == Complex64 {
            return num;
        } else {
            return c64(cast(float64)num);
            // res : Complex64;
            // res.real = cast(float64) num;
            // res.imag = 0.0;
            // return res;
        }
    } else {
        #if N == Complex64 {
            return cast(T) num.real;
        } else
            return cast(T) num;
    }
}

zero :: inline ($T: Type) -> T {
    return ncast(T, 0);
}

one :: inline ($T: Type) -> T {
    return ncast(T, 1);
}



DELTA :: 1.0e-6;
approx :: (a: $T, b: T) -> bool {
    // #if T == CVector64 || T == Vector64 {
    //     return norm_2(a - b) < DELTA;
    // } else #if T == CMatrix64 || T == Matrix64 {
    //     return norm_frobenius(a - b) < DELTA;
    // } else
        return abs(a - b) < DELTA;
}




#scope_module




array_view_range :: inline (arr: []$T, from: s64, to: s64) -> []T {
    #if CHECKS {
        assert(to < arr.count, "Index out of bounds: to = %, arr.count = %", to, arr.count);
        assert(from < to, "Negative range: from = %, to = %", from, to);
    }
    res := arr;
    res.data += from;
    res.count = (to - from);
    return res;
}

array_view_count :: inline (arr: []$T, from: s64, count: s64) -> []T {
    #if CHECKS {
        assert(from + count <= arr.count, "Index out of bounds: from+count = %, arr.count = %", from+count, arr.count);
        assert(count > 0, "Negative/zero range: count = %", count);
    }

    res := arr;
    res.data += from;
    res.count = count;
    return res;
}





f64 :: inline (num: $T) -> float64 {
    return ncast(float64, num);
}




// Super :: struct (T1: Type, T2: Type) {
//     T :: #run super_type(T1, T2);
// }


super_type :: inline ($a: Type, $b: Type) -> Type {
    #if a == Complex64 || b == Complex64 {
        return Complex64;
    } else
        return float64;
}

is_minor_type :: inline ($base: Type, $b: Type) -> bool {
    #if base == Complex64 {
        return true;
    } else {
        #if b == Complex64 {
            return false;
        } else {
            return true;
        }
    }
}




// Copied over from Hash_Table.jai
next_power_of_two :: inline (x : int) -> int {
    assert(x != 0);
    p := 1;
    while x > p p += p;
    return p;
}





COMBI_V_V : [4][2]string : .[
    .["Vector(float64)", "Vector(float64)"],
    .["Vector(float64)", "Vector(Complex64)"],
    .["Vector(Complex64)", "Vector(float64)"],
    .["Vector(Complex64)", "Vector(Complex64)"],
];
COMBI_V_X : [4][2]string : .[
    .["Vector(float64)", "float64"],
    .["Vector(float64)", "Complex64"],
    .["Vector(Complex64)", "float64"],
    .["Vector(Complex64)", "Complex64"] 
];
COMBI_V : [4][1]string : .[
    .["Vector(float64)"],
    .["Vector(Complex64)"],
    .["Vector(Complex64)"],
    .["Vector(Complex64)"],
];
COMBI_X : [4][1]string : .[
    .["float64"],
    .["Complex64"],
    .["Complex64"],
    .["Complex64"],
];



// meant to be used in the form of
//     #insert #run specialize(...);
// generates functions
//     fun_name :: inline (x0: itypes[i][0], x1: itypes[i][1],...) -> otypes[i][0], otypes[i][1],... {...}
// by accessing the generic function
//     _fun_name :: (x0: $A, x1: $B, ..., xn: *$NA, xn+1: *$NB, ...) {...}
// writing specialized versions with return types;
specialize :: (backend_name: string, fun_name: string, itypes: [][]string, otypes: [][]string) -> string {

    assert(itypes.count == otypes.count);

    builder: String_Builder;
    defer free_buffers(*builder);

    for f: 0..itypes.count-1 {
        n := itypes[f].count;

        print_to_builder(*builder, "% :: inline (", fun_name);
        for i: 0..itypes[f].count-1 {
            print_to_builder(*builder, "x%: %", i, itypes[f][i]);
            if i != itypes[f].count-1 {
                append(*builder, ", ");
            }
        }
        append(*builder, ") -> ");
        for i: 0..otypes[f].count-1 {
            print_to_builder(*builder, "%", otypes[f][i]);
            if i != otypes[f].count-1 {
                append(*builder, ", ");
            }
        }
        
        append(*builder, " {\n");


        for i: 0..otypes[f].count-1 {
            print_to_builder(*builder, 
                "    x% : % = ---;\n", i+n, otypes[f][i]
            );
        }


        print_to_builder(*builder, "    %(", backend_name);
        for i: 0..itypes[f].count-1 {
            print_to_builder(*builder, "x%", i);
            if i != itypes[f].count-1 {
                append(*builder, ", ");
            } else if otypes[f].count != 0 {
                append(*builder, ", ");
            }
        }
        for i: 0..otypes[f].count-1 {
            print_to_builder(*builder, "*x%", i+n);
            if i != otypes[f].count-1 {
                append(*builder, ", ");
            }
        }

        append(*builder, ");\n");
        
        append(*builder, "    return ");
        for i: 0..otypes[f].count-1 {
            print_to_builder(*builder, "x%", i+n);
            if i != otypes[f].count-1 {
                append(*builder, ", ");
            }
        }
        append(*builder, ";\n");

        append(*builder, "}\n");
    }

    s := builder_to_string(*builder);
    print("Functions created: \n%\n", s);
    return s;
}










test_utils :: () {
    print("\n\nUtils -------------------------------\n\n");
    assert(sign(2.0) == 1.0);
    assert(sign(-2.0) == -1.0);
    assert(sign(-2) == -1);

    assert(binomial(42, 42) == 1);
    assert(binomial(42, 0) == 1);
    assert(binomial(5, 3) == 10);
    assert(binomial(11, 6) == 462);
    // print("%\n", binomial(67, 30));

    {
        a := float64.[0,1,2,3,4,5,6];
        v := array_view_range(a, 2, 4);
        print("%\n", a);
        print("%\n", v);
    }

    {
        assert(Complex64 == super_type(Complex64, u8));
        assert(float64 == super_type(float32, float64));
    }

    // {
    //     a : float64 = ---;
    //     M : CMatrix64 = ---;
    //     calculation(#code {
    //         x := 42.0; // 8 bytes
    //         a = x;

    //         copy(*M, cmat64_id(4,4)); // 2 * 8 bytes [Complex] * 16 = 256 bytes
    //     }); // Pool allocated 264 bytes!
    //     assert(a == 42.0);
    //     print("M = %\n", M);
    //     print("%\n", cmat64_id(4,4));
    // }
}