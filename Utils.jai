





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







is_a_matrix :: inline (T: Type) -> bool {
    return T == Matrix(Complex64)
        || T == Matrix(float64)
        || T == Matrix(float32)
        || T == Matrix(float)
        || T == Matrix(int)
        || T == Matrix(s64)
        || T == Matrix(s32)
        || T == Matrix(s16)
        || T == Matrix(u64)
        || T == Matrix(u32)
        || T == Matrix(u16)
        || T == Matrix(u8)
        ;
}
is_a_vector :: inline (T: Type) -> bool {
    return T == Vector(Complex64)
        || T == Vector(float64)
        || T == Vector(float32)
        || T == Vector(float)
        || T == Vector(int)
        || T == Vector(s64)
        || T == Vector(s32)
        || T == Vector(s16)
        || T == Vector(u64)
        || T == Vector(u32)
        || T == Vector(u16)
        || T == Vector(u8)
        ;
}
is_a_scalar :: inline (T: Type) -> bool {
    return T == Complex64
        || T == float64
        || T == float32
        || T == float
        || T == int
        || T == s64
        || T == s32
        || T == s16
        || T == u64
        || T == u32
        || T == u16
        || T == u8
        ;
}







ncast :: inline ($T: Type, num: $N) -> T {
    #if CHECKS {
        assert(is_a_scalar(T) && is_a_scalar(N), "cannot cast between non-scalars! T = %, N = %\n", T, N);
    }

    // #if #run !is_a_scalar(T) || !is_a_scalar(N) {
    //     print("cannot cast between non-scalars! T = %, N = %\n", T, N);
    //     assert(false);
    // }

    #if T == N {
        return num;
    } else #if T == Complex64 {
        return c64(cast(float64)num);
    } else #if #run is_a_scalar(N) && N != Complex64 {
        // print("num = %, N = %\n", num, N);
        return cast,no_check(T) num; // TODO other types?
    } else {
        assert(N != Complex64, "cannot cast complex to real number!");
        return zero(T);
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
    #if #run is_a_vector(T) {
        return norm_2(a - b) < DELTA;
    } else #if #run is_a_matrix(T) {
        return norm_frobenius(a - b) < DELTA;
    } else
        return abs(a - b) < DELTA;
}




#scope_module


array_view_range :: inline (arr: []$T, from: s64, to: s64) -> []T {
    #if CHECKS {
        assert(to < arr.count, "Index out of bounds: to = %, arr.count = %" , to, arr.count); 
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



COMBI_V_V_V :: string.[
    "Vector(float64)   & Vector(float64)   => Vector(float64)",
    "Vector(float64)   & Vector(Complex64) => Vector(Complex64)",
    "Vector(Complex64) & Vector(float64)   => Vector(Complex64)",
    "Vector(Complex64) & Vector(Complex64) => Vector(Complex64)",
];
COMBI_V_V_M :: string.[
    "Vector(float64)   & Vector(float64)   => Matrix(float64)",
    "Vector(float64)   & Vector(Complex64) => Matrix(Complex64)",
    "Vector(Complex64) & Vector(float64)   => Matrix(Complex64)",
    "Vector(Complex64) & Vector(Complex64) => Matrix(Complex64)",
];
COMBI_V_X_V :: string.[
    "Vector(float64)   & float64   => Vector(float64)",
    "Vector(float64)   & Complex64 => Vector(Complex64)",
    "Vector(Complex64) & float64   => Vector(Complex64)",
    "Vector(Complex64) & Complex64 => Vector(Complex64)",
];
COMBI_V_V_X :: string.[
    "Vector(float64)   & Vector(float64)   => float64",
    "Vector(float64)   & Vector(Complex64) => Complex64",
    "Vector(Complex64) & Vector(float64)   => Complex64",
    "Vector(Complex64) & Vector(Complex64) => Complex64",
];


COMBI_M_M_M :: string.[
    "Matrix(float64)   & Matrix(float64)   => Matrix(float64)",
    "Matrix(float64)   & Matrix(Complex64) => Matrix(Complex64)",
    "Matrix(Complex64) & Matrix(float64)   => Matrix(Complex64)",
    "Matrix(Complex64) & Matrix(Complex64) => Matrix(Complex64)",
];
COMBI_M_V_V :: string.[
    "Matrix(float64)   & Vector(float64)   => Vector(float64)",
    "Matrix(float64)   & Vector(Complex64) => Vector(Complex64)",
    "Matrix(Complex64) & Vector(float64)   => Vector(Complex64)",
    "Matrix(Complex64) & Vector(Complex64) => Vector(Complex64)",
];
COMBI_V_M_V :: string.[
    "Vector(float64)   & Matrix(float64)   => Vector(float64)",
    "Vector(float64)   & Matrix(Complex64) => Vector(Complex64)",
    "Vector(Complex64) & Matrix(float64)   => Vector(Complex64)",
    "Vector(Complex64) & Matrix(Complex64) => Vector(Complex64)",
];
COMBI_M_X_M :: string.[
    "Matrix(float64)   & float64   => Matrix(float64)",
    "Matrix(float64)   & Complex64 => Matrix(Complex64)",
    "Matrix(Complex64) & float64   => Matrix(Complex64)",
    "Matrix(Complex64) & Complex64 => Matrix(Complex64)",
];



COMBI_M_M_b :: string.[
    "Matrix(float64)   & Matrix(float64)   => bool",
    "Matrix(float64)   & Matrix(Complex64) => bool",
    "Matrix(Complex64) & Matrix(float64)   => bool",
    "Matrix(Complex64) & Matrix(Complex64) => bool",
];
COMBI_V_V_b :: string.[
    "Vector(float64)   & Vector(float64)   => bool",
    "Vector(float64)   & Vector(Complex64) => bool",
    "Vector(Complex64) & Vector(float64)   => bool",
    "Vector(Complex64) & Vector(Complex64) => bool",
];

COMBI_M_M :: string.[
    "Matrix(float64)   => Matrix(float64)",
    "Matrix(Complex64) => Matrix(Complex64)",
];
COMBI_V_V :: string.[
    "Vector(float64)   => Vector(float64)",
    "Vector(Complex64) => Vector(Complex64)",
];


// meant to be used in the form of
//     #insert #run specialize(...);
// generates functions
//     fun_name :: inline (x0: itypes[i][0], x1: itypes[i][1],...) -> otypes[i][0], otypes[i][1],... {...}
// by accessing the generic function
//     _fun_name :: (x0: $A, x1: $B, ..., xn: *$NA, xn+1: *$NB, ...) {...}
// writing specialized versions with return types;
specialize :: (backend_name: string, fun_name: string, iotypes: []string, o_opts: []string = .[]) -> string {


    builder: String_Builder;
    defer free_buffers(*builder);

    for f: 0..iotypes.count-1 {

        itypes : []string = ---;
        otypes : []string = ---;
        {
            tmp := String.split(iotypes[f], "=>");
            // print("splitted at => : %\n", tmp);

            itmp := tmp[0];
            otmp := tmp[1];

            itypes = String.split(itmp, "&");
            otypes = String.split(otmp, "&");

            // print("itypes : %\n", itypes);
            // print("otypes : %\n", otypes);
        }

        n := itypes.count;

        print_to_builder(*builder, "% :: inline (", fun_name);
        for i: 0..itypes.count-1 {
            print_to_builder(*builder, "x%: %", i, itypes[i]);
            if i != itypes.count-1 {
                append(*builder, ", ");
            }
        }
        append(*builder, ") ");

        if otypes.count != 0 {
            append(*builder, "-> ");
            for i: 0..otypes.count-1 {
                print_to_builder(*builder, "%", otypes[i]);
                
                // options for output types such as #symmetric, #must, etc.
                if i < o_opts.count {
                    print_to_builder(*builder, " %", o_opts[i]);
                }
                
                if i != otypes.count-1 {
                    append(*builder, ", ");
                }
            }
        }
        
        append(*builder, " {\n");


        for i: 0..otypes.count-1 {
            print_to_builder(*builder, 
                "    x% : % = ---;\n", i+n, otypes[i]
            );
        }


        print_to_builder(*builder, "    %(", backend_name);
        for i: 0..itypes.count-1 {
            print_to_builder(*builder, "x%", i);
            if i != itypes.count-1 {
                append(*builder, ", ");
            } else if otypes.count != 0 {
                append(*builder, ", ");
            }
        }
        for i: 0..otypes.count-1 {
            print_to_builder(*builder, "*x%", i+n);
            if i != otypes.count-1 {
                append(*builder, ", ");
            }
        }

        append(*builder, ");\n");
        
        append(*builder, "    return ");
        for i: 0..otypes.count-1 {
            print_to_builder(*builder, "x%", i+n);
            if i != otypes.count-1 {
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