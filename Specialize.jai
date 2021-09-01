







V_V_V :: string.[
    "Vector(float64, $N)   & Vector(float64, N)   => Vector(float64, N)",
    "Vector(Complex64, $N) & Vector(float64, N)   => Vector(Complex64, N)",
    "Vector(float64, $N)   & Vector(Complex64, N) => Vector(Complex64, N)",
    "Vector(Complex64, $N) & Vector(Complex64, N) => Vector(Complex64, N)",
];

// M_M_M :: string.[
//     "Matrix(float64, $N, $M)   & Matrix(float64, N)   => Matrix(float64, N)",
//     "Matrix(Complex64, $N, $M) & Matrix(float64, N)   => Matrix(Complex64, N)",
//     "Matrix(float64, $N, $M)   & Matrix(Complex64, N) => Matrix(Complex64, N)",
//     "Matrix(Complex64, $N, $M) & Matrix(Complex64, N) => Matrix(Complex64, N)",
// ];






// meant to be used in the form of
//     #insert #run specialize(...);
// generates functions
//     fun_name :: inline (x0: itypes[i][0], x1: itypes[i][1],...) -> otypes[i][0], otypes[i][1],... {...}
// by accessing the generic function
//     _fun_name :: (x0: $A, x1: $B, ..., xn: *$NA, xn+1: *$NB, ...) {...}
// writing specialized versions with return types;
// E.g.
// _foo :: (a: $A, b: int, c: *$C) {...}
// #insert #run specialize("_foo", "foo", .[
//     "float & int => float",
//     "int   & int => int"
// ]);
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

