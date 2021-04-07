#import "Basic";




// #insert "a :: 42;";

// fct :: #string FCT

// foo :: (a: float64) -> float64 {
//     return 2*a;
// }

// FCT
// #insert fct;


get_variants :: (types: ..Type) -> string {
    builder: String_Builder;
    defer free_buffers(*builder);

    for types {
        print("encountered type %\n", it);

        msg := sprint("    print(\"This is of type %\\n\", );\n", it);
        
        s := sprint("foo :: (a: %) -> % {\n%    return 2*a; \n}\n", it, it, msg);
        print("generated function: \n%\n", s);
        print_to_builder(*builder, s);
    }

    return builder_to_string(*builder);
}

#insert #run get_variants(u8, float64);



// bar :: (a: float64) -> #insert "float64" {
//     return a;
// }



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
                "    x% : %;\n", i+n, otypes[f][i]
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



#insert #run specialize("_bar", "bar",.[.["int"], .["float64"]], .[.["int"], .["float64"]]);

_bar :: (a: $A, b: *$B) {
    <<b = cast(B) ( cast(A)(2) * a );
}



main :: () {
    print("%\n", foo(2));
    print("%\n", foo(2.0));

    x : float64;
    _bar(4.0, *x);
    print("x = %\n", x);

    print("%\n", bar(4.0));
}