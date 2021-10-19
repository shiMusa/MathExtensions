/*
 ######  ######## ########
##    ##    ##    ##     ##
##          ##    ##     ##
 ######     ##    ########
      ##    ##    ##   ##
##    ##    ##    ##    ##
 ######     ##    ##     ##
*/


pstr :: inline (v: $V/VectorType) -> string {
    return str(v, print_options.offset, print_options.shift, "│", print_options.indent_charactoer_every);
}

str :: (v: $V/VectorType, offset: int = 0, indent: int = 4, indent_character: string = " ", indent_character_every: int = 4) -> string {
    builder: String_Builder;
    defer free_buffers(*builder);

    init_string_builder(*builder);
    append(*builder, "\n");

    shift :: (ind: int = 0) #expand {
        for 0..offset-1 {
            append(*builder, " ");
        }
        for 0..ind-1 {
            if it % indent_character_every == 0 {
                append(*builder, indent_character);
            } else {
                append(*builder, " ");
            }
        }
    }

    shift(indent+4);
    append(*builder, "(");
    for v {
        append(*builder, str(it));
        if it_index != dim(v)-1 {
            append(*builder, ", ");
        }
    }
    append(*builder, ")\n");

    shift(indent);
    append(*builder, "└╴");

    return builder_to_string(*builder);
}



/*
#### ##    ## #### ######## ####    ###    ##       #### ######## ########
 ##  ###   ##  ##     ##     ##    ## ##   ##        ##       ##  ##
 ##  ####  ##  ##     ##     ##   ##   ##  ##        ##      ##   ##
 ##  ## ## ##  ##     ##     ##  ##     ## ##        ##     ##    ######
 ##  ##  ####  ##     ##     ##  ######### ##        ##    ##     ##
 ##  ##   ###  ##     ##     ##  ##     ## ##        ##   ##      ##
#### ##    ## ####    ##    #### ##     ## ######## #### ######## ########
*/

// * Assuming `to` is initialized to zeroes.
// * User has to know, what type the result has to be!
cast_to :: (to: *$W/VectorType, from: $V/VectorType) {
    for raw(*from) {
        set(to, it_index, ncast(W.data_type, it));
    }
}

// * Assuming `to` is initialized to zeroes.
// * User has to know, what type the result has to be!
copy :: (to: *$V/VectorType, from: V) {
    for raw(*from) {
        set(to, it_index, it);
    }
}


vec_init :: inline (v: *$V/VectorType, values: ..$T) {
    #if CHECKS {
        #assert(is_minor_number_type(V.data_type, T));
        assert(values.count <= dim(v), "too many values (%) for vector of % entries\n", values.count, dim(v));
    }
    for values {
        set(v, it_index, ncast(V.data_type, it));
    }
}

vec_init_all :: inline (v: *$V/VectorType, value: $T) {
    #if CHECKS {
        #assert(is_minor_number_type(V.data_type, T));
    }
    for raw(v) {
        set_raw(v, it_raw_index, ncast(V.data_type, value));
    }
}



///////////////////////////////////////////////////////////////////////////////////////
/*
      ######## ########  ######  ########  ######
           ##    ##       ##    ##    ##    ##    ##
           ##    ##       ##          ##    ##
           ##    ######    ######     ##     ######
           ##    ##             ##    ##          ##
           ##    ##       ##    ##    ##    ##    ##
           ##    ########  ######     ##     ######
*/
////////////////////////////////////////////////////////////////////////////////////////

#scope_module
test_vector :: () {

    println_push("Vector", color = .FG_WHITE);

    {
        println_push("DenseVector", color = .FG_GREEN);    
        {
            v := make(DenseVector(float64, 5));
            println("v = %", pstr(v));
            vec_init(*v, 1.0, 2.0, 3.0, 4.0, 5.0);
            println("v = %", pstr(v));
            // print("%, type_of(v) = %, Vector(T) = Vector(%)\n", v, T.vector_type, T.data_type);
            a := v + v;
            println("a = %", pstr(v));
            a = a + v;
            println("a = %", pstr(v));
        }
    }

    {
        println_push("DenseHeapVector", color = .FG_GREEN);    
        {
            v := make(DenseHeapVector(float64, 5));
            defer free(v);
            println("v = %", pstr(v));
            vec_init(*v, 1.0, 2.0, 3.0, 4.0, 5.0);
            println("v = %", pstr(v));
            // print("%, type_of(v) = %, Vector(T) = Vector(%)\n", v, T.vector_type, T.data_type);
            a := v + v;
            println("a = %", pstr(v));
            a = a + v;
            println("a = %", pstr(v));
        }
    }
}

