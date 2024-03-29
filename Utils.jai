


clamp :: inline (mi: $T, val: T, ma: T) -> T {
    return min(ma, max(mi, val));
}


min :: (values: ..$T) -> T {
    res := values[0];
    for i: 1..values.count-1 {
        if values[i] < res {
            res = values[i];
        }
    }
    return res;
}

max :: (values: ..$T) -> T {
    res := values[0];
    for i: 1..values.count-1 {
        if values[i] > res {
            res = values[i];
        }
    }
    return res;
}


array :: ($N: int, values: ..$T) -> [N]T {
    #if CHECKS {
        assert(N == values.count);
    }
    res : [N]T = ---;
    for values {
        res[it_index] = it;
    }
    return res;
}

array_equals :: (x: []$T, y: []T) -> bool {
    if x.count != y.count then return false;
    for x {
        if it != y[it_index] then return false;
    }
    return true;
}



pool :: (code: Code) #expand {
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

// Copied over from Hash_Table.jai
next_power_of_two :: inline (x : int) -> int {
    assert(x != 0);
    p := 1;
    while x > p p += p;
    return p;
}

str :: inline (val: $T) -> string #modify {
    if !is_intrinsic_number(T) { T = null; return; }
} {
    return sprint("%", val);
}





f64 :: inline (x: $T) -> float64 #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    return ncast(float64, x);
}
f32 :: inline (x: $T) -> float64 #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    return ncast(float32, x);
}



swap_code :: (a: Code, b: Code) #expand {
    t := (#insert a);
    (#insert a) = (#insert b);
    (#insert b) = t;
}


_if :: ($$condition: bool, body_true: Code, body_false: Code = #code {}) #expand {
    #if is_constant(condition) {
        #if condition {
            #insert body_true;
        } else {
            #insert body_false;
        }
    } else {
        if condition {
            #insert body_true;
        } else {
            #insert body_false;
        }
    }
}

for_unroll :: (from: int, to: int, body: Code, reverse: bool = false) #expand {
  #insert -> string {
    builder: String_Builder;
    init_string_builder(*builder);
    defer free_buffers(*builder);

    append(*builder, "{\n");
    append(*builder, "    `it: int;\n");
    if reverse {
        for <i:a..b {
          print_to_builder(*builder, "    it = %;\n", i);
          append(*builder, "    #insert body;\n");
        }
    } else {
        for a..b {
          print_to_builder(*builder, "    it = %;\n", it);
          append(*builder, "    #insert body;\n");
        }
    }
    append(*builder, "}\n");
    return builder_to_string(*builder);
  }
}




read_file_and_split :: (file_url: string) -> [][]string {
    data, ok := File.read_entire_file(file_url);
    if !ok {
        assert(false, "Could not read entire file: %\n", file_url);
    }

    lines := String.split(data, "\n");

    res := NewArray(lines.count, []string);

    for lines {
        res[it_index] = String.split(it, " ");
    }

    return res;
}




test_utils :: () {
    println_push("Utils", color = .FG_WHITE);

    {
        println("test 1 : assert sign", color = .FG_GREEN);
        assert(sign(2.0) == 1.0);
        assert(sign(-2.0) == -1.0);
        assert(sign(-2) == -1);
    }

    {
        println_push("test 2 : assert binomial", color = .FG_GREEN);
        println("assert 1");
        assert(binomial(42, 42) == 1);
        println("assert 2");
        assert(binomial(42, 0) == 1);
        println("assert 3");
        assert(binomial(5, 3) == 10);
        println("assert 4");
        assert(binomial(11, 6) == 462);
        // println("print result");
        // print("%\n", binomial(67, 30));
    }
}