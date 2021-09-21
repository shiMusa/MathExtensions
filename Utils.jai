





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







// is_a_matrix :: inline (T: Type) -> bool {
//     return T == Matrix(Complex64)
//         || T == Matrix(float64)
//         || T == Matrix(float32)
//         || T == Matrix(float)
//         || T == Matrix(int)
//         || T == Matrix(s64)
//         || T == Matrix(s32)
//         || T == Matrix(s16)
//         || T == Matrix(u64)
//         || T == Matrix(u32)
//         || T == Matrix(u16)
//         || T == Matrix(u8)
//         ;
// }
// is_a_vector :: inline (T: Type) -> bool {
//     return T == Vector(Complex64)
//         || T == Vector(float64)
//         || T == Vector(float32)
//         || T == Vector(float)
//         || T == Vector(int)
//         || T == Vector(s64)
//         || T == Vector(s32)
//         || T == Vector(s16)
//         || T == Vector(u64)
//         || T == Vector(u32)
//         || T == Vector(u16)
//         || T == Vector(u8)
//         ;
// }


// ######## ##    ## ########  ########  ######
//    ##     ##  ##  ##     ## ##       ##    ##
//    ##      ####   ##     ## ##       ##
//    ##       ##    ########  ######    ######
//    ##       ##    ##        ##             ##
//    ##       ##    ##        ##       ##    ##
//    ##       ##    ##        ########  ######


NUMBER_TYPES :: Type.[
    float64, float32, s64, s32, s16, s8, u64, u32, u16, u8, bool
];



is_intrinsic_number :: ($T: Type) -> bool {
    for NUMBER_TYPES {
        if T == it then return true;
    }
    return false;
}

is_intrinsic_number :: (n: *Type_Info) -> bool {
    return n.type == .FLOAT
        || n.type == .INTEGER
        || n.type == .BOOL;
}


is_complex_number :: ($T: Type) -> bool {

    #insert -> string {
        builder: String_Builder;

        for NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == Complex(%) then return true;\n",
                it
            );
        }

        return builder_to_string(*builder);
    }

    return false;
}

is_complex_number :: (n: *Type_Info) -> bool {
    if n.type != .STRUCT then return false;
    s := cast(*Type_Info_Struct)n;

    if s.name == "Complex" then return true;

    return false;
}



is_number :: inline ($T: Type) -> bool {
    return is_intrinsic_number(T) || is_complex_number(T);
}

is_number :: inline (n: *Type_Info) -> bool {
    return is_intrinsic_number(n) || is_complex_number(n);
}




// ? assuming we already know it's a Complex(T) number
get_complex_inner_type :: ($T: Type) -> Type {
    #insert -> string {
        builder: String_Builder;

        for NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == Complex(%) then return %;\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return Any;
}

get_complex_inner_type :: (T: *Type_Info) -> *Type_Info {
    s := cast(*Type_Info_Struct)T;
    return s.members[0].type;
}


get_complex_type_info_from_intrinsic_type_info :: (T: *Type_Info) -> *Type_Info {
    #insert -> string {
        builder: String_Builder;

        for NUMBER_TYPES {
            print_to_builder(*builder,
                "if cast(*Type_Info)type_info(%) == T then return type_info(Complex(%));\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return null;
}



// ##    ##  ######     ###     ######  ########
// ###   ## ##    ##   ## ##   ##    ##    ##
// ####  ## ##        ##   ##  ##          ##
// ## ## ## ##       ##     ##  ######     ##
// ##  #### ##       #########       ##    ##
// ##   ### ##    ## ##     ## ##    ##    ##
// ##    ##  ######  ##     ##  ######     ##


// * I don't need to use #modify since this is only used in terms of numbers anyways
ncast :: inline ($T: Type, num: $N) -> T {
    #if T == N {
        return num;
    } else #if #run is_complex_number(T) {
        X :: #run get_complex_inner_type(T);
        return cpx(cast(X)num);
    } else #if #run !is_complex_number(N) {
        return cast,no_check(T) num;
    } else {
        assert(type_info(N).name != "Complex", "cannot cast complex to real number!");
        return zero(T);
    }
}

zero :: inline ($T: Type) -> T {
    return ncast(T, 0);
}

one :: inline ($T: Type) -> T {
    return ncast(T, 1);
}

nan :: inline ($T: Type) -> T {
    #if T == float64 {
        return Math.FLOAT64_NAN;
    } else #if T == float32 {
        return Math.FLOAT32_NAN;
    } else #if #run is_complex_number(T) {
        X :: #run get_complex_inner_type(T);
        n :: #run nan(X);
        return #run Complex(X).{n,n};
    } else return null;
}





#scope_module


array_copy :: (from: [$N]$T, to: *[N]T) {
    for 0..N-1 {
        (<<to)[it] = from[it];
    }
}


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


// ##     ## ########   ######     ###     ######  ########
// ##     ## ##     ## ##    ##   ## ##   ##    ##    ##
// ##     ## ##     ## ##        ##   ##  ##          ##
// ##     ## ########  ##       ##     ##  ######     ##
// ##     ## ##        ##       #########       ##    ##
// ##     ## ##        ##    ## ##     ## ##    ##    ##
//  #######  ##         ######  ##     ##  ######     ##


upcast :: (A: *Type_Info, B: *Type_Info) -> *Type_Info {
    if !is_number(A) || !is_number(B) then return null;

    println_push("upcast (%, %)", A.type, B.type, color = C_FG_WHITE);

    if is_complex_number(A) && is_complex_number(B) {
        a := get_complex_inner_type(A);
        b := get_complex_inner_type(B);
        c := upcast(a,b);
        println("c = %", c.type);
        return get_complex_type_info_from_intrinsic_type_info(c);
    }

    if is_complex_number(A) {
        a := get_complex_inner_type(A);
        c := upcast(a,B);
        println("c = %", c.type);
        return get_complex_type_info_from_intrinsic_type_info(c);
    }

    if is_complex_number(B) {
        b := get_complex_inner_type(B);
        c := upcast(A,b);
        println("c = %", c.type);
        return get_complex_type_info_from_intrinsic_type_info(c);
    }

    #insert -> string {
        builder : String_Builder;

        for NUMBER_TYPES {
            print_to_builder(*builder,
                "if A == type_info(%) || B == type_info(%) then return type_info(%);",
                it, it, it
            );
        }

        return builder_to_string(*builder);
    }

    return null;
}


upcast :: inline ($A: Type, $B: Type) -> Type {
    is_A_complex :: #run is_complex_number(A);
    is_B_complex :: #run is_complex_number(B);

    #if is_A_complex && is_B_complex {
        a :: #run get_complex_inner_type(A);
        b :: #run get_complex_inner_type(B);
        T :: #run upcast(a, b);
        return Complex(T);
    } else #if is_A_complex {
        a :: #run get_complex_inner_type(A);
        T :: #run upcast(a, B);
        return Complex(T);
    } else #if is_B_complex {
        b :: #run get_complex_inner_type(B);
        T :: #run upcast(A,b);
        return Complex(T);
    } else {
        // ? the types in NUMBER_TYPES are sorted in order from mathematical superset to subset, e.g. ℂ > ℝ > ℤ > ℕ 
        for NUMBER_TYPES {
            if A == it || B == it then return it;
        }
        return Any;
    }
}


UpCast :: struct(A: Type, B: Type) {
    T :: #run upcast(A, B);
}



is_minor_number :: inline ($major: Type, $minor: Type) -> bool {
    i_major := -1;
    i_minor := -1;
    for NUMBER_TYPES {
        if it == major then i_major = it_index;
        if it == minor then i_minor = it_index;
    }
    return i_major <= i_minor; // * higher-order types come first!
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




// ########  ########  #### ##    ## ########
// ##     ## ##     ##  ##  ###   ##    ##
// ##     ## ##     ##  ##  ####  ##    ##
// ########  ########   ##  ## ## ##    ##
// ##        ##   ##    ##  ##  ####    ##
// ##        ##    ##   ##  ##   ###    ##
// ##        ##     ## #### ##    ##    ##




PrintOptions :: struct {
    offset := 0;
    shift := 0;
    indent_every := 4;
    indent_character := "│";
};
print_options : PrintOptions;

print_indent :: (levels: int) {
    print_options.shift += 4 * levels;
}


C_FG_RESET   :: "\e[0m";
C_FG_BLACK   :: "\e[30m";
C_FG_RED     :: "\e[31m";
C_FG_GREEN   :: "\e[32m";
C_FG_YELLOW  :: "\e[33m";
C_FG_BLUE    :: "\e[34m";
C_FG_MAGENTA :: "\e[35m";
C_FG_CYAN    :: "\e[36m";
C_FG_WHITE   :: "\e[37m";


fg_color :: inline (txt: string, color: string) -> string {
    return sprint("%1%2%3", color, txt, C_FG_RESET);
}


println_push :: (text: string = "", args: ..Any, color: string = C_FG_RESET) #expand {
    c_text := ifx color != C_FG_RESET then sprint("%1%2%3", color, text, C_FG_RESET) else text;
    println(c_text, ..args);
    print_indent(1);
    `defer println("└╴");
    `defer print_indent(-1);
}

println :: (text: string = "", args: ..Any, color: string = C_FG_RESET) {

    c_text := ifx color != C_FG_RESET then sprint("%1%2%3", color, text, C_FG_RESET) else text;

    for 0..print_options.offset-1 {
        print(" ");
    }

    for 0..print_options.shift-1 {
        if it % print_options.indent_every == 0 {
            print(print_options.indent_character);
        } else
            print(" ");
    }
    print(c_text, ..args);
    print("\n");
}

println :: (arg: Any, color: string = C_FG_RESET) {
    text := sprint("%", arg);
    c_text := ifx color != C_FG_RESET then sprint("%1%2%3", color, text, C_FG_RESET) else text;
    
    for 0..print_options.offset-1 {
        print(" ");
    }

    for 0..print_options.shift-1 {
        if it % print_options.indent_every == 0 {
            print(print_options.indent_character);
        } else
            print(" ");
    }
    print("%\n", c_text);
}




test_utils :: () {
    println_push("Utils", color = C_FG_WHITE);

    {
        println("test 1 : assert sign", color = C_FG_GREEN);
        assert(sign(2.0) == 1.0);
        assert(sign(-2.0) == -1.0);
        assert(sign(-2) == -1);
    }

    {
        println_push("test 2 : assert binomial", color = C_FG_GREEN);
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

    {
        println_push("test 3 : array_view_range", color = C_FG_GREEN);
        a := float64.[0,1,2,3,4,5,6];
        v := array_view_range(a, 2, 4);
        println("%", a);
        println("%", v);
    }

    // {
    //     assert(Complex64 == super_type(Complex64, u8));
    //     assert(float64 == super_type(float32, float64));
    // }

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