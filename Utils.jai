



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


// ######## ##    ## ########  ########  ######
//    ##     ##  ##  ##     ## ##       ##    ##
//    ##      ####   ##     ## ##       ##
//    ##       ##    ########  ######    ######
//    ##       ##    ##        ##             ##
//    ##       ##    ##        ##       ##    ##
//    ##       ##    ##        ########  ######


// ? this is used for matrices and vectors!
make :: inline ($T: Type, $initialize_values: bool = true) -> T {
    res : T = ---;
    init(*res, initialize_values);
    return res;
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





#scope_module

INTRINSIC_NUMBER_TYPES :: Type.[
    float64, float32, s64, s32, s16, s8, u64, u32, u16, u8, bool
];



is_float_number :: inline ($T: Type) -> bool {
    #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        return is_float_number(X);
    } else #if #run is_quaternion_number(T) {
        X :: #run get_number_inner_type(T);
        return is_float_number(X);
    } else #if T == float64 || T == float32 {
        return true;
    } else return false;
}



is_intrinsic_number :: ($T: Type) -> bool {
    for INTRINSIC_NUMBER_TYPES {
        if T == it then return true;
    }
    return false;
}

is_intrinsic_number :: (T: *Type_Info) -> bool {
    return T.type == .FLOAT
        || T.type == .INTEGER
        || T.type == .BOOL;
}



is_non_real_number :: inline ($T: Type) -> bool {
    return is_complex_number(T) || is_quaternion_number(T);
}
is_non_real_number :: inline (T: *Type_Info) -> bool {
    return is_complex_number(T) || is_quaternion_number(T);
}


is_complex_number :: ($T: Type) -> bool {

    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "#if T == Complex(%) then return true;\n",
                it
            );
        }

        return builder_to_string(*builder);
    }

    return false;
}

is_complex_number :: (T: *Type_Info) -> bool {
    if T.type != .STRUCT then return false;
    S := cast(*Type_Info_Struct)T;

    if S.name == "Complex" then return true;

    return false;
}


is_quaternion_number :: ($T: Type) -> bool {

    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "#if T == Quaternion(%) then return true;\n",
                it
            );
        }

        return builder_to_string(*builder);
    }

    return false;
}

is_quaternion_number :: (N: *Type_Info) -> bool {
    if N.type != .STRUCT then return false;
    S := cast(*Type_Info_Struct)N;

    if S.name == "Quaternion" then return true;

    return false;
}


is_number :: inline ($T: Type) -> bool {
    return is_intrinsic_number(T) || is_complex_number(T) || is_quaternion_number(T);
}

is_number :: inline (T: *Type_Info) -> bool {
    return is_intrinsic_number(T) || is_complex_number(T) || is_quaternion_number(T);
}






get_number_inner_type :: ($T: Type) -> Type {
    #if #run is_complex_number(T) {
        return #run get_complex_inner_type(T);
    } else #if #run is_quaternion_number(T) {
        return #run get_quaternion_inner_type(T);
    } else {
        return T;
    }
}
get_number_inner_type :: (T: *Type_Info) -> *Type_Info {
    if is_intrinsic_number(T) then return T;

    S := cast(*Type_Info_Struct)T;
    return S.members[0].type;
}



// ? assuming we already know it's a Complex(T) number
get_complex_inner_type :: ($T: Type) -> Type {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "#if T == Complex(%) then return %;\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return Any;
}


get_complex_type_info_from_intrinsic_type_info :: (T: *Type_Info) -> *Type_Info {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if cast(*Type_Info)type_info(%) == T then return type_info(Complex(%));\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return null;
}






// ? assuming we already know it's a Complex(T) number
get_quaternion_inner_type :: ($T: Type) -> Type {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "#if T == Quaternion(%) then return %;\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return Any;
}


get_quaternion_type_info_from_intrinsic_type_info :: (T: *Type_Info) -> *Type_Info {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if cast(*Type_Info)type_info(%) == T then return type_info(Quaternion(%));\n",
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


// * I don't need to use #modify since this is only used in terms of numbers anyways!?
#scope_export
ncast :: inline ($T: Type, num: $N) -> T {
    #if T == N {

        return num;

    } else #if #run is_quaternion_number(T) {
        
        X :: #run get_number_inner_type(T);

        #if #run is_quaternion_number(N) {
            return qtn(cast(X)num.r, cast(X)num.i, cast(X)num.j, cast(X)num.k);
        } else #if #run is_complex_number(N) {
            return qtn(cast(X)num.real, cast(X)num.imag, 0, 0);
        } else {
            return qtn(cast(X)num, 0, 0, 0);
        }

    } else #if #run is_complex_number(T) {
        
        #if #run is_quaternion_number(N) {
            assert(false, "cannot ncast % -> %", N, T);
        } else {

            X :: #run get_number_inner_type(T);
            #if #run is_complex_number(N) {
                return cpx(cast(X)num.real, cast(X)num.imag);
            } else {
                return cpx(cast(X)num, zero(X));
            }
        }

    } else { // real numbers

        #if #run is_quaternion_number(N) || is_complex_number(N) {
            assert(false, "cannot ncast % -> %", N, T);
        } else {

            return cast(T)num;

        }
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
        X :: #run get_number_inner_type(T);
        n :: #run nan(X);
        return #run Complex(X).{n,n};
    } else #if #run is_quaternion_number(T) {
        X :: #run get_number_inner_type(T);
        n :: #run nan(X);
        return #run Quaternion(X).{n,n,n,n};
    } else return null;
}






// ##     ## ########   ######     ###     ######  ########
// ##     ## ##     ## ##    ##   ## ##   ##    ##    ##
// ##     ## ##     ## ##        ##   ##  ##          ##
// ##     ## ########  ##       ##     ##  ######     ##
// ##     ## ##        ##       #########       ##    ##
// ##     ## ##        ##    ## ##     ## ##    ##    ##
//  #######  ##         ######  ##     ##  ######     ##


#scope_module
upcast :: (A: *Type_Info, B: *Type_Info) -> *Type_Info {
    if !is_number(A) || !is_number(B) then return null;

    if is_quaternion_number(A) && is_quaternion_number(B) {
        a := get_number_inner_type(A);
        b := get_number_inner_type(B);
        c := upcast(a,b);
        println("c = %", c.type);
        return get_quaternion_type_info_from_intrinsic_type_info(c);
    }

    if is_quaternion_number(A) {
        a := get_number_inner_type(A);
        c := upcast(a,B);
        println("c = %", c.type);
        return get_quaternion_type_info_from_intrinsic_type_info(c);
    }

    if is_quaternion_number(B) {
        b := get_number_inner_type(B);
        c := upcast(A,b);
        println("c = %", c.type);
        return get_quaternion_type_info_from_intrinsic_type_info(c);
    }

    if is_complex_number(A) && is_complex_number(B) {
        a := get_number_inner_type(A);
        b := get_number_inner_type(B);
        c := upcast(a,b);
        println("c = %", c.type);
        return get_complex_type_info_from_intrinsic_type_info(c);
    }

    if is_complex_number(A) {
        a := get_number_inner_type(A);
        c := upcast(a,B);
        println("c = %", c.type);
        return get_complex_type_info_from_intrinsic_type_info(c);
    }

    if is_complex_number(B) {
        b := get_number_inner_type(B);
        c := upcast(A,b);
        println("c = %", c.type);
        return get_complex_type_info_from_intrinsic_type_info(c);
    }

    #insert -> string {
        builder : String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if A == type_info(%) || B == type_info(%) then return type_info(%);\n",
                it, it, it
            );
        }

        return builder_to_string(*builder);
    }

    return null;
}


upcast :: inline ($A: Type, $B: Type) -> Type {
    is_A_quaternion :: #run is_quaternion_number(A);
    is_B_quaternion :: #run is_quaternion_number(B);
    is_A_complex    :: #run is_complex_number(A);
    is_B_complex    :: #run is_complex_number(B);

    #if is_A_quaternion && is_B_quaternion {
        a :: #run get_number_inner_type(A);
        b :: #run get_number_inner_type(B);
        T :: #run upcast(a, b);
        return Quaternion(T);
    } else #if is_A_quaternion {
        a :: #run get_number_inner_type(A);
        T :: #run upcast(a, B);
        return Quaternion(T);
    } else #if is_B_quaternion {
        b :: #run get_number_inner_type(B);
        T :: #run upcast(A,b);
        return Quaternion(T);

    } else #if is_A_complex && is_B_complex {
        a :: #run get_number_inner_type(A);
        b :: #run get_number_inner_type(B);
        T :: #run upcast(a, b);
        return Complex(T);
    } else #if is_A_complex {
        a :: #run get_number_inner_type(A);
        T :: #run upcast(a, B);
        return Complex(T);
    } else #if is_B_complex {
        b :: #run get_number_inner_type(B);
        T :: #run upcast(A,b);
        return Complex(T);

    } else {
        // ? the types in INTRINSIC_NUMBER_TYPES are sorted in order from mathematical superset to subset, e.g. ℂ > ℝ > ℤ > ℕ 
        for INTRINSIC_NUMBER_TYPES {
            if A == it || B == it then return it;
        }
        assert(false, "cannot upcast % -> %", A, B);
    }
}


UpCast :: struct(A: Type, B: Type) {
    T :: #run upcast(A, B);
}



is_minor_number_type :: inline ($major: Type, $minor: Type) -> bool {
    return upcast(major, minor) == major;
}






// ########  ########  #### ##    ## ########
// ##     ## ##     ##  ##  ###   ##    ##
// ##     ## ##     ##  ##  ####  ##    ##
// ########  ########   ##  ## ## ##    ##
// ##        ##   ##    ##  ##  ####    ##
// ##        ##    ##   ##  ##   ###    ##
// ##        ##     ## #### ##    ##    ##

#scope_module

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






COLOR :: enum {
    FG_RESET;
    FG_BLACK;
    FG_RED;
    FG_GREEN;
    FG_YELLOW;
    FG_BLUE;
    FG_MAGENTA;
    FG_CYAN;
    FG_WHITE;
}

ansi_color :: (c: COLOR = .FG_RESET) -> string {
    if c == {
        case .FG_RESET  ; return "\e[0m";
        case .FG_BLACK  ; return "\e[30m";
        case .FG_RED    ; return "\e[31m";
        case .FG_GREEN  ; return "\e[32m";
        case .FG_YELLOW ; return "\e[33m";
        case .FG_BLUE   ; return "\e[34m";
        case .FG_MAGENTA; return "\e[35m";
        case .FG_CYAN   ; return "\e[36m";
        case .FG_WHITE  ; return "\e[37m";
    }
    return ansi_color();
}





fg_color :: inline (txt: string, color: COLOR = .FG_RESET) -> string {
    return sprint("%1%2%3", ansi_color(color), txt, ansi_color());
}


println_push :: (text: string = "", args: ..Any, color: COLOR = .FG_RESET) #expand {
    c_text := ifx color != .FG_RESET then fg_color(text, color) else text;
    println(c_text, ..args);
    print_indent(1);
    `defer println("└╴");
    `defer print_indent(-1);
}

println :: (text: string = "", args: ..Any, color: COLOR = .FG_RESET) {

    c_text := ifx color != .FG_RESET then fg_color(text, color) else text;

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

println :: (arg: Any, color: COLOR = .FG_RESET) {
    text := sprint("%", arg);
    c_text := ifx color != .FG_RESET then fg_color(text, color) else text;
    
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