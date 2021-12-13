
#scope_module



// ######## ##    ## ########  ########  ######
//    ##     ##  ##  ##     ## ##       ##    ##
//    ##      ####   ##     ## ##       ##
//    ##       ##    ########  ######    ######
//    ##       ##    ##        ##             ##
//    ##       ##    ##        ##       ##    ##
//    ##       ##    ##        ########  ######


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

is_intrinsic_float_number :: inline (T: *Type_Info) -> bool {
    return T.type == .FLOAT;
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



// ? include Octonions later
is_commuting :: ($T: Type) -> bool {
    return !is_quaternion_number(T);
}
is_commuting :: (T: *Type_Info) -> bool {
    return !is_quaternion_number(T);
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







is_minor_number_type :: inline ($major: Type, $minor: Type) -> bool {
    // print("is_minor_number_type(%,%)\n", major, minor);
    T := upcast(major, minor);
    // print("is_minor_number_type: T = %\n", T);
    // print("T == major -> %\n", T == major);
    return T == major;
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
icast :: inline ($T: Type, num: Quaternion($N)) -> Quaternion(T) {
    res : Quaternion(T) = ---;
    res.r = cast(T)num.r;
    res.i = cast(T)num.i;
    res.j = cast(T)num.j;
    res.k = cast(T)num.k;
    return res;
}
icast :: inline ($T: Type, num: Complex($N)) -> Complex(T) {
    res : Complex(T) = ---;
    res.real = cast(T)num.real;
    res.imag = cast(T)num.imag;
    return res;
}
icast :: inline ($T: Type, num: $N) -> T #modify {
    if !is_intrinsic_number(N) {N = null; return;}
} {
    return cast(T)num;
}

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
    // print("upcast(%,%)\n", A, B);
    #if A == B {
        return A;
    } else {
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
                // print("checking intrinsic type % against (%,%)\n", it, A, B);
                if A == it || B == it {
                    // print(" > found intrinsic %\n", it);
                    return it;
                }
            }
            assert(false, "cannot upcast % & %", A, B);
        }
    }
}


UpCast :: struct(A: Type, B: Type) {
    T :: #run upcast(A, B);
}



