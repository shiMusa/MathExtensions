
#scope_module



// ######## ##    ## ########  ########  ######
//    ##     ##  ##  ##     ## ##       ##    ##
//    ##      ####   ##     ## ##       ##
//    ##       ##    ########  ######    ######
//    ##       ##    ##        ##             ##
//    ##       ##    ##        ##       ##    ##
//    ##       ##    ##        ########  ######


INTRINSIC_NUMBER_TYPES :: Type.[
    float64, float32, s64, s32, s16, s8, u64, u32, u16, u8, bool, int, float
];

INTRINSIC_FLOATING_NUMBER_TYPES :: Type.[
    float64, float32, float
];



is_float_number :: (T: Type) -> bool {
    if is_complex_number(T) || is_quaternion_number(T) {
        X := get_number_inner_type(T);
        return is_float_number(X);
    } 

    for INTRINSIC_FLOATING_NUMBER_TYPES
        if T == it return true;

    return false;
}


is_intrinsic_number :: (T: Type) -> bool {
    for INTRINSIC_NUMBER_TYPES {
        if T == it then return true;
    }
    return false;
}


is_non_real_number :: inline (T: Type) -> bool {
    return is_complex_number(T) || is_quaternion_number(T);
}



complex_with_inner_type :: (T: Type) -> Type {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == % then return Complex(%);\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    // dummy
    return Complex(float64);
}



is_complex_number :: (T: Type) -> bool {

    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == Complex(%) then return true;\n",
                it
            );
        }

        return builder_to_string(*builder);
    }

    return false;
}


quaternion_with_inner_type :: (T: Type) -> Type {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == % then return Quaternion(%);\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    // dummy
    return Quaternion(float64);
}


is_quaternion_number :: (T: Type) -> bool {

    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == Quaternion(%) then return true;\n",
                it
            );
        }

        return builder_to_string(*builder);
    }

    return false;
}



// ? include Octonions later
is_commuting :: (T: Type) -> bool {
    return !is_quaternion_number(T);
}



is_number :: inline (T: Type) -> bool {
    return is_intrinsic_number(T) || is_complex_number(T) || is_quaternion_number(T);
}




get_number_inner_type :: (T: Type) -> Type {
    if is_complex_number(T) {
        return get_complex_inner_type(T);
    } else if is_quaternion_number(T) {
        return get_quaternion_inner_type(T);
    } else {
        return T;
    }
}

// ? assuming we already know it's a Complex(T) number
get_complex_inner_type :: (T: Type) -> Type, bool {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == Complex(%) then return %, true;\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return T, false;
}

get_complex_type_info_from_intrinsic_type_info :: (T: Type) -> Type, bool {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if % == T then return Complex(%), true;\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return T, false;
}





// ? assuming we already know it's a Complex(T) number
get_quaternion_inner_type :: (T: Type) -> Type, bool {
    #insert -> string {
        builder: String_Builder;

        for INTRINSIC_NUMBER_TYPES {
            print_to_builder(*builder,
                "if T == Quaternion(%) then return %, true;\n",
                it, it
            );
        }

        return builder_to_string(*builder);
    }
    return T, false;
}





is_minor_number_type :: inline (major: Type, minor: Type) -> bool {
    T := upcast(major, minor);
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
    return is_intrinsic_number(N);
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
            #run assert(false, "cannot ncast % -> %", N, T);
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
            #run assert(false, "cannot ncast % -> %", N, T);
        } else {
            
            return cast(T)num;

        }
    }
}

zero :: #bake_arguments ncast(num=0);
one  :: #bake_arguments ncast(num=1);

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
upcast :: (A: Type, B: Type) -> Type {
    // print("upcast(%,%)\n", A, B);
    if A == B {
        return A;
    } else {
        is_A_quaternion := is_quaternion_number(A);
        is_B_quaternion := is_quaternion_number(B);
        is_A_complex    := is_complex_number(A);
        is_B_complex    := is_complex_number(B);

        if is_A_quaternion && is_B_quaternion {
            a := get_number_inner_type(A);
            b := get_number_inner_type(B);
            T := upcast(a, b);
            return quaternion_with_inner_type(T);
        } else if is_A_quaternion {
            a := get_number_inner_type(A);
            T := upcast(a, B);
            return quaternion_with_inner_type(T);
        } else if is_B_quaternion {
            b := get_number_inner_type(B);
            T := upcast(A,b);
            return quaternion_with_inner_type(T);

        } else if is_A_complex && is_B_complex {
            a := get_number_inner_type(A);
            b := get_number_inner_type(B);
            T := upcast(a, b);
            return complex_with_inner_type(T);
        } else if is_A_complex {
            a := get_number_inner_type(A);
            T := upcast(a, B);
            return complex_with_inner_type(T);
        } else if is_B_complex {
            b := get_number_inner_type(B);
            T := upcast(A,b);
            return complex_with_inner_type(T);

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
            return int; // ? to get rid of warnings.
        }
    }
}


UpCast :: struct(A: Type, B: Type) {
    T :: #run upcast(A, B);
}



