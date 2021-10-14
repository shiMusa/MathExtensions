




str :: (q: Quaternion($T)) -> string {
    builder : String_Builder;
    defer free_buffers(*builder);
    init_string_builder(*builder);

    print_to_builder(*builder, "⟨%", q.r);

    if q.i >= 0 {
        print_to_builder(*builder, " + %i", q.i);
    } else {
        print_to_builder(*builder, " - %i", -q.i);
    }

    if q.j >= 0 {
        print_to_builder(*builder, " + %j", q.j);
    } else {
        print_to_builder(*builder, " - %j", -q.j);
    }

    if q.k >= 0 {
        print_to_builder(*builder, " + %k", q.k);
    } else {
        print_to_builder(*builder, " - %k", -q.k);
    }

    append(*builder, "⟩");

    return builder_to_string(*builder);
}


qR     :: Q64.{1.0, 0.0, 0.0, 0.0};
qI     :: Q64.{0.0, 1.0, 0.0, 0.0};
qJ     :: Q64.{0.0, 0.0, 1.0, 0.0};
qK     :: Q64.{0.0, 0.0, 0.0, 1.0};
Q_ZERO :: Q64.{0.0, 0.0, 0.0, 0.0};
Q_NAN  :: Q64.{
    Math.FLOAT64_NAN,
    Math.FLOAT64_NAN,
    Math.FLOAT64_NAN,
    Math.FLOAT64_NAN
};


Q64 :: Quaternion(float64);
Q32 :: Quaternion(float32);
Quaternion :: struct (T: Type) {
    r: T;
    i: T;
    j: T;
    k: T;
}


q64 :: inline ($$r: $T = 0.0, $$i: $U = 0.0, $$j: $V = 0.0, $$k: $W = 0.0) -> Q64 #modify {
    if !is_intrinsic_number(T)
    || !is_intrinsic_number(U)
    || !is_intrinsic_number(V)
    || !is_intrinsic_number(W) {
        T = null;
        U = null;
        V = null;
        W = null;
        return;
    }
} {
    return qtn(
        ncast(float64, r),
        ncast(float64, i),
        ncast(float64, j),
        ncast(float64, k)
    );
}


q32 :: inline ($$r: $T = 0.0, $$i: $U = 0.0, $$j: $V = 0.0, $$k: $W = 0.0) -> Q64 #modify {
    if !is_intrinsic_number(T)
    || !is_intrinsic_number(U)
    || !is_intrinsic_number(V)
    || !is_intrinsic_number(W) {
        T = null;
        U = null;
        V = null;
        W = null;
        return;
    }
} {
    return qtn(
        ncast(float64, r),
        ncast(float64, i),
        ncast(float64, j),
        ncast(float64, k)
    );
}



qtn :: inline ($$r : $T = 0, $$v: $V/VectorType) -> UpCast(T, V.data_type).T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    #run assert(is_intrinsic_number(V.data_type), "cannot use non-real numbers for this method: V.data_type = %\n", V.data_type);
    
    X :: UpCast(T, V.data_type).T;
    #if is_constant(r) && is_constant(v) {
        return Quaternion(X).{r, get(v, 0), get(v, 1), get(v, 2)};
    } else {
        res : Quaternion(X) = ---;
        res.r = r;
        res.i = get(v, 0);
        res.j = get(v, 1);
        res.k = get(v, 2);
        return res;
    }
}


qtn :: inline ($$r: $T = 0, $$i: $U = 0, $$j: $V = 0, $$k: $W = 0) -> Quaternion(T) #modify {
    if !is_intrinsic_number(T)
    || !is_intrinsic_number(U)
    || !is_intrinsic_number(V)
    || !is_intrinsic_number(W) {
        T = null;
        U = null;
        V = null;
        W = null;
        return;
    }
} {
    #if is_constant(r) 
     && is_constant(i)
     && is_constant(j)
     && is_constant(k) {
        return Quaternion(T).{
            r,
            #run ncast(T, i),
            #run ncast(T, j),
            #run ncast(T, k)
        };
    } else {
        res : Quaternion(T) = ---;
        res.r = r;
        res.i = ncast(T, i);
        res.j = ncast(T, j);
        res.k = ncast(T, k);
        return res;
    }
}

qtn_r :: inline ($T: Type) -> Quaternion(T) {
    #run assert(is_intrinsic_number(T), "can only use intrinsic numbers to generate quaternion: Quaternion(T), T = % does not work.", T);
    return qtn(one(T), 0, 0, 0);
}
qtn_i :: inline ($T: Type) -> Quaternion(T) {
    #run assert(is_intrinsic_number(T), "can only use intrinsic numbers to generate quaternion: Quaternion(T), T = % does not work.", T);
    return qtn(zero(T), one(T), zero(T), zero(T));
}
qtn_j :: inline ($T: Type) -> Quaternion(T) {
    #run assert(is_intrinsic_number(T), "can only use intrinsic numbers to generate quaternion: Quaternion(T), T = % does not work.", T);
    return qtn(zero(T), 0, 1, 0);
}
qtn_k :: inline ($T: Type) -> Quaternion(T) {
    #run assert(is_intrinsic_number(T), "can only use intrinsic numbers to generate quaternion: Quaternion(T), T = % does not work.", T);
    return qtn(zero(T), 0, 0, 1);
}

qtn_vec_is_zero :: inline (q: Quaternion($T)) -> bool {
    return q.i == zero(T) && q.j == zero(T) && q.k == zero(T);
}

qtn_vec :: inline (q: Quaternion($T)) -> Quaternion(T) {
    return qtn(zero(T), q.i, q.j, q.k);
}

qtn_vec_abs_sq :: inline (q: Quaternion($T)) -> T {
    return q.i*q.i + q.j*q.j + q.k*q.k;
}

qtn_vec_abs :: inline (q: Quaternion($T)) -> T {
    return Math.sqrt(qtn_vec_abs_sq(q));
}

qtn_inv :: inline (q: Quaternion($T)) -> Quaternion(T) {
    q_abs_sq := abs_sq(q);
    return qtn(q.r/q_abs_sq,-q.i/q_abs_sq,-q.j/q_abs_sq,-q.k/q_abs_sq);
}







qtn_rot :: inline (angle: $U, axis: $T/VectorType, $normalize: bool = false) -> UpCast(T.data_type, U).T #modify {
    if !is_intrinsic_number(U) {U = null; return;}    
} {
    #run assert(is_intrinsic_number(T.data_type), "this method only works for real numbers, but axis.data_type = %\n", T.data_type);
    #if normalize {
        ax := axis/abs_sq(axis);
        return qtn_rot(angle, ax);
    } else {
        return qtn(Math.cos(angle/2), axis*Math.sin(angle/2));
    }
}

mat_3_from_qtn :: (q: Quaternion($T)) -> DenseMatrix(T, 3, 3) {
    res : DenseMatrix(T, 3, 3);
    s := one(T)/abs_sq(q); // ! is this correct? Double check!
    
    mat_init(*res, 
        1-2*s*(q.j*q.j+q.k*q.k),   2*s*(q.i*q.j-q.k*q.r),   2*s*(q.i*q.k+q.j*q.r),
          2*s*(q.i*q.j+q.k*q.r), 1-2*s*(q.i*q.i+q.k*q.k),   2*s*(q.j*q.k-q.i*q.r),
          2*s*(q.i*q.k-q.j*q.r),   2*s*(q.j*q.k+q.i*q.k), 1-2*s*(q.i*q.i+q.j*q.j)
    );

    return res;
}





//  #######  ########  ######## ########     ###    ########  #######  ########
// ##     ## ##     ## ##       ##     ##   ## ##      ##    ##     ## ##     ##
// ##     ## ##     ## ##       ##     ##  ##   ##     ##    ##     ## ##     ##
// ##     ## ########  ######   ########  ##     ##    ##    ##     ## ########
// ##     ## ##        ##       ##   ##   #########    ##    ##     ## ##   ##
// ##     ## ##        ##       ##    ##  ##     ##    ##    ##     ## ##    ##
//  #######  ##        ######## ##     ## ##     ##    ##     #######  ##     ##



//   ##
//   ##
// ###### #######
//   ##
//   ##



operator + :: inline (a: Quaternion($A), b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) {
    return qtn(
        a.r + b.r,
        a.i + b.i,
        a.j + b.j,
        a.k + b.k
    );
}
operator - :: inline (a: Quaternion($A), b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) {
    return qtn(
        a.r - b.r,
        a.i - b.i,
        a.j - b.j,
        a.k - b.k
    );
}


operator + :: inline (a: Quaternion($A), b: Complex($B)) -> Quaternion(UpCast(A,B).T) #symmetric {
    return qtn(
        a.r + b.real,
        a.i + b.imag,
        a.j,
        a.k
    );
}
operator - :: inline (a: Quaternion($A), b: Complex($B)) -> Quaternion(UpCast(A,B).T) {
    return qtn(
        a.r - b.real,
        a.i - b.imag,
        a.j,
        a.k
    );
}
operator - :: inline (a: Complex($A), b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) {
    return qtn(
        b.real - a.r,
        b.imag - a.i,
        - a.j,
        - a.k
    );
}

operator + :: inline (a: Quaternion($A), b: $B) -> Quaternion(UpCast(A,B).T) #modify {
    if !is_intrinsic_number(B) {B = null; return;}
} #symmetric {
    return qtn(
        a.r + b,
        a.i,
        a.j,
        a.k
    );
}
operator - :: inline (a: Quaternion($A), b: $B) -> Quaternion(UpCast(A,B).T) #modify {
    if !is_intrinsic_number(B) {B = null; return;}
} {
    return qtn(
        a.r - b,
        a.i,
        a.j,
        a.k
    );
}
operator - :: inline (a: $A, b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) #modify {
    if !is_intrinsic_number(A) {A = null; return;}
} {
    return qtn(
        a - b.r,
          - b.i,
          - b.j,
          - b.k
    );
}


operator - :: inline (q: Quaternion($T)) -> Quaternion(T) {
    return qtn(-q.r, -q.i, -q.j, -q.k);
}


//                 ##
//  ##   ##       ##
//   ## ##       ##
// #########    ##
//   ## ##     ##
//  ##   ##   ##
//           ##



operator * :: inline (q: Quaternion($T), v: $V) -> Quaternion(UpCast(T,V).T) #modify {
    if !is_intrinsic_number(V) {V = null; return;}    
} #symmetric {
    return qtn(
        q.r * v,
        q.i * v,
        q.j * v,
        q.k * v
    );
}
operator / :: inline (q: Quaternion($T), v: $V) -> Quaternion(UpCast(T,V).T) #modify {
    if !is_intrinsic_number(V) {V = null; return;}    
} {
    return qtn(
        q.r / v,
        q.i / v,
        q.j / v,
        q.k / v
    );
}
// wikipedia, https://www.wikiwand.com/en/Quaternion
operator / :: inline (v: $V, q: Quaternion($T)) -> Quaternion(UpCast(T,V).T) #modify {
    if !is_intrinsic_number(V) {V = null; return;}    
} {
    abs := q.r*q.r + q.i*q.i + q.j*q.j + q.k*q.k;
    return qtn(
         v * q.r / abs,
        -v * q.i / abs,
        -v * q.j / abs,
        -v * q.k / abs
    );
}



operator * :: inline (a: Quaternion($A), b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) {
    return qtn(
        a.r*b.r - a.i*b.i - a.j*b.j - a.k*b.k,
        a.r*b.i + a.i*b.r + a.j*b.k - a.k*b.j,
        a.r*b.j - a.i*b.k + a.j*b.r + a.k*b.i,
        a.r*b.k + a.i*b.j - a.j*b.i + a.k*b.r
    );
}
operator * :: inline (a: Quaternion($A), b: Complex($B)) -> Quaternion(UpCast(A,B).T) #symmetric {
    return qtn(
        a.r*b.real - a.i*b.imag,
        a.r*b.imag + a.i*b.real,
        a.j*b.real + a.k*b.imag,
        a.j*b.imag + a.k*b.real
    );
}


operator / :: (a: Quaternion($A), b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) {
    return a * qtn_inv(b);
}
operator / :: (a: Quaternion($A), b: Complex($B)) -> Quaternion(UpCast(A,B).T) {
    return a * cpx_inv(b);
}
operator / :: (a: Complex($A), b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) {
    return a * qtn_inv(b);
}


// ##### #####

// ##### #####


operator == :: inline (a: Quaternion($T), b: Quaternion($U)) -> bool {
    // TODO UpCast here?
    return a.r == b.r && a.i == b.i && a.j == b.j && a.k == b.k;
}

operator == :: inline (a: Quaternion($T), b: Complex($U)) -> bool #symmetric {
    // TODO UpCast here?
    return a.r == b.real && a.i == b.imag && a.j == 0 && a.k == 0;
}

operator == :: inline (a: Quaternion($T), b: $U) -> bool #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} #symmetric {
    return a.r == b && a.i == 0 && a.j == 0 && a.k == 0;
}




// ######## ########  ######  ########  ######
//    ##    ##       ##    ##    ##    ##    ##
//    ##    ##       ##          ##    ##
//    ##    ######    ######     ##     ######
//    ##    ##             ##    ##          ##
//    ##    ##       ##    ##    ##    ##    ##
//    ##    ########  ######     ##     ######

#scope_module
test_quaternion :: () {
    println_push("Quaternion", color = .FG_WHITE);

    a := q64(1.0, 0.0, 0.0, 0.0);
    println("a = %", str(a));
    b := q64(0.0, 1.0, 0.0, 0.0);
    println("b = %", str(b));
    c := q64(0.0, 0.0, 1.0, 0.0);
    println("c = %", str(c));
    d := q64(0.0, 0.0, 0.0, 1.0);
    println("d = %", str(d));
    
    assert( a != b );

    println("testing negative zero");
    assert( q64(0.0, 0.0, 0.0, 0.0) == q64(-0.0, -0.0, -0.0, -0.0) );

    println("testing + - operations");
    assert( 3 + 2*I + q64(2, 0, -6, 0) - 3 * q64(1, 3, -2, 0) == q64(2, -7, 0, 0) );


    {
        println_push("testing * / operations", color = .FG_GREEN);

        println("interaction with complex numbers");
        assert( b * b     == -1.0 * R       );
        assert( a * a     ==  1.0 * R       );
        assert( a/b       == -1.0 * I       );
        assert( 2.0*a+b   == c64(2.0,  1.0) );
        assert( 2.0*a-3*b == c64(2.0, -3.0) );

        println("quaternions");
        println("% : %", str(qtn(-16, 0, -1, 25)), type_of(qtn(-16, 0, -1, 25)) );
        assert( qtn(2, 0, -6, 3) * qtn(1, 3, -2, 2) == qtn(-16, 0, -1, 25) );    // ? non-commutative!
        assert( qtn(1, 3, -2, 2) * qtn(2, 0, -6, 3) == qtn(-16, 12, -19, -11) ); // ? non-commutative! 

        x := q64(2, 0, -6, 3);
        y := q64(1, 3, -2, 2);
        z := q64(10.0/9.0, 0, -11.0/18.0, -19.0/18.0);
        println("% / %", str(x), str(y));
        println("    = %", str(x/y));
        println("    ? %", str(z));
        // assert( x / y == z );
    }
}