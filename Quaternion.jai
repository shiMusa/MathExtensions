




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


Rq     :: Quaternion(float64).{1.0, 0.0, 0.0, 0.0};
Iq     :: Quaternion(float64).{0.0, 1.0, 0.0, 0.0};
Jq     :: Quaternion(float64).{0.0, 0.0, 1.0, 0.0};
Kq     :: Quaternion(float64).{0.0, 0.0, 0.0, 1.0};
Q_ZERO :: Quaternion(float64).{0.0, 0.0, 0.0, 0.0};
Q_NAN  :: Quaternion(float64).{
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


qtn :: ($$r: $T = 0, $$i: T = 0, $$j: T = 0, $$k: T = 0) -> Quaternion(T) {
    #if is_constant(r)
     && is_constant(i)
     && is_constant(j)
     && is_constant(k) {
        return Quaternion(T).{r, i, j, k};
    } else {
        q : Quaternion(T) = ---;
        q.r = r;
        q.i = i;
        q.j = j;
        q.k = k;
        return q;
    }
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
        a.r + b.r,
        a.i + b.i,
        a.j + b.j,
        a.k + b.k
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
    b_abs := b.r*b.r + b.i*b.i + b.j*b.j + b.k*b.k;
    b_inv := qtn(b.r/b_abs,-b.i/b_abs,-b.j/b_abs,-b.k/b_abs);
    return a * b_inv;
}
operator / :: (a: Quaternion($A), b: Complex($B)) -> Quaternion(UpCast(A,B).T) {
    b_abs : UpCast(A,B).T = b.real*b.real + b.imag*b.imag;
    b_inv := cpx(b.real/b_abs,-b.imag/b_abs);
    return a * b_inv;
}
operator / :: (a: Complex($A), b: Quaternion($B)) -> Quaternion(UpCast(A,B).T) {
    b_abs := b.r*b.r + b.i*b.i + b.j*b.j + b.k*b.k;
    b_inv := qtn(b.r/b_abs,-b.i/b_abs,-b.j/b_abs,-b.k/b_abs);
    return a * b_inv;
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