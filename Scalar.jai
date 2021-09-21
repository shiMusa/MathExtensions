


Scalar :: struct(T: Type) {
    value: T;
}

scl :: inline (val: $T) -> Scalar(T) {
    res : Scalar(T) = ---;
    res.value = val;
    return res;
}


operator + :: inline(a: Scalar($A), b: Scalar($B)) -> Scalar($C) #modify { C = upcast(A,B); } {
    return scl(a.value + b.value);
}
operator + :: inline(a: Scalar($A), b: $B) -> Scalar($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A,B);
} #symmetric {
    return scl(a.value + b);
}

operator - :: inline(a: Scalar($A), b: Scalar($B)) -> Scalar($C) #modify { C = upcast(A,B); } {
    return scl(a.value - b.value);
}
operator - :: inline(a: Scalar($A), b: $B) -> Scalar($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A,B);
} {
    return scl(a.value - b);
}
operator - :: inline(b: $B, a: Scalar($A)) -> Scalar($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A,B);
} {
    return scl(b - a.value);
}

operator * :: inline(a: Scalar($A), b: Scalar($B)) -> Scalar($C) #modify { C = upcast(A,B); } {
    return scl(a.value * b.value);
}
operator * :: inline(a: Scalar($A), b: $B) -> Scalar($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A,B);
    // print("operator * (Scalar(%), %) -> %\n", A.type, B.type, C);
} #symmetric {
    return scl(a.value * b);
}

operator / :: inline(a: Scalar($A), b: Scalar($B)) -> Scalar($C) #modify { C = upcast(A,B); } {
    return scl(a.value / b.value);
}
operator / :: inline(a: Scalar($A), b: $B) -> Scalar($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A,B);
} {
    return scl(a.value/b);
}
operator / :: inline(b: $B, a: Scalar($A)) -> Scalar($C) #modify {
    if !is_number(B) { C = null; return; }
    C = upcast(A,B);
} {
    return scl(b/a.value);
}

operator - :: inline(a: Scalar($A)) -> Scalar(A) {
    return scl(-a.value);
}


operator == :: inline(a: Scalar($A), b: Scalar($B)) -> bool {
    return a.value == b.value;
}
operator == :: inline(a: Scalar($A), b: $B) -> bool #modify {
    if !is_number(B) { return; }
} #symmetric {
    return a.value == b;
}




////////////////////////////////////////////////////////////////////////////////////////
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
test_scalar :: () {

    println_push("Scalar", color = .FG_WHITE);
    {
        s := scl(2.0);
        r := s * 42.0;
        assert(r == 84.0);
    }
}