#import "Basic";
#import "Tagged_Union"(DEBUG=true);

Foo :: struct(T: Type) {
    data: T;
}


OPERATIONS :: enum_flags u32 {
    NONE    :: 0;
    LEVEL_1 :: (1 << 1);
    LEVEL_2 :: (1 << 2);
    LEVEL_3 :: (1 << 3);
    ADD_FOO_FOO :: (1 << 4) + xx LEVEL_1;
}


// Node :: Tagged_Union(Data, Operation);

// Data :: struct(T: Type) {
//     data: T;
// }

// data :: (value: $T) -> Node {
//     val: Data(T) = ---;
//     val.data = value;

//     res: Node;
//     set(*res, val);
//     return res;
// }

// Operation :: struct(T: Type) {
//     left: Node;
//     right: Node;
//     operation: OPERATIONS;
// }

// operation :: (op: OPERATIONS, left: Node, right: Node) -> Node {
//     restype :: left.T; // TODO insert function that determines major-type
//     val: Operation = ---;
//     val.left = left;
//     val.right = right;
//     val.operation = op;

//     res: Node;
//     set(*res, val);
//     return res;
// }


NODE_TYPE :: enum u8 {
    OPERATION;
    DATA;
}
Node :: struct {
    type: NODE_TYPE;
    left: *Node = null;
    right: *Node = null;
    operation: OPERATIONS = .NONE;
    data: Any;
    resType: Type;
}

operation :: (op: OPERATIONS, left: *Node, right: *Node, resType: Type) -> Node {
    res : Node;
    res.type = .OPERATION;
    res.left = left;
    res.right = right;
    res.operation = op;
    res.data = null;
    res.resType = resType;
    return res;
}

data :: (value: $T) -> Node {
    res : Node;
    res.data = value;
    res.type = .DATA;
    res.resType = T;
    return res;
}


Engine :: struct {}
SimpleEngine :: struct {
    using e: Engine;
}


add_foo_foo :: (se: SimpleEngine, a: Foo($T), b: Foo(T)) -> Foo(T) {
    res : Foo(T) = ---;
    res.data = a.data + b.data;
    return res;
}



// I can not define $LN/$RN as Node($L)/Node($R), so I have to distinguish via #if
eval :: inline (engine: $E/Engine, op: Node ) -> Any {

    if op.type == .OPERATION {
        if op.operation & .ADD_FOO_FOO > 0 {
            return add_foo_foo(engine, xx op.left.data, xx op.right.data);
        } else {
            res : T = ---;
            return res;
        }
    }
}



operator + :: (a: Foo($T), b: Foo(T)) -> Node {
    return operation(.ADD_FOO_FOO, *data(a), *data(b), Foo(T));
}


main :: () {
    engine : SimpleEngine;

    a := Foo(int).{42};
    b := Foo(int).{12};
    ops := a + b;
    print("%: %\n", ops.operation, ops);
    res := eval(engine, ops);
    print("res = %\n", res);
}