Vector :: struct (T: Type, N: int) {
    data: [N]T;
}
dim :: inline (v: *Vector($T, $N)) -> int {
    return N;
}
get :: inline (v: *Vector($T, $N), i: int) -> T {
    #if CHECKS {
        assert(i < dim(v), "index out of dimensions: i = %, dim = %\n", i, dim(v));
    }
    return v.data[i];
}
set :: inline (v: *Vector($T, $N), i: int, val: T) {
    #if CHECKS {
        assert(i < dim(v), "index out of dimensions: i = %, dim = %\n", i, dim(v));
    }
    v.data[i] = val;
}


main :: () {
    print("\n\ntesting new vector definition\n\n");

    v : Vector(float64, 4) = ---;

    println("%", v);
    set(*v, 1, 42.0);

    println("%", v);

}



CHECKS :: true;


#scope_module

#import "Basic";
Math   :: #import "Math";
Pool   :: #import "Pool";
String :: #import "String";


pshift := 0;
println_shift :: inline (levels: int) {
    pshift += 4 * levels;
}

println :: inline (msg: string, args: ..Any) {
    for 0..pshift {
        print(" ");
    }
    print(msg, ..args);
    print("\n");
}

println :: inline (arg: Any) {
    for 0..pshift {
        print(" ");
    }
    print("%\n", arg);
}