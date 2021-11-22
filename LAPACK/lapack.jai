
#load "xerbla.jai";
#load "Level1/_axpy.jai";
#load "Level1/_copy.jai";
#load "Level1/_dotc.jai";
#load "Level1/_dotu.jai";
#load "Level1/_scal.jai";
#load "Level1/_rot.jai";
#load "Level1/_swap.jai";
#load "Level1/_asum.jai";
#load "Level1/_abs.jai";
#load "Level1/drotm.jai";
#load "Level1/drotmg.jai";

#load "Level2/_trsv.jai";

// ? tests if val equals any of lower_case_chars
char_eq :: (val: s64, $lower_case_chars: ..s64) -> bool {
    #if lower_case_chars.count == 1 {
        return val == lower_case_chars[0] || val == lower_case_chars[0]-32;
    } else {
        for lower_case_chars {
            if val == it || val == it-32 then return true;
        }
        return false;
    }
}


test_lapack :: () {
    println_push("LAPACK", color = .FG_WHITE);

    {
        println_push("char_eq", color = .FG_GREEN);
        val := #char "a";
        VAL := #char "A";
        println("a = %", val);
        println("A = %", VAL);
        println("\"a\" == \"a\" == \"A\": %, %", char_eq(val, #char"a"), char_eq(VAL, #char"a"));
    }

    test_lapack_axpy();
    test_lapack_copy();
    test_lapack_dotc();
    test_lapack_dotu();
    test_lapack_scal();
    test_lapack_rot();
    test_lapack_swap();
    test_lapack_asum();
}