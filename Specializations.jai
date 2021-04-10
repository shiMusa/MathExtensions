






#insert #run specialize("_equal", "operator ==", COMBI_M_M_b);
#insert #run specialize("_equal", "operator ==", COMBI_V_V_b);
_equal :: inline (a: $A, b: $B, res: *bool) {
    #if #run is_a_matrix(A) {
        if a.rows != b.rows || a.cols != b.cols {
                <<res = false;
                return;
        }
        for a.data {
            if it != b.data[it_index] {
                <<res = false;
                return;
            }
        }
        <<res = true;
    } else #if #run is_a_vector(A) {
        if a.dim != b.dim {
            <<res = false;
        }
        for a.data {
            if it != b.data[it_index] {
                <<res = false;
                return;
            }
        }
        <<res = true;
    } else {
        print("types not considered.")
        <<res = false;
        assert(false);
    }
}









#insert #run specialize("_mul", "operator *", COMBI_M_M_M);
#insert #run specialize("_mul", "operator *", COMBI_M_V_V);
#insert #run specialize("_mul", "operator *", COMBI_V_M_V);
#insert #run specialize("_mul", "operator *", COMBI_M_X_M, .["#symmetric"]);
#insert #run specialize("_mul", "operator *", COMBI_V_V_X);
#insert #run specialize("_mul", "operator *", COMBI_V_X_V, .["#symmetric"]);
_mul :: (a: $A, b: $B, res: *$C) {

    // M * M -> M
    //
    // Scientific Computing, Vol I: Linear and nonlinear equations
    // Texts in computational science and engineering 18
    // Springer
    // 2.11 Case Study: Matrix-Matrix Multiplication, page 118
    // According to Table 2.12, the best performing loop order is j/k/i with
    // the least number of cache misses.
    #if #run is_a_matrix(A) && is_a_matrix(B) && is_a_matrix(C) {
        // print("M_M_M\n");
        #if CHECKS {
            assert(a.cols == b.rows, "Matrix dimensions don't match: a.cols = %, b.rows = %", a.cols, b.rows);
        }

        <<res = mat_zeros(C.T, a.rows, b.cols);
        for j: 0..b.cols-1 {
            for k: 0..a.cols-1 {
                tmp := ncast(C.T, b[k][j]);
                for i: 0..a.rows-1 {
                    (<<res)[i][j] += ncast(C.T, a[i][k]) * tmp;
                }
            }
        }

    // M * V -> V
    } else #if #run is_a_matrix(A) && is_a_vector(B) && is_a_vector(C) {
        // print("M_V_V\n");
        #if CHECKS {
            assert(a.cols == b.dim, "Dimensions mismatch: a.cols = %, b.dim = %", a.cols, b.dim);
        }

        <<res = vec_zeros(C.T, a.rows);
        for i: 0..a.rows-1 {
            for j: 0..a.cols-1 {
                res.data[i] += ncast(C.T, a[i][j]) * ncast(C.T, b[j]);
            }
        }

    // V * M -> V
    } else #if #run is_a_vector(A) && is_a_matrix(B) && is_a_vector(C) {
        // print("V_M_V\n");
        #if CHECKS {
            assert(a.dim == b.rows, "Dimensions mismatch: a.dim = %, b.rows = %", a.dim, b.rows);
        }

        <<res = vec_zeros(C.T, b.cols);
        for j: 0..b.rows-1 {
            for i: 0..a.dim-1 {
                res.data[j] += ncast(C.T, a[i]) * ncast(C.T, b[i][j]);
            }
        }

    // M * X -> M
    } else #if #run is_a_matrix(A) && is_a_scalar(B) && is_a_matrix(C) {
        // print("M_X_M\n");
        <<res = mat(C.T, a.rows, a.cols);
        for 0..a.data.count-1 {
            res.data[it] = ncast(C.T, a.data[it]) * ncast(C.T, b);
        }

    // V * V -> X
    } else #if #run is_a_vector(A) && is_a_vector(B) && is_a_scalar(C) {
        <<res = zero(C);
        for 0..a.dim-1 {
            <<res == ncast(C, a[it] * b[it]);
        }

    // V * X -> V
    } else #if #run is_a_vector(A) && is_a_scalar(B) && is_a_vector(C) {
        <<res = vec(C.T, a.dim);
        for 0..a.dim-1 {
            res.data[it] = ncast(C.T, a[it]) * ncast(C.T, b);
        }
    }
    
    else {
        print("no matrix operation found!\na = %\n, b = %\n, res = %\n", a, b, <<res);
    }
}










// I have to write it this way, because even the M * X -> M case needs an generic
// parameter $B for float64/Complex64/... but that can be interpreted as any
// other type such as Matrix(T)/Vector(T) etc. as well.
// I would have to create a struc Scalar(T) but that seems kinda useless and
// overly complicated?
// April 10th, 2021 - Felix Fehse
#insert #run specialize("_div", "operator /", COMBI_M_M_M);
#insert #run specialize("_div", "operator /", COMBI_V_M_V);
#insert #run specialize("_div", "operator /", COMBI_M_X_M);
#insert #run specialize("_div", "operator /", COMBI_V_V_M);
#insert #run specialize("_div", "operator /", COMBI_V_X_V, .["#symmetric"]);
_div :: (a: $A, b: $B, res: *$C) {

    // M / M -> M
    #if #run is_a_matrix(A) && is_a_matrix(B) && is_a_matrix(C) {
        // print("M_M_M\n");
        <<res = inverse(b) * a;

    // V / M -> V
    } else #if #run is_a_vector(A) && is_a_matrix(B) && is_a_vector(C) {
        // print("V_M_V\n");
        <<res = inverse(b) * a;

    // M / X -> M
    } else #if #run is_a_matrix(A) && is_a_scalar(B) && is_a_matrix(C) {
        // print("M_X_M\n");
        assert(C == Matrix);
        <<res = mat(C.T, a.rows, a.cols);
        for 0..a.data.count-1 {
            res.data[it] = a.data[it] / b;
        }

    // V / V -> M
    } else #if #run is_a_vector(A) && is_a_vector(B) && is_a_scalar(C) {
        <<res == outer_product(a, b);

    // V / X -> V
    } else #if #run is_a_vector(A) && is_a_scalar(B) && is_a_vector(C) {
        <<res = vec(C.T, a.dim);
        for 0..a.dim-1 {
            res.data[it] = a[it] / b;
        }
    }
    
    else {
        print("no operator/ found!\na = %\n, b = %\n, res = %\n", a, b, <<res);
    }
}













#insert #run specialize("_outer_product", "outer_product", COMBI_V_V_M);
_outer_product :: (a: Vector($A), b: Vector($B), res: *Matrix($C)) {
    <<res = mat(C, a.dim, b.dim);
    for i: 0..a.dim-1 {
        for j: 0..b.dim-1 {
            (<<res)[i][j] = a[i] * b[j];
        }
    }
}











#insert #run specialize("_cross", "cross", COMBI_V_V_V);
_cross :: (a: Vector($A), b: Vector($B), res: *Vector($C)) {
    #if CHECKS {
        assert(a.dim == 3 && b.dim == 3);
    }
    <<res = vec(C, 
        a[1]*b[2] - a[2]*b[1], 
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}










// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.10 Apply Elementary Reflector, page 170
#insert #run specialize("_reflect", "reflect", COMBI_V_V_V);
_reflect :: (vec: Vector($A), normal: Vector($B), res: *Vector($C)) {
    tau := 2.0 / (conjugate(normal) * normal);
    alpha := tau * conjugate(normal) * vec;
    <<res = vec - normal * alpha;
}












#insert #run specialize("_add", "operator +", COMBI_M_M_M);
#insert #run specialize("_add", "operator +", COMBI_V_V_V);
_add :: (a: $A, b: $B, res: *$C) {

    // M + M -> M
    #if #run is_a_matrix(A) && is_a_matrix(B) && is_a_matrix(C) {
        #if CHECKS {
            assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
        }

        <<res = mat(C.T, a.rows, a.cols);
        for 0..a.data.count-1 {
            res.data[it] = a.data[it] + b.data[it];
        }
    
    // V + V -> V
    } else #if #run is_a_vector(A) && is_a_vector(B) && is_a_vector(C) {
        #if CHECKS {
            assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
        }
        <<res = vec(C.T, a.dim);
        for 0..a.dim-1 {
            res.data[it] = a[it] + b[it];
        }
    }
}













#insert #run specialize("_sub", "operator -", COMBI_M_M_M);
#insert #run specialize("_sub", "operator -", COMBI_V_V_V);
_sub :: (a: $A, b: $B, res: *$C) {

    // M + M -> M
    #if #run is_a_matrix(A) && is_a_matrix(B) && is_a_matrix(C) {
        #if CHECKS {
            assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
        }

        <<res = mat(C.T, a.rows, a.cols);
        for 0..a.data.count-1 {
            res.data[it] = a.data[it] - b.data[it];
        }
    
    // V + V -> V
    } else #if #run is_a_vector(A) && is_a_vector(B) && is_a_vector(C) {
        #if CHECKS {
            assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
        }
        <<res = vec(C.T, a.dim);
        for 0..a.dim-1 {
            res.data[it] = a[it] - b[it];
        }
    }
}











#insert #run specialize("_neg", "operator -", COMBI_M_M);
#insert #run specialize("_neg", "operator -", COMBI_V_V);
_neg :: (a: $A, res: *$C) {
    // M -> M
    #if #run is_a_matrix(A) && is_a_matrix(C) {
        #if CHECKS {
            assert(a.rows == b.rows && a.cols == b.cols, "Matrices have different dimensions: a: %x%, b: %x%", a.rows, a.cols, b.rows, b.cols);
        }

        <<res = mat(C.T, a.rows, a.cols);
        for 0..a.data.count-1 {
            res.data[it] = -a.data[it];
        }
    
    // V -> V
    } else #if #run is_a_vector(A) && is_a_vector(C) {
        #if CHECKS {
            assert(a.dim == b.dim, "Vectors have different dimensions: a.dim = %, b.dim = %", a.dim, b.dim);
        }
        <<res = vec(C.T, a.dim);
        for 0..a.dim-1 {
            res.data[it] = -a[it];
        }
    }
}