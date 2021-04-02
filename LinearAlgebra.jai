






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.1 2x2 Linear System Solve, page 153
// solves the system A x = b with 2x2 matrix A
sl2 :: solve_linear_2x2;
solve_linear_2x2 :: (A: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(A.rows == 2 && A.cols == 2, "A is not 2x2 matrix: % rows, % columns", A.rows, A.cols);
        assert(b.dim == 2, "b is not a 2-dim vector: dim = %", b.dim);
    }

    sigma := max(abs(A[0][0]), abs(A[0][1]), abs(A[1][0]), abs(A[1][1]));

    mu := A / sigma;
    det := mu[0][0] * mu[1][1] - mu[0][1] * mu[1][0];
    g := b/sigma;
    return cmat64(2, 2, mu[1][1], -mu[0][1], -mu[1][0], mu[0][0]) * g / det;
}









// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.2 Forward Linear Recurrence, page 155
// solve the system Lx = b where L is a unit left-triangluar matrix
// that means, all diagonals are 1 and all L_ij = 0 for j > i.
sllt :: solve_linear_left_triangular;
solve_linear_left_triangular :: (L: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(is_quadratic(L), "L is not quadratic! % rows, % columns", L.rows, L.cols);
        assert(is_diagonal_unit(L), "L is not unit: Diagonal entries are not 1");
        assert(is_left_triangular(L), "L is not left-triangular");
    }

    res := cvec64(b.dim);
    for i: 0..b.dim-1 {
        sum := CZERO;
        for j: 0..i-1 {
            sum += L[i][j]*res[j];
        }
        res.data[i] = b[i] - sum;
    }
    return res;
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.4 Backward Linear Recurrence, page 158
// solve the system Rx = b where R is a right-triangluar matrix
// that means all R_ij = 0 for j < i.
slrt :: solve_linear_right_triangular;
solve_linear_right_triangular :: (R: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(is_quadratic(R), "R is not quadratic! % rows, % columns", R.rows, R.cols);
        assert(is_right_triangular(R), "R is not right-triangular");
    }

    res := cvec64(b.dim);
    i : int;
    for 0..b.dim-1 {
        i = b.dim-1 - it;
        sum := CZERO;
        for j: i+1..b.dim-1 {
            sum += R[i][j]*res[j];
        }
        res.data[i] = (b[i] - sum)/R[i][i];
    }
    return res;
}






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.6 Solve Right Quasi-Triangular System, page 164
// solve the system Rx=b where R is a right(upper) quasi-triangular matrix.
slrqt :: solve_linear_right_quasi_triangular;
solve_linear_right_quasi_triangular :: (R: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(is_quadratic(R), "R is not quadratic! % rows, % columns", R.rows, R.cols);
        assert(is_right_quasi_triangular(R), "R is not right quasi-triangular!");
    }
    res := cvec64(b.dim);

    n := b.dim;
    i := n-1; // start with last row and wrap from bottom to up
    while i >= 0 {
        // check if 2x2 block by looking at left column
        // the term `|| i == 0` was missing
        if (i > 0 && R[i][i-1] == CZERO) || i == 0 {
            sum := CZERO;
            for j: i+1..n-1 {
                sum += R[i][j] * res[j];
            }
            res.data[i] = (b[i] - sum)/R[i][i];
            i -= 1;
        } else if i > 0 { // so now we have a 2x2 block where R[i][i] is the bottom-right element
            A := cmat64(2, 2, 
                R[i-1][i-1], R[i-1][i],
                R[i][i-1], R[i][i]
            );
            sum1 := CZERO;
            sum2 := CZERO;
            for j: i+1..n-1 {
                sum1 += R[i-1][j] * res[j];
                sum2 += R[i][j] * res[j];
            }
            c := cvec64(
                b[i-1] - sum1,
                b[i] - sum2
            );
            tmp := solve_linear_2x2(A, c);
            res.data[i-1] = tmp[0];
            res.data[i]   = tmp[1];
            i -= 2;
        }
    }
    return res;
}






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.8 Orthogonal Projection, page 167
solve_linear_unitary :: solve_linear_orthogonal_projection;
slop :: solve_linear_orthogonal_projection;
solve_linear_orthogonal_projection :: (Q: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(is_unitary(Q));
    }

    return dagger(Q) * b;
}


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.9 Successive Orthogonal Projection, page 168
// not sure if this algo is useful. seems to be slower than the previous one.
slsop :: solve_linear_successive_orthogonal_projection;
solve_linear_successive_orthogonal_projection :: (Q: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(is_unitary(Q));
    }
    bb := copy(b);

    res := cvec64(bb.dim);
    for i: 0..bb.dim-1 {
        q := row(Q, i);
        res.data[i] = conjugate(q) * bb; // not sure if this is performant...
        bb = bb - q * res[i];
    }

    return res;
}



solve_linear_left_trapezoidal :: (L: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(L.cols == b.dim);
        assert(is_left_trapezoidal(L));
    }

    L2 := cmat64_submatrix(L, 0, L.cols, 0, L.cols);
    b2 := cvec64_subvector(b, 0, L.cols);

    res := solve_linear_left_triangular(L2, b2);

    #if CHECKS {
        L3 := cmat64_submatrix(L, L.cols+1, L.rows, 0, L.cols);
        b3 := cvec64_subvector(b, L.cols+1, L.rows);
        assert( approx(L3 * res, b3) );
    }

    return res;
}

solve_linear_right_trapezoidal :: (R: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(R.cols == b.dim);
        assert(is_right_trapezoidal(R));
    }

    R2 := cmat64_submatrix(R, 0, R.rows, 0, R.rows);
    b2 := cvec64_subvector(b, 0, R.rows);

    res := solve_linear_right_triangular(R2, b2);

    #if CHECKS {
        R3 := cmat64_submatrix(R, 0, R.rows, R.rows+1, R.cols);
        b3 := cvec64_subvector(b, R.rows+1, R.cols);
        assert( approx(R3 * res, b3) );
    }

    return res;
}











test_linear_algebra :: () {
    print("\n\nLinear Algebra -------------------------------\n\n");

    {
        A := cmat64(2, 2, R, R, 2.0 - I, I);
        b := cvec64(I, c64(2.0, 4.0));
        res := solve_linear_2x2(A, b);
        print("%\n", res);
        assert(res == cvec64(
            -0.25 + 1.75*I,
            0.25 - 0.75*I
        ));
    }

    {
        L := cmat64(3, 3,
            R, CZERO, CZERO,
            3.0*R, R, CZERO,
            I, R+I, R
        );
        b := cvec64(3.0*R, -I, R-2.0*I);
        res := solve_linear_left_triangular(L, b);
        print("% = % ?\n", L * res, b);
        assert(L * res == b);
    }

    {
        M := cmat64(3, 3,
            2.0*R, I, R-2.0*I,
            CZERO, R, 3.0*I,
            CZERO, CZERO, -3.0*R
        );
        b := cvec64(3.0*R, -I, R-2.0*I);
        res := solve_linear_right_triangular(M, b);

        test := norm_2(M * res - b);
        print("%\n", test); 
        assert(test < 1.0e-10, "res = %, M * res = % =? b = %", res, M*res, b);
    }

    {
        m := mat64(4, 4,
            1, 3, 4, 5,
            0, 2, 3, 4,
            0, 1, 2, 3,
            0, 0, 0, 2
        );
        cm := cmat64(m);
        
        b := vec64(2,3,4,5);
        cb := cvec64(b);

        res := solve_linear_right_quasi_triangular(cm,cb);
        print("res = %\n", res);

        test := norm_2(cm * res - cb);
        assert(test < 1.0e-10);
    }
}