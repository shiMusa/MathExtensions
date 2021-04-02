






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
sllta :: solve_linear_left_triangular;
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
slrta :: solve_linear_right_triangular;
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
slrqta :: solve_linear_right_quasi_triangular;
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


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 3.4.5.1 Left Trapezoidal, page 159
slltz :: solve_linear_left_trapezoidal;
solve_linear_left_trapezoidal :: (L: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(L.rows == b.dim);
        assert(is_left_trapezoidal(L));
    }

    L2 := cmat64_submatrix(L, 0, L.cols, 0, L.cols);
    b2 := cvec64_subvector(b, 0, L.cols);

    // print("given \n    L = %\n    b = %\n", L, b);
    // print("supsystem \n    L2 = %\n    b2 = %\n", L2, b2);

    res := solve_linear_left_triangular(L2, b2);

    #if CHECKS {
        if L.cols == L.rows {
            assert( approx( L2 * res, b2 ) );
        } else {
            L3 := cmat64_submatrix(L, L.cols, L.rows, 0, L.cols);
            b3 := cvec64_subvector(b, L.cols, L.rows);
            assert( approx(L3 * res, b3) );
        }
    }

    return res;
}

// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 3.4.5.2 Right Trapezoidal, page 160
// There are degrees of freedom left. These are set to 0 here.
slrtz :: solve_linear_right_trapezoidal;
solve_linear_right_trapezoidal :: (R: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(R.rows == b.dim);
        assert(is_right_trapezoidal(R));
    }

    R2 := cmat64_submatrix(R, 0, R.rows, 0, R.rows);

    res := solve_linear_right_triangular(R2, b);

    rvec := cvec64_zeros(R.cols);
    for 0..res.dim-1 {
        rvec.data[it] = res[it];
    }

    return rvec;
}


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7.1 Gaussian Factorization, page 199
gaussian_factorization :: (M: CMatrix64) -> CMatrix64 {
    A := copy(M);
    m := M.rows;
    n := M.cols;
    for k: 0..min(m,n)-1 {
        assert(A[k][k] != CZERO, "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        for i: k+1..m-1 {
            A[i][k] = A[i][k]/A[k][k];
        }
        for j: k+1..n-1 {
            for i: k+1..m-1 {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }
    return A;
}


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7.1 Gaussian Factorization, page 199
decompose_lr :: (A: CMatrix64) -> L: CMatrix64, R: CMatrix64 {
    M := gaussian_factorization(A);

    L := cmat64_id(M.rows, M.cols);
    R := cmat64_zeros(M.rows, M.cols);

    L[0][0] = c64(1.0, 0.0);
    for i: 0..M.rows-1 {
        for j: 0..M.cols-1 {
            if j >= i {
                R[i][j] = M[i][j];
            } else {
                L[i][j] = M[i][j];
            }
        }
    }

    return L, R;
}




// // Scientific Computing, Vol I: Linear and nonlinear equations
// // Texts in computational science and engineering 18
// // Springer
// // Algorithm 3.7.2 Solve Linear System via Gaussian Factorization, page 199
slgf :: solve_linear_gaussian_factorization;
solve_linear_gaussian_factorization :: (A: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(A.rows == b.dim);
    }

    L, R := decompose_lr(A);

    y := solve_linear_left_trapezoidal(L, b);
    x := solve_linear_right_trapezoidal(R, y);

    return x;
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

    {
        M := cmat64(3, 2, 
            R, CZERO, 
            -R, R, 
            -R, -R
        );
        b := cvec64(R, 2.0*R, -4.0*R);

        res := solve_linear_left_trapezoidal(M, b);
        print("sllt: %\n", res);
        assert( approx(M*res, b) );
    }

    {
        M := cmat64(2, 3, 
            R, R, CZERO, 
            CZERO, R, CZERO 
        );
        b := cvec64(4*R, 6*R);

        res := solve_linear_right_trapezoidal(M, b);
        print("slrt: %\n", res);
        assert( approx(M*res, b) );
    }
    {
        M := cmat64(2,2,
            2*R, 3*R,
            -2*R, R
        );
        res := cmat64(2,2,
            2*R, 3*R,
            -1*R, 4*R
        );
        A := gaussian_factorization(M);
        print("gaussian factorization:\n    %\n -> %\n", M, A);
        assert( A == res , "\n% != %!", A, res);
    }

    {
        M := cmat64(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := cvec64(4*R, 6*I);

        res := solve_linear_gaussian_factorization(M, b);
        print("solve linear system via Gaussian factorization: %\n", res);
        assert(res == cvec64(2.5-4.5*I, 4.5-1.5*I), "\nCould not solve linear system via Gaussain factorization!");
    }
}