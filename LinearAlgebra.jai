






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.1 2x2 Linear System Solve, page 153
// solves the system A x = b with 2x2 matrix A
solve_linear_2x2 :: (A: CMatrix64, b: CVector64, $check: bool = false) -> CVector64 {
    #if check {
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
solve_linear_left_triangular :: (L: CMatrix64, b: CVector64, $check: bool = false) -> CVector64 {
    #if check {
        assert(L.cols == L.rows, "L is not quadratic! % rows, % columns", L.rows, L.cols);
        for i: 0..L.cols-1 {
            assert(L[i][i] == R, "L is not unit: Diagonal entries are not 1: L[%][%] = %", i, i, L[i][i]);
        }
        for r: 0..L.rows-2 {
            for c: r+1..L.cols-1 {
                assert(L[r][c] == CZERO, "L is not left-triangular: L[%][%] = %", r, c, L[r][c]);
            }
        }
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
solve_linear_right_triangular :: (R: CMatrix64, b: CVector64, $check: bool = false) -> CVector64 {
    #if check {
        assert(R.cols == R.rows, "R is not quadratic! % rows, % columns", R.rows, R.cols);
        for r: 0..R.rows-2 {
            for c: 0..r-1 {
                assert(R[r][c] == CZERO, "R is not right-triangular: R[%][%] = %", r, c, R[r][c]);
            }
        }
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















test_linear_algebra :: () {
    print("\n\nLinear Algebra -------------------------------\n\n");

    {
        A := cmat64(2, 2, R, R, 2.0 - I, I);
        b := cvec64(I, c64(2.0, 4.0));
        res := solve_linear_2x2(A, b, true);
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
        res := solve_linear_left_triangular(L, b, true);
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
        res := solve_linear_right_triangular(M, b, true);

        test := norm_2(M * res - b);
        print("%\n", test);
        assert(test < 1.0e-10, "res = %, M * res = % =? b = %", res, M*res, b);
    }
}