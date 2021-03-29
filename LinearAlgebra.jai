






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.1 2x2 Linear System Solve, page 153
// solves the system A x = b with 2x2 matrix A
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
solve_linear_left_triangular :: (L: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(cmat64_is_quadratic(L), "L is not quadratic! % rows, % columns", L.rows, L.cols);
        assert(cmat64_is_diagonal_unit(L), "L is not unit: Diagonal entries are not 1");
        assert(cmat64_is_left_triangular(L), "L is not left-triangular");
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
solve_linear_right_triangular :: (R: CMatrix64, b: CVector64) -> CVector64 {
    #if CHECKS {
        assert(cmat64_is_quadratic(R), "R is not quadratic! % rows, % columns", R.rows, R.cols);
        assert(cmat64_is_right_triangular(R), "R is not right-triangular");
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
// solve_linear_right_quasi_triangular :: (R: CMatrix64, b: CVector64) -> CVector64 {
//     #if CHECKS {
//         assert(cmat64_is_quadratic(R), "R is not quadratic! % rows, % columns", R.rows, R.cols);
//         assert(cmat64_is_right_quasi_triangular(R), "R is not right quasi-triangular!");
//         // for i: 0..R.rows-1 {
//         //     for j: 0..R.cols-1 {
//         //         if i>j+1 then assert(R[i][j] == CZERO, "R is not right quasi-triangular: R[%][%] = % != 0", i,j,R[i][j]);
//         //     }
//         //     if i>=2 {
//         //         // diagonal blocks
//         //         if R[i-1][i] != CZERO {
//         //             assert(R[i-2][i-1] == CZERO && R[i][i+1] == CZERO, "R is not right quasi-triangular");
//         //         }
//         //         // 1x1 diagonal blocks invertible
//         //         if R[i-1][i] == CZERO && R[i+1][i] == CZERO {
//         //             assert(R[i][i] != CZERO, "R is not right quasi-triangular");)
//         //         }
//         //         // 2x2 diagonal blocks are invertible
//         //         if R[i][i-1] != CZERO {
//         //             assert(R[i-1][i-1]*R[i][i] - R[i-1][i]*R[i][i-1] != CZERO, "R is not right quasi-triangular");
//         //         }
//         //     }
//         // }
//     }


// }










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
}