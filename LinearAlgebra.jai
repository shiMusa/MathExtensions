






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.1 2x2 Linear System Solve, page 153
// solves the system A x = b with 2x2 matrix A
solve_2x2 :: (A: CMatrix64, b: CVector64) -> CVector64 {
    assert(A.rows == 2 && A.cols == 2);

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
// solve the system Lx = b where L is a left-triangluar matrix
solve_left_triangular :: (L: CMatrix64, b: CVector64) -> CVector64 {
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















test_linear_algebra :: () {
    print("\n\nLinear Algebra -------------------------------\n\n");

    {
        A := cmat64(2, 2, R, R, 2.0 - I, I);
        b := cvec64(I, c64(2.0, 4.0));
        res := solve_2x2(A, b);
        print("%\n", res);
        assert(res == cvec64(
            -0.25 + 1.75*I,
            0.25 - 0.75*I
        ));
    }
}