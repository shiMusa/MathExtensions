






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.4.1 2x2 Linear System Solve, page 153
// solves the system A x = b with 2x2 matrix A
sl2 :: solve_linear_2x2;
solve_linear_2x2 :: (A: CMatrix64, b: CVector64) -> CVector64 #must {
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
solve_linear_left_triangular :: (L: CMatrix64, b: CVector64) -> CVector64 #must {
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
solve_linear_right_triangular :: (R: CMatrix64, b: CVector64) -> CVector64 #must {
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
solve_linear_right_quasi_triangular :: (R: CMatrix64, b: CVector64) -> CVector64 #must {
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
solve_linear_orthogonal_projection :: (Q: CMatrix64, b: CVector64) -> CVector64 #must {
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
solve_linear_successive_orthogonal_projection :: (Q: CMatrix64, b: CVector64) -> CVector64 #must {
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
solve_linear_left_trapezoidal :: (L: CMatrix64, b: CVector64) -> CVector64 #must {
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

solve_linear_left_trapezoidal :: (L: CMatrix64, M: CMatrix64) -> CMatrix64 #must {
    #if CHECKS {
        assert(L.rows == M.rows);
    }

    res := cmat64(L.cols, M.cols);

    for j: 0..M.cols-1 {
        b := column(M, j);
        tmp := solve_linear_left_trapezoidal(L, b);
        for tmp.data {
            res[it_index][j] = it;
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
solve_linear_right_trapezoidal :: (R: CMatrix64, b: CVector64) -> CVector64 #must {
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


solve_linear_right_trapezoidal :: (R: CMatrix64, M: CMatrix64) -> CMatrix64 #must {
    #if CHECKS {
        assert(R.rows == M.rows);
    }
    
    res := cmat64(R.cols, M.cols);

    for j: 0..M.cols-1 {
        b := column(M, j);
        tmp := solve_linear_right_trapezoidal(R, b);
        for tmp.data {
            res[it_index][j] = it;
        }
    }

    return res;
}






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7.1 Gaussian Factorization - No Pivoting, page 199
decompose_LR :: (M: CMatrix64) -> L: CMatrix64 #must, R: CMatrix64 #must {

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







// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7.1 Gaussian Factorization - No Pivoting, page 199
gaussian_factorization_no_pivot :: (M: CMatrix64) -> CMatrix64 #must {
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




// // Scientific Computing, Vol I: Linear and nonlinear equations
// // Texts in computational science and engineering 18
// // Springer
// // Algorithm 3.7.2 Solve Linear System via Gaussian Factorization - No Pivoting, page 199
slgf_np :: solve_linear_gaussian_factorization_no_pivot;
solve_LRnp :: solve_linear_gaussian_factorization_no_pivot;
solve_LUnp :: solve_linear_gaussian_factorization_no_pivot;
solve_linear_gaussian_factorization_no_pivot :: (A: CMatrix64, b: CVector64) -> CVector64 #must {
    #if CHECKS {
        assert(A.rows == b.dim);
    }

    M := gaussian_factorization_no_pivot(A);
    L, R := decompose_LR(M);

    y := solve_linear_left_trapezoidal(L, b);
    x := solve_linear_right_trapezoidal(R, y);

    return x;
}








// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7.3 Gaussian Factorization with Full Pivoting, page 204
gaussian_factorization_full_pivot :: (M: CMatrix64) -> 
    CMatrix64 #must, rows: []int #must, cols: []int #must, 
    row_permutations: int, col_permutations: int {

    A := copy(M);
    m := A.rows;
    n := A.cols;

    // rows : [..]int;
    // cols : [..]int;
    // array_resize(*rows, A.rows, false);
    // array_resize(*cols, A.cols, false);
    rows := NewArray(A.rows, int, false);
    cols := NewArray(A.cols, int, false);
    row_permutations := 0;
    col_permutations := 0;

    for 0..A.rows-1 {
        rows[it] = it;
    }
    for 0..A.cols-1 {
        cols[it] = it;
    }


    for k: 0..min(m,n)-1 {
        ii := k;
        jj := k;
        rho := abs(A[k][k]);

        for j: k..n-1 {
            for i: k..m-1 {
                tmp := abs(A[i][j]);
                if tmp > rho {
                    ii = i;
                    jj = j;
                    rho = tmp;
                }
            }
        }

        assert(rho != 0.0, "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        // instead of swapping, use the permutation arrays rows and cols to access the 
        // entries of A. Should be faster?
        if ii > k {
            swap_rows(A, k, ii);
            rows[k], rows[ii] = swap(rows[k], rows[ii]);
            row_permutations += 1;
        }
        if jj > k {
            swap_columns(A, k, jj);
            cols[k], cols[jj] = swap(cols[k], cols[jj]);
            col_permutations += 1;
        }

        for i: k+1..m-1 {
            A[i][k] = A[i][k]/A[k][k];
        }
        for j: k+1..n-1 {
            for i: k+1..m-1 {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }

    return A, rows, cols, row_permutations, col_permutations;
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7.4 Solve Linear System via Gaussian Factorization with Full Pivoting, page 206
slgf_fp :: solve_linear_gaussian_factorization_full_pivot;
solve_LRfp :: solve_linear_gaussian_factorization_full_pivot;
solve_LUfp :: solve_linear_gaussian_factorization_full_pivot;
solve_linear_gaussian_factorization_full_pivot :: (A: CMatrix64, b: CVector64) -> CVector64 #must {
    #if CHECKS {
        assert(A.rows == b.dim);
    }

    M, rows, cols := gaussian_factorization_full_pivot(A);
    L, R := decompose_LR(M);

    y := solve_linear_left_trapezoidal(L, permute(b, ..rows));
    w := solve_linear_right_trapezoidal(R, y);
    x := permute(w, ..cols);

    return x;
}










// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7. Gaussian Factorization with Partial Pivoting, page 210
gaussian_factorization_row_pivot :: (M: CMatrix64) -> CMatrix64 #must, rows: []int, permutations: int {
    A := copy(M);
    m := A.rows;
    n := A.cols;

    // rows : [..]int;
    // array_resize(*rows, A.rows, false);
    rows := NewArray(A.rows, int, false);
    permutations := 0;

    for 0..A.rows-1 {
        rows[it] = it;
    }


    for k: 0..min(m,n)-1 {
        ii := k;
        rho := abs(A[k][k]);

        for i: k..m-1 {
            tmp := abs(A[i][k]);
            if tmp > rho {
                ii = i;
                rho = tmp;
            }
        }

        assert(rho != 0.0, "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        // instead of swapping, use the permutation arrays rows and cols to access the 
        // entries of A. Should be faster?
        if ii > k {
            swap_rows(A, k, ii);
            rows[k], rows[ii] = swap(rows[k], rows[ii]);
            permutations += 1;
        }

        for i: k+1..m-1 {
            A[i][k] = A[i][k]/A[k][k];
        }
        for j: k+1..n-1 {
            for i: k+1..m-1 {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }

    return A, rows, permutations;
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.7.7 Solve Linear System via Gaussian Factorization with Partial Pivoting, page 215
slgf_rp :: solve_linear_gaussian_factorization_row_pivot;
solve_LRrp :: solve_linear_gaussian_factorization_row_pivot;
solve_LUrp :: solve_linear_gaussian_factorization_row_pivot;
solve_linear_gaussian_factorization_row_pivot :: (A: CMatrix64, b: CVector64) -> CVector64 #must {
    #if CHECKS {
        assert(A.rows == b.dim);
    }

    M, rows := gaussian_factorization_row_pivot(A);
    L, R := decompose_LR(M);

    y := solve_linear_left_trapezoidal(L, permute(b, ..rows));
    w := solve_linear_right_trapezoidal(R, y);

    return w;
}






// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 3.10.1 Matrix Inverse, page 261
inverse :: (A: CMatrix64) -> CMatrix64 {
    #if CHECKS {
        assert(is_quadratic(A));
    }
    M, rows := gaussian_factorization_row_pivot(A);
    L, R := decompose_LR(M);

    Y := solve_linear_left_trapezoidal(L, cmat64_id(A.rows, A.cols));
    X := solve_linear_right_trapezoidal(R, Y);
    return permute_columns(X, ..rows);
}





// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 3.11 Determinants, page 263
//
// it is not recommended to calculate the dterminant of large matrices
det :: determinant;
determinant :: (A: CMatrix64) -> Complex64 {
    #if CHECKS {
        assert(is_quadratic(A));
    }

    if A.rows == 2 {
        return A[0][0]*A[1][1] - A[0][1]*A[1][0];
    } // implement optimized formula for 3x3 matrix

    M, rows, permutations := gaussian_factorization_row_pivot(A);
    L, R := decompose_LR(M);

    p := c64(1.0);
    for 0..R.rows-1 {
        p *= R[it][it];
    }

    // +1 for even, -1 for odd num of permutations
    s := 2.0 * cast(float64)(permutations % 2 == 0) - 1.0;

    return s * p;
}








test_linear_algebra :: () {
    print("\n\nLinear Algebra -------------------------------\n\n");

    {
        A := cmat64(2, 2, R, R, 2.0 - I, I);
        b := cvec64(I, c64(2.0, 4.0));
        res := solve_linear_2x2(A, b);
        print("%\n", str(res));
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
        print("% = % ?\n", str(L * res), str(b, ""));
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
        print("res = %\n", str(res));

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
        print("sllt: %\n", str(res));
        assert( approx(M*res, b) );
    }

    {
        M := cmat64(2, 3, 
            R, R, CZERO, 
            CZERO, R, CZERO 
        );
        b := cvec64(4*R, 6*R);

        res := solve_linear_right_trapezoidal(M, b);
        print("slrt: %\n", str(res));
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
        A := gaussian_factorization_no_pivot(M);
        print("gaussian factorization: % -> %\n", str(M), str(A));
        assert( A == res , "%\n != %\n!", str(A), str(res));
    }

    {
        M := cmat64(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := cvec64(4*R, 6*I);

        res := solve_LRnp(M, b);
        print("solve linear system via Gaussian factorization (no pivoting): %\n", str(res));
        assert(res == cvec64(2.5-4.5*I, 4.5-1.5*I), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        M := cmat64(3,3,
            CZERO, 4*R, -3*R,
            2*R, CZERO, 4*R,
            3*R, -5*R, 3*R
        );

        A, rows, cols := gaussian_factorization_full_pivot(M);

        print("A = %\n", str(A));
    }

    {
        M := cmat64(3,3,
            CZERO, 4*R, -3*R,
            2*R, CZERO, 4*R,
            3*R, -5*R, 3*R
        );

        A, rows, cols := gaussian_factorization_full_pivot(M);
        L, R := decompose_LR(A);

        print("L = %\nR = %\n", str(L), str(R));
    }

    {
        M := cmat64(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := cvec64(4*R, 6*I);

        res := solve_LRfp(M, b);
        print("solve linear system via Gaussian factorization (full pivoting): %\n", str(res));
        assert(approx(res, cvec64(2.5-4.5*I, 4.5-1.5*I)), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        M := cmat64(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := cvec64(4*R, 6*I);

        res := solve_LRrp(M, b);
        print("solve linear system via Gaussian factorization (row pivoting): %\n", str(res));
        assert(approx(res, cvec64(2.5-4.5*I, 4.5-1.5*I)), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        L := cmat64(3,3,
            R, CZERO, CZERO,
            I, R, CZERO,
            2*R, 3+I, R
        );
        one := cmat64_id(3,3);
        res := solve_linear_left_trapezoidal(L,one);
        print("K = %\n", str(res));
        assert(approx(res, cmat64(3,3,
            R, CZERO, CZERO,
            -I, R, CZERO, 
            -3+3*I, -3-I, R
        )));
    }

    {
        M := cmat64(3,3,
            2*R, I, -2*R,
            CZERO, -I, CZERO,
            CZERO, CZERO, 1+I
        );
        one := cmat64_id(3,3);
        res := solve_linear_right_trapezoidal(M,one);
        print("K = %\n", str(res));
        assert(approx(res, cmat64(3,3,
            0.5*R, 0.5*R, 0.5-0.5*I,
            CZERO, I, CZERO,
            CZERO, CZERO, 0.5-0.5*I
        )));
    }

    {
        M := cmat64(3, 3,
            2*R, I, -2*R,
            CZERO, -I, 2*I,
            R, CZERO, R+I
        );
        iM := inverse(M);
        print("iM = %\n", str(iM));
        assert(approx(iM, cmat64(3,3,
            c64(0.25, 0.25), c64(0.25, 0.25), c64(0.5, -0.5),
            -0.5*R, -0.5*R+I, R,
            -0.25*R, -0.25*R, 0.5*R
        )));
    }

    {
        M := cmat64(3, 3,
            2*R, I, -2*R,
            CZERO, -I, 2*I,
            R, CZERO, R+I
        );
        d := determinant(M);
        print("det(M) = %\n", d);
        assert(d == -4*I);
    }
}