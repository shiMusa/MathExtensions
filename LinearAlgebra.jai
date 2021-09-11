






// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.1 2x2 Linear System Solve, page 153
// ? solves the system A x = b with 2x2 matrix A
sl2 :: solve_linear_2x2;
solve_linear_2x2 :: (A: Matrix($T, $N, N), b: Vector($B, N)) -> Vector(UpCast(T,B).T, N) {
    #run assert(N == 2, "dimension mismatch: 2 != N = %\n", N);
    U :: UpCast(T,B).T;

    sigma := scl(max(abs(A[0][0]), abs(A[0][1]), abs(A[1][0]), abs(A[1][1])));

    mu := A / sigma;
    det := scl(ncast(U, mu[0][0] * mu[1][1] - mu[0][1] * mu[1][0]));
    g := b / sigma;
    return mat(2, 2, mu[1][1], -mu[0][1], -mu[1][0], mu[0][0]) * g / det;
}







// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.2 Forward Linear Recurrence, page 155
// ? solve the system Lx = b where L is a unit left-triangluar matrix
// ? that means, all diagonals are 1 and all L_ij = 0 for j > i.
sllta :: solve_linear_left_triangular;
solve_linear_left_triangular :: (L: Matrix($T, $N, N), b: Vector($V, N)) -> Vector(UpCast(T,V).T, N) {
    #if CHECKS {
        assert(is_diagonal_unit(L), "L is not unit: Diagonal entries are not 1");
        assert(is_left_triangular(L), "L is not left-triangular");
    }
    U :: UpCast(T,V).T;

    res : Vector(U, N);
    for i: 0..N-1 {
        sum := zero(U);
        for j: 0..i-1 {
            sum += L[i][j] * res.data[j];
        }
        res.data[i] = b[i] - sum;
    }
    return res;
}




// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.4 Backward Linear Recurrence, page 158
// ? solve the system Rx = b where R is a right-triangluar matrix
// ? that means all R_ij = 0 for j < i.
slrta :: solve_linear_right_triangular;
solve_linear_right_triangular :: (R: Matrix($T, $N, N), b: Vector($V, N)) -> Vector(UpCast(T,V).T, N) {
    #if CHECKS {
        assert(is_right_triangular(R), "R is not right-triangular");
    }
    U :: UpCast(T,V).T;
    
    res : Vector(U, N);
    i : int;
    for 0..N-1 {
        i = N-1 - it;
        sum := zero(U);
        for j: i+1..N-1 {
            sum += R[i][j]*res.data[j];
        }
        res.data[i] = (b[i] - sum)/R[i][i];
    }
    return res;
}







// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.6 Solve Right Quasi-Triangular System, page 164
// ? solve the system Rx=b where R is a right(upper) quasi-triangular matrix.
slrqta :: solve_linear_right_quasi_triangular;
solve_linear_right_quasi_triangular :: (R: Matrix($T, $N, N), b: Vector($V, N)) -> Vector(UpCast(T,V).T, N) {
    #if CHECKS {
        assert(is_right_quasi_triangular(R), "R is not right quasi-triangular!");
    }
    U :: UpCast(T,V).T;
    res : Vector(UpCast(T,V).T, N) = ---;

    i := N-1; // start with last row and wrap from bottom to up
    while i >= 0 {
        // check if 2x2 block by looking at left column
        // the term `|| i == 0` was missing
        if (i > 0 && R[i][i-1] == zero(U)) || i == 0 {
            sum := zero(U);
            for j: i+1..N-1 {
                sum += R[i][j] * res.data[j];
            }
            res.data[i] = (b[i] - sum)/R[i][i];
            i -= 1;
        } else if i > 0 { // so now we have a 2x2 block where R[i][i] is the bottom-right element
            A := mat(2, 2, 
                R[i-1][i-1], R[i-1][i],
                R[i][i-1], R[i][i]
            );
            sum1 := zero(U);
            sum2 := zero(U);
            for j: i+1..N-1 {
                sum1 += R[i-1][j] * res.data[j];
                sum2 += R[i][j] * res.data[j];
            }
            c := vec(2, U,
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







// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.8 Orthogonal Projection, page 167
solve_linear_unitary :: solve_linear_orthogonal_projection;
slop :: solve_linear_orthogonal_projection;
solve_linear_orthogonal_projection :: (Q: Matrix($T, $N, $M), b: Vector($V, N)) -> Vector(UpCast(T,V).T, M) {
    #if CHECKS {
        assert(is_unitary(Q));
    }

    return dagger(Q) * b;
}






// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.9 Successive Orthogonal Projection, page 168
// ? not sure if this algo is useful. seems to be slower than the previous one.
slsop :: solve_linear_successive_orthogonal_projection;
solve_linear_successive_orthogonal_projection :: (Q: Matrix($T, $N, $M), b: Vector($V, M)) -> Vector(UpCast(T,V).T, M) {
    #if CHECKS {
        assert(is_unitary(Q));
    }
    bb := copy(b);

    res : Vector(UpCast(T,V).T, M) = ---;
    for i: 0..M-1 {
        q := row(Q, i);
        res.data[i] = conjugate(q) * bb; // not sure if this is performant...
        bb = bb - q * res.data[i];
    }
    return res;
}





// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 3.4.5.1 Left Trapezoidal, page 159
slltz :: solve_linear_left_trapezoidal;
solve_linear_left_trapezoidal :: (L: Matrix($T, $N, $M), b: Vector($V, N)) -> Vector(UpCast(T,V).T, M) {
    #if CHECKS {
        assert(is_left_trapezoidal(L)); // more rows than columns N >= M
    }

    L2 := mat_submatrix(L, 0, M, 0, M);
    b2 := vec_subvector(b, 0, M);

    res := solve_linear_left_triangular(L2, b2);

    #if CHECKS {
        #if N == M {
            assert( approx( L2 * res, b2 ) );
        } else {
            L3 := mat_submatrix(L, M, N, 0, M);
            b3 := vec_subvector(b, M, N);
            assert( approx(L3 * res, b3) );
        }
    }

    return res;
}


solve_linear_left_trapezoidal :: (l: Matrix($T, $N, $M), m: Matrix($U, N, $K)) -> Matrix(UpCast(T,U).T, M, K) {
    res : Matrix(UpCast(T,U).T, N, K) = ---;

    for j: 0..K-1 {
        b := column(m, j);
        tmp := solve_linear_left_trapezoidal(l, b);
        for tmp.data {
            res[it_index][j] = it;
        }
    }

    return res;
}




// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 3.4.5.2 Right Trapezoidal, page 160
// ? There are degrees of freedom left. These are set to 0 here.
slrtz :: solve_linear_right_trapezoidal;
solve_linear_right_trapezoidal :: (R: Matrix($T, $N, $M), b: Vector($V, N)) -> Vector(UpCast(T,V).T, M) {
    #if CHECKS {
        assert(is_right_trapezoidal(R), "Matrix is not right trapezoidal: M(%) >= N(%)?\n", M, N);
    }

    R2 := mat_submatrix(R, 0, N, 0, N);

    tmp := solve_linear_right_triangular(R2, b);

    res : Vector(UpCast(T,V).T, M);
    for 0..tmp.N-1 {
        res.data[it] = tmp[it];
    }
    return res;
}


solve_linear_right_trapezoidal :: (r: Matrix($T, $N, $M), m: Matrix($V, N, $K)) -> Matrix(UpCast(T, V).T, M, K) {
    res : Matrix(UpCast(T, V).T, M, K) = ---;

    for j: 0..K-1 {
        b := column(m, j);
        tmp := solve_linear_right_trapezoidal(r, b);
        for tmp.data {
            res[it_index][j] = it;
        }
    }

    return res;
}






// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.1 Gaussian Factorization - No Pivoting, page 199
decompose_LR :: (m: Matrix($T, $N, $M)) -> L: Matrix(T, N, M) #must, R: Matrix(T, N, M) #must {

    L := mat_id(T, N, M);
    R : Matrix(T, N, M) = ---;
    // L[0][0] = one(T); // this should be unnecessary

    for i: 0..N-1 {
        for j: 0..M-1 {
            if j >= i {
                R[i][j] = m[i][j];
            } else {
                L[i][j] = m[i][j];
            }
        }
    }

    return L, R;
}







// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.1 Gaussian Factorization - No Pivoting, page 199
gaussian_factorization_no_pivot :: (m: Matrix($T, $N, $M)) -> Matrix(T, N, M) #must {
    A := copy(m);

    for k: 0..#run min(N,M)-1 {
        assert(abs(A[k][k]) != zero(T), "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        for i: k+1..N-1 {
            A[i][k] = A[i][k]/A[k][k];
        }
        for j: k+1..M-1 {
            for i: k+1..N-1 {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }
    return A;
}




// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.2 Solve Linear System via Gaussian Factorization - No Pivoting, page 199
slgf_np :: solve_linear_gaussian_factorization_no_pivot;
solve_LRnp :: solve_linear_gaussian_factorization_no_pivot;
solve_LUnp :: solve_linear_gaussian_factorization_no_pivot;
solve_linear_gaussian_factorization_no_pivot :: (A: Matrix($T, $N, $M), b: Vector($V, N)) -> Vector(UpCast(T,V).T, M) {
    M := gaussian_factorization_no_pivot(A);
    L, R := decompose_LR(M);
    y := solve_linear_left_trapezoidal(L, b);
    return solve_linear_right_trapezoidal(R, y);
}








// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.3 Gaussian Factorization with Full Pivoting, page 204
gaussian_factorization_full_pivot :: (m: Matrix($T, $N, $M))
 -> Matrix(T, N, M) #must, 
    rows: [N]int #must, 
    cols: [M]int #must, 
    row_permutations: int, 
    col_permutations: int 
{
    A := copy(m);
    rows : [N]int;
    cols : [M]int;
    row_permutations := 0;
    col_permutations := 0;

    for 0..N-1 {
        rows[it] = it;
    }
    for 0..M-1 {
        cols[it] = it;
    }


    for k: 0..min(N,M)-1 {
        ii := k;
        jj := k;
        rho := abs(A[k][k]);

        for j: k..M-1 {
            for i: k..N-1 {
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
            swap_rows(*A, k, ii);
            rows[k], rows[ii] = swap(rows[k], rows[ii]);
            row_permutations += 1;
        }
        if jj > k {
            swap_columns(*A, k, jj);
            cols[k], cols[jj] = swap(cols[k], cols[jj]);
            col_permutations += 1;
        }

        for i: k+1..N-1 {
            A[i][k] = A[i][k]/A[k][k];
        }
        for j: k+1..M-1 {
            for i: k+1..N-1 {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }

    return A, rows, cols, row_permutations, col_permutations;
}








// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.4 Solve Linear System via Gaussian Factorization with Full Pivoting, page 206
slgf_fp :: solve_linear_gaussian_factorization_full_pivot;
solve_LRfp :: solve_linear_gaussian_factorization_full_pivot;
solve_LUfp :: solve_linear_gaussian_factorization_full_pivot;
solve_linear_gaussian_factorization_full_pivot :: (A: Matrix($T, $N, $M), b: Vector($V, N)) -> Vector(UpCast(T,V).T, N) {
    M, rows, cols := gaussian_factorization_full_pivot(A);
    L, R := decompose_LR(M);

    y := solve_linear_left_trapezoidal(L, permute(b, ..rows));
    w := solve_linear_right_trapezoidal(R, y);
    return permute(w, ..cols);
}










// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7. Gaussian Factorization with Partial Pivoting, page 210
gaussian_factorization_row_pivot :: (m: Matrix($T, $N, $M)) -> Matrix(T, N, M) #must, rows: [N]int, permutations: int {
    A := copy(m);

    rows : [N]int = ---;
    permutations := 0;

    for 0..N-1 {
        rows[it] = it;
    }

    for k: 0..min(N,M)-1 {
        ii := k;
        rho := abs(A[k][k]);

        for i: k..N-1 {
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
            swap_rows(*A, k, ii);
            rows[k], rows[ii] = swap(rows[k], rows[ii]);
            permutations += 1;
        }

        for i: k+1..N-1 {
            A[i][k] = A[i][k]/A[k][k];
        }
        for j: k+1..M-1 {
            for i: k+1..N-1 {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }

    return A, rows, permutations;
}



// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.7 Solve Linear System via Gaussian Factorization with Partial Pivoting, page 215
slgf_rp :: solve_linear_gaussian_factorization_row_pivot;
solve_LRrp :: solve_linear_gaussian_factorization_row_pivot;
solve_LUrp :: solve_linear_gaussian_factorization_row_pivot;
solve_linear_gaussian_factorization_row_pivot :: (A: Matrix($T, $N, $M), b: Vector($V, N)) -> Vector(UpCast(T,V).T, N) {
    M, rows := gaussian_factorization_row_pivot(A);
    L, R := decompose_LR(M);

    y := solve_linear_left_trapezoidal(L, permute(b, ..rows));
    return solve_linear_right_trapezoidal(R, y);
}






// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.10.1 Matrix Inverse, page 261
inverse :: (A: Matrix($T, $N, N)) -> Matrix(T, N, N) {
    M, rows := gaussian_factorization_row_pivot(A);
    L, R := decompose_LR(M);

    Y := solve_linear_left_trapezoidal(L, mat_id(T, N, N));
    X := solve_linear_right_trapezoidal(R, Y);

    return permute_columns(X, ..rows);
}





// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 3.11 Determinants, page 263
// ?
// ? it is not recommended to calculate the dterminant of large matrices
det :: determinant;
determinant :: (A: Matrix($T, $N, N)) -> T {

    // implement optimized formula for 3x3 matrix
    #if N == 2 {
        return A[0][0]*A[1][1] - A[0][1]*A[1][0];
    } else {    
        M, rows, permutations := gaussian_factorization_row_pivot(A);
        L, R := decompose_LR(M);
    
        p := one(T);
        for 0..R.N-1 {
            p *= R[it][it];
        }

        // +1 for even, -1 for odd num of permutations
        s : float64 = 2.0 * cast(float64)(permutations % 2 == 0) - 1.0;
    
        return s * p;
    }
}































//* //////////////////////////////////////////////////////////////////////////////////////
/*
      ######## ########  ######  ########  ######
           ##    ##       ##    ##    ##    ##    ##
           ##    ##       ##          ##    ##
           ##    ######    ######     ##     ######
           ##    ##             ##    ##          ##
           ##    ##       ##    ##    ##    ##    ##
           ##    ########  ######     ##     ######
*/
//* //////////////////////////////////////////////////////////////////////////////////////



test_linear_algebra :: () {
    print("\n\nLinear Algebra -------------------------------\n\n");

    {
        A := mat(2, 2, R, R, 2.0 - I, I);
        b := vec(2, I, c64(2.0, 4.0));
        res := solve_linear_2x2(A, b);
        print("%\n", str(res));
        assert(res == vec(2, 
            -0.25 + 1.75*I,
            0.25 - 0.75*I
        ));
    }

    {
        L := mat(3, 3,
            R, CZERO, CZERO,
            3.0*R, R, CZERO,
            I, R+I, R
        );
        b := vec(3, 3.0*R, -I, R-2.0*I);
        res := solve_linear_left_triangular(L, b);
        print("% = % ?\n", str(L * res), str(b, ""));
        assert(L * res == b);
    }

    {
        M := mat(3, 3,
            2.0*R, I, R-2.0*I,
            CZERO, R, 3.0*I,
            CZERO, CZERO, -3.0*R
        );
        b := vec(3, 3.0*R, -I, R-2.0*I);
        res := solve_linear_right_triangular(M, b);

        test := norm_2(M * res - b);
        print("%\n", test); 
        assert(test < 1.0e-10, "res = %, M * res = % =? b = %", res, M*res, b);
    }

    {
        m := mat(4, 4,
            1, 3, 4, 5,
            0, 2, 3, 4,
            0, 1, 2, 3,
            0, 0, 0, 2
        );
        cm := mat_cast(Complex64, m);
        
        b := vec(4, float64, 2,3,4,5);
        cb := vec_cast(Complex64, b);

        res := solve_linear_right_quasi_triangular(cm,cb);
        print("res = %\n", str(res));

        test := norm_2(cm * res - cb);
        assert(test < 1.0e-10);
    }

    {
        M := mat(3, 2, 
            R, CZERO, 
            -R, R, 
            -R, -R
        );
        b := vec(3,   R, 2.0*R, -4.0*R);

        res := solve_linear_left_trapezoidal(M, b);
        print("sllt: %\n", str(res));
        assert( approx(M*res, b) );
    }

    {
        M := mat(2, 3, 
            R, R, CZERO, 
            CZERO, R, CZERO 
        );
        b := vec(2,   4*R, 6*R);

        res := solve_linear_right_trapezoidal(M, b);
        print("slrt: %\n", str(res));
        assert( approx(M*res, b) );
    }
    {
        M := mat(2,2,
            2*R, 3*R,
            -2*R, R
        );
        res := mat(2,2,
            2*R, 3*R,
            -1*R, 4*R
        );
        A := gaussian_factorization_no_pivot(M);
        print("gaussian factorization: % -> %\n", str(M), str(A));
        assert( A == res , "%\n != %\n!", str(A), str(res));
    }

    {
        M := mat(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := vec(2, Complex64, 4*R, 6*I);

        res := solve_LRnp(M, b);
        print("solve linear system via Gaussian factorization (no pivoting): %\n", str(res));
        assert(res == vec(2, Complex64, 2.5-4.5*I, 4.5-1.5*I), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        M := mat(3,3,
            CZERO, 4*R, -3*R,
            2*R, CZERO, 4*R,
            3*R, -5*R, 3*R
        );

        A, rows, cols := gaussian_factorization_full_pivot(M);

        print("A = %\n", str(A));
    }

    {
        M := mat(3,3,
            CZERO, 4*R, -3*R,
            2*R, CZERO, 4*R,
            3*R, -5*R, 3*R
        );

        A, rows, cols := gaussian_factorization_full_pivot(M);
        L, R := decompose_LR(A);

        print("L = %\nR = %\n", str(L), str(R));
    }

    {
        M := mat(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := vec(2,   4*R, 6*I);

        res := solve_LRfp(M, b);
        print("solve linear system via Gaussian factorization (full pivoting): %\n", str(res));
        assert(approx(res, vec(2,   2.5-4.5*I, 4.5-1.5*I)), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        M := mat(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := vec(2,   4*R, 6*I);

        res := solve_LRrp(M, b);
        print("solve linear system via Gaussian factorization (row pivoting): %\n", str(res));
        assert(approx(res, vec(2,   2.5-4.5*I, 4.5-1.5*I)), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        L := mat(3,3,
            R, CZERO, CZERO,
            I, R, CZERO,
            2*R, 3+I, R
        );
        one := mat_id(Complex64, 3, 3);
        res := solve_linear_left_trapezoidal(L,one);
        print("K = %\n", str(res));
        assert(approx(res, mat(3,3,
            R, CZERO, CZERO,
            -I, R, CZERO, 
            -3+3*I, -3-I, R
        )));
    }

    {
        M := mat(3,3,
            2*R, I, -2*R,
            CZERO, -I, CZERO,
            CZERO, CZERO, 1+I
        );
        one := mat_id(Complex64, 3, 3);
        res := solve_linear_right_trapezoidal(M,one);
        print("K = %\n", str(res));
        assert(approx(res, mat(3,3,
            0.5*R, 0.5*R, 0.5-0.5*I,
            CZERO, I, CZERO,
            CZERO, CZERO, 0.5-0.5*I
        )));
    }

    {
        M := mat(3, 3,
            2*R, I, -2*R,
            CZERO, -I, 2*I,
            R, CZERO, R+I
        );
        iM := inverse(M);
        print("iM = %\n", str(iM));
        assert(approx(iM, mat(3,3,
            c64(0.25, 0.25), c64(0.25, 0.25), c64(0.5, -0.5),
            -0.5*R, -0.5*R+I, R,
            -0.25*R, -0.25*R, 0.5*R
        )));
    }

    {
        M := mat(3, 3,
            2*R, I, -2*R,
            CZERO, -I, 2*I,
            R, CZERO, R+I
        );
        d := determinant(M);
        print("det(M) = %\n", d);
        assert(d == -4*I);
    }
}