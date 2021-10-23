
//  ######     ##       #### ##    ##     #######  ##     ##  #######
// ##    ##    ##        ##  ###   ##    ##     ##  ##   ##  ##     ##
// ##          ##        ##  ####  ##           ##   ## ##          ##
//  ######     ##        ##  ## ## ##     #######     ###     #######
//       ##    ##        ##  ##  ####    ##          ## ##   ##
// ##    ##    ##        ##  ##   ###    ##         ##   ##  ##
//  ######     ######## #### ##    ##    ######### ##     ## #########

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.1 2x2 Linear System Solve, page 153
// ? solves the system `a x = b` with 2x2 matrix `a` and vectors `x`, `b`
sl2 :: solve_linear_2x2;
solve_linear_2x2 :: inline (a: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_2x2_default(a,b,x);
}

#scope_file
solve_linear_2x2_default :: (a: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(rows(a) == 2, "dimension mismatch: 2 != A.N = %\n", A.N);
        assert(cols(a) == 2, "dimension mismatch: 2 != A.M = %\n", A.M);
        assert(dim(b) == 2, "dimension mismatch: 2 != B.N = %\n", B.N);
        assert(dim(x) == 2, "dimension mismatch: 2 != C.N = %\n", C.N);
    }

    sigma := max(abs(get(a,0,0)), abs(get(a,0,1)), abs(get(a,1,0)), abs(get(a,1,1)));

    mu := a / sigma;
    det := get(mu,0,0) * get(mu,1,1) - get(mu,0,1) * get(mu,1,0);
    g := b / sigma / det;

    m := dmat(2,2, get(mu,1,1), -get(mu,0,1), -get(mu,1,0), get(mu,0,0));
    mat_mul_vec(m, g, x);
}
#scope_export





//  ######     ##       #### ##    ##    ##       ########    ######## ########  ####
// ##    ##    ##        ##  ###   ##    ##       ##             ##    ##     ##  ##
// ##          ##        ##  ####  ##    ##       ##             ##    ##     ##  ##
//  ######     ##        ##  ## ## ##    ##       ######         ##    ########   ##
//       ##    ##        ##  ##  ####    ##       ##             ##    ##   ##    ##
// ##    ##    ##        ##  ##   ###    ##       ##             ##    ##    ##   ##
//  ######     ######## #### ##    ##    ######## ########       ##    ##     ## ####

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.2 Forward Linear Recurrence, page 155
// ? solve the system Lx = b where L is a unit left-triangluar matrix
// ? that means, all diagonals are 1 and all L_ij = 0 for j > i.
sllta :: solve_linear_left_triangular;
solve_linear_left_triangular :: inline (L: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {    
    solve_linear_left_triangular_default(L,b,x);
}

#scope_file
solve_linear_left_triangular_default :: (L: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(rows(L) == cols(L), "L is not quadratic L(%,%)\n", rows(L), cols(L));
        assert(rows(L) == dim(b), "dimensions mismatch L(%,%), b(%)\n", rows(L), cols(L), dim(b));
        assert(cols(L) == dim(x), "dimensions mismatch L(%,%), x(%)\n", rows(L), cols(L), dim(x));
        assert(is_diagonal_unit(L), "L is not unit: Diagonal entries are not 1\n");
        assert(is_left_triangular(L), "L is not left-triangular\n");
    }
    
    for i: 0..dim(b)-1 {
        sum := zero(C.data_type);
        for j: 0..i-1 {
            sum += get(L,i,j) * get(x,j);
        }
        set(x,i, get(b,i) - sum);
    }
}
#scope_export


//  ######     ##       #### ##    ##    ########  ####    ######## ########  ####
// ##    ##    ##        ##  ###   ##    ##     ##  ##        ##    ##     ##  ##
// ##          ##        ##  ####  ##    ##     ##  ##        ##    ##     ##  ##
//  ######     ##        ##  ## ## ##    ########   ##        ##    ########   ##
//       ##    ##        ##  ##  ####    ##   ##    ##        ##    ##   ##    ##
// ##    ##    ##        ##  ##   ###    ##    ##   ##        ##    ##    ##   ##
//  ######     ######## #### ##    ##    ##     ## ####       ##    ##     ## ####

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.4 Backward Linear Recurrence, page 158
// ? solve the system Rx = b where R is a right-triangluar matrix
// ? that means all R_ij = 0 for j < i.
slrta :: solve_linear_right_triangular;
solve_linear_right_triangular :: inline (R: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_right_triangular_default(R,b,x);
}

#scope_file
solve_linear_right_triangular_default :: (R: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(rows(R) == cols(R), "L is not quadratic L(%,%)\n", rows(R), cols(R));
        assert(rows(R) == dim(b), "dimensions mismatch L(%,%), b(%)\n", rows(R), cols(R), dim(b));
        assert(cols(R) == dim(x), "dimensions mismatch L(%,%), x(%)\n", rows(R), cols(R), dim(x));
        assert(is_right_triangular(R), "R is not right-triangular\n");
    }

    i : int;
    for 0..dim(b)-1 {
        i = dim(b)-1 - it;
        sum := zero(x.data_type);
        for j: i+1..dim(b)-1 {
            sum += get(R,i,j)*get(x,j);
        }
        set(x,i, (get(b,i) - sum)/get(R,i,i));
    }
}
#scope_export





//  ######     ##       #### ##    ##    ########  ####     #######  ######## ########  ####
// ##    ##    ##        ##  ###   ##    ##     ##  ##     ##     ##    ##    ##     ##  ##
// ##          ##        ##  ####  ##    ##     ##  ##     ##     ##    ##    ##     ##  ##
//  ######     ##        ##  ## ## ##    ########   ##     ##     ##    ##    ########   ##
//       ##    ##        ##  ##  ####    ##   ##    ##     ##  ## ##    ##    ##   ##    ##
// ##    ##    ##        ##  ##   ###    ##    ##   ##     ##    ##     ##    ##    ##   ##
//  ######     ######## #### ##    ##    ##     ## ####     ##### ##    ##    ##     ## ####

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.6 Solve Right Quasi-Triangular System, page 164
// ? solve the system `Rx=b` where `R` is a right(upper) quasi-triangular matrix and `x`, `b` are vectors.
slrqta :: solve_linear_right_quasi_triangular;
solve_linear_right_quasi_triangular :: inline (R: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_right_quasi_triangular_default(R,b,x);
}

#scope_file
solve_linear_right_quasi_triangular_default :: (R: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(rows(R) == cols(R), "L is not quadratic L(%,%)\n", rows(R), cols(R));
        assert(rows(R) == dim(b), "dimensions mismatch L(%,%), b(%)\n", rows(R), cols(R), dim(b));
        assert(cols(R) == dim(x), "dimensions mismatch L(%,%), x(%)\n", rows(R), cols(R), dim(x));
        assert(is_right_quasi_triangular(R), "R is not right quasi-triangular!\n");
    }
    U :: C.data_type;

    i := dim(b)-1; // start with last row and wrap from bottom to up
    while i >= 0 {
        // check if 2x2 block by looking at left column
        // the term `|| i == 0` was missing
        if (i > 0 && get(R,i,i-1) == zero(U)) || i == 0 {
            sum := zero(U);
            for j: i+1..dim(b)-1 {
                sum += get(R,i,j) * get(x,j);
            }
            set(x,i, (get(b,i) - sum)/get(R,i,i));
            i -= 1;
        } else if i > 0 { // so now we have a 2x2 block where get(R,i,i) is the bottom-right element
            a : DenseMatrix(U, 2, 2) = ---;
            mat_init(*a, 
                get(R,i-1,i-1), get(R,i-1,i),
                get(R,i,i-1), get(R,i,i)
            );
            sum1 := zero(U);
            sum2 := zero(U);
            for j: i+1..dim(b)-1 {
                sum1 += get(R,i-1,j) * get(x,j);
                sum2 += get(R,i,j) * get(x,j);
            }
            c := dvec(2, get(b,i-1)-sum1, get(b,i)-sum2);
            tmp : DenseVector(U, 2, 2);
            solve_linear_2x2(a, c, *tmp);
            set(x, i-1, get(tmp,0));
            set(x, i,   get(tmp,1));
            i -= 2;
        }
    }
}
#scope_export





//  ######     ##       #### ##    ##    ##     ## ##    ## ####
// ##    ##    ##        ##  ###   ##    ##     ## ###   ##  ##
// ##          ##        ##  ####  ##    ##     ## ####  ##  ##
//  ######     ##        ##  ## ## ##    ##     ## ## ## ##  ##
//       ##    ##        ##  ##  ####    ##     ## ##  ####  ##
// ##    ##    ##        ##  ##   ###    ##     ## ##   ###  ##
//  ######     ######## #### ##    ##     #######  ##    ## ####

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.8 Orthogonal Projection, page 167
solve_linear_unitary :: solve_linear_orthogonal_projection;
slop                 :: solve_linear_orthogonal_projection;
solve_linear_orthogonal_projection :: inline (Q: *$A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_orthogonal_projection_default(Q,b,x);
}

#scope_file
solve_linear_orthogonal_projection_default :: (Q: *$A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(cols(Q) == dim(b), "dimensions mismatch (%,%), (%)\n", rows(Q), cols(Q), dim(b));
        assert(cols(Q) == dim(x), "dimensions mismatch (%,%), (%)\n", rows(Q), cols(Q), dim(x));
        assert(is_unitary(Q));
    }
    dagger(Q);
    mat_mul_vec(<<Q, b, x);
}
#scope_export




//  ######     ##       #### ##    ##     ######  ##     ##  ######   ######      #######  ########  ######## ##     ##    ########  ########        ##
// ##    ##    ##        ##  ###   ##    ##    ## ##     ## ##    ## ##    ##    ##     ## ##     ##    ##    ##     ##    ##     ## ##     ##       ##
// ##          ##        ##  ####  ##    ##       ##     ## ##       ##          ##     ## ##     ##    ##    ##     ##    ##     ## ##     ##       ##
//  ######     ##        ##  ## ## ##     ######  ##     ## ##       ##          ##     ## ########     ##    #########    ########  ########        ##
//       ##    ##        ##  ##  ####          ## ##     ## ##       ##          ##     ## ##   ##      ##    ##     ##    ##        ##   ##   ##    ##
// ##    ##    ##        ##  ##   ###    ##    ## ##     ## ##    ## ##    ##    ##     ## ##    ##     ##    ##     ##    ##        ##    ##  ##    ##
//  ######     ######## #### ##    ##     ######   #######   ######   ######      #######  ##     ##    ##    ##     ##    ##        ##     ##  ######

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.9 Successive Orthogonal Projection, page 168
// ? not sure if this algo is useful. seems to be slower than the previous one.
slsop :: solve_linear_successive_orthogonal_projection;
solve_linear_successive_orthogonal_projection :: inline (Q: $A/MatrixType, b: *$B/VectorType, x: *$C/VectorType) {
    solve_linear_successive_orthogonal_projection(Q,b,x);
}

#scope_file
solve_linear_successive_orthogonal_projection_default :: (Q: $A/MatrixType, b: *$B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(cols(Q) == dim(b), "dimensions mismatch (%,%), (%)\n", rows(Q), cols(Q), dim(b));
        assert(cols(Q) == dim(x), "dimensions mismatch (%,%), (%)\n", rows(Q), cols(Q), dim(x));
        assert(is_unitary(Q));
    }

    for i: 0..rows(Q)-1 {

        q := mat_row_heap_view(*Q, i, cols(Q));
        cq_dot_b := q * b; // ? the conjugation is on by default for the dot product
        set(x, i, cq_dot_b);
        saxpy(b, q, -cq_dot_b);

        // q := row(Q, i);
        // set(x,i, conjugate(q) * b); // not sure if this is performant...
        // b = b - q * get(x,i);
    }
}
#scope_export



//  ######     ##       #### ##    ##    ##       ########    ######## ########  ########  ########
// ##    ##    ##        ##  ###   ##    ##       ##             ##    ##     ## ##     ##      ##
// ##          ##        ##  ####  ##    ##       ##             ##    ##     ## ##     ##     ##
//  ######     ##        ##  ## ## ##    ##       ######         ##    ########  ########     ##
//       ##    ##        ##  ##  ####    ##       ##             ##    ##   ##   ##          ##
// ##    ##    ##        ##  ##   ###    ##       ##             ##    ##    ##  ##         ##
//  ######     ######## #### ##    ##    ######## ########       ##    ##     ## ##        ########

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 3.4.5.1 Left Trapezoidal, page 159
slltz :: solve_linear_left_trapezoidal;

solve_linear_left_trapezoidal :: inline (L: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_left_trapezoidal(L,b,x);
}

#scope_file
solve_linear_left_trapezoidal_default :: (L: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(rows(L) == dim(b), "dimensions mismatch (%,%), (%)\n", rows(L), cols(L), dim(b));
        assert(cols(L) == dim(x), "dimensions mismatch (%,%), (%)\n", rows(L), cols(L), dim(x));
        assert(is_left_trapezoidal(L)); // more rows than columns N >= M
    }

    L2 := mat_heap_view(*L, cols(L), cols(L));
    b2 := vec_heap_view(*b, cols(L));
    solve_linear_left_triangular(L2, b2, x);
}
#scope_export


solve_linear_left_trapezoidal :: inline (l: $A/MatrixType, m: $B/MatrixType, x: *$C/MatrixType) {
    solve_linear_left_trapezoidal_default(l,m,x);
}

#scope_file
solve_linear_left_trapezoidal_default :: (l: $A/MatrixType, m: $B/MatrixType, x: *$C/MatrixType) {
    #if CHECKS {
        assert(rows(l) == rows(m), "dimensions mismatch l(%,%), m(%,%)\n", rows(l), cols(l), rows(b), cols(b));
        assert(cols(l) == rows(x), "dimensions mismatch l(%,%), x(%,%)\n", rows(l), cols(l), rows(x), cols(x));
        assert(cols(m) == cols(x), "dimensions mismatch m(%,%), x(%,%)\n", rows(m), cols(m), rows(x), cols(x));
    }

    for j: 0..cols(m)-1 {
        b := mat_col_heap_view(*m, j, rows(m));
        y := mat_col_heap_view( x, j, rows(x));
        solve_linear_left_trapezoidal(l, b, *y);
    }
}
#scope_export


//  ######     ##       #### ##    ##    ########  ####    ######## ########  ########  ########
// ##    ##    ##        ##  ###   ##    ##     ##  ##        ##    ##     ## ##     ##      ##
// ##          ##        ##  ####  ##    ##     ##  ##        ##    ##     ## ##     ##     ##
//  ######     ##        ##  ## ## ##    ########   ##        ##    ########  ########     ##
//       ##    ##        ##  ##  ####    ##   ##    ##        ##    ##   ##   ##          ##
// ##    ##    ##        ##  ##   ###    ##    ##   ##        ##    ##    ##  ##         ##
//  ######     ######## #### ##    ##    ##     ## ####       ##    ##     ## ##        ########

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 3.4.5.2 Right Trapezoidal, page 160
// ? There are degrees of freedom left. These are set to 0 here.
slrtz :: solve_linear_right_trapezoidal;
solve_linear_right_trapezoidal :: inline (R: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_right_trapezoidal_default(R,b,x);
}

#scope_file
solve_linear_right_trapezoidal_default :: (R: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(rows(R) == dim(b), "dimensions mismatch (%,%), (%)\n", rows(R), cols(R), dim(b));
        assert(cols(R) == dim(x), "dimensions mismatch (%,%), (%)\n", rows(R), cols(R), dim(x));
        assert(is_right_trapezoidal(R), "Matrix is not right trapezoidal: M(%) >= N(%)?\n", cols(R), rows(R));
    }

    R2 := mat_heap_view(*R, rows(R), rows(R));
    b2 := vec_heap_view(*b, rows(R));
    x2 := vec_heap_view(x, rows(R));
    solve_linear_right_triangular(R2,b2,*x2);
}
#scope_export


solve_linear_right_trapezoidal :: inline (r: $A/MatrixType, m: $B/MatrixType, x: *$C/MatrixType) {
    solve_linear_right_trapezoidal_default(r,m,x);
}

#scope_file
solve_linear_right_trapezoidal_default :: (r: $A/MatrixType, m: $B/MatrixType, x: *$C/MatrixType) {
    #if CHECKS {
        assert(rows(r) == rows(m), "dimensions mismatch (%,%), (%,%)\n", rows(r), cols(r), rows(m), cols(m));
        assert(cols(r) == rows(x), "dimensions mismatch (%,%), (%,%)\n", rows(r), cols(r), rows(x), cols(x));
        assert(cols(x) == cols(m), "dimensions mismatch (%,%), (%,%)\n", rows(x), cols(x), rows(m), cols(m));
    }


    for j: 0..cols(m)-1 {
        b := mat_col_heap_view(*m, j, rows(m));
        y := mat_col_heap_view( x, j, rows(x));
        solve_linear_right_trapezoidal(r, b, *y);
    }
}
#scope_export




// ##       ########
// ##       ##     ##
// ##       ##     ##
// ##       ########
// ##       ##   ##
// ##       ##    ##
// ######## ##     ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.1 Gaussian Factorization - No Pivoting, page 199
decompose_LR :: (m: $A/MatrixType) -> L: A #must, R: A #must {
    L : A;
    mat_init_unit(*L);
    R : A;

    for m {
        if it_column >= it_row {
            set(*R, it_row, it_column, it);
        } else {
            set(*L, it_row, it_column, it);
        }
    }

    return L, R;
}







// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.1 Gaussian Factorization - No Pivoting, page 199
gaussian_factorization_no_pivot :: (m: $A/MatrixType) -> A #must {
    a := copy(m); // TODO isn't the matrix passed by value anyways??

    for k: 0..min(A.N,A.M)-1 {
        assert(abs(get(a,k,k)) != zero(A.data_type), "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        for i: k+1..A.N-1 {
            set(*a,i,k, get(a,i,k)/get(a,k,k));
        }
        for j: k+1..A.M-1 {
            for i: k+1..A.N-1 {
                set(*a,i,j, get(a,i,j) - get(a,i,k) * get(a,k,j));
            }
        }
    }
    return a;
}




// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.2 Solve Linear System via Gaussian Factorization - No Pivoting, page 199
slgf_np :: solve_linear_gaussian_factorization_no_pivot;
solve_LRnp :: solve_linear_gaussian_factorization_no_pivot;
solve_LUnp :: solve_linear_gaussian_factorization_no_pivot;
solve_linear_gaussian_factorization_no_pivot :: (a: $A/MatrixType, b: $B/VectorType) -> VectorHelper(A,B,A.M).T {
    M := gaussian_factorization_no_pivot(a);
    L, R := decompose_LR(M);
    y := solve_linear_left_trapezoidal(L, b);
    return solve_linear_right_trapezoidal(R, y);
}








// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.3 Gaussian Factorization with Full Pivoting, page 204
gaussian_factorization_full_pivot :: (m: $A/MatrixType)
 -> A #must, 
    rows: [A.N]int #must, 
    cols: [A.M]int #must, 
    row_permutations: int, 
    col_permutations: int 
{
    // println("gaussian_factorization_full_pivot");
    a := copy(m);
    rows : [A.N]int;
    cols : [A.M]int;
    row_permutations := 0;
    col_permutations := 0;

    for 0..A.N-1 {
        rows[it] = it;
    }
    for 0..A.M-1 {
        cols[it] = it;
    }

    // println("    entering for k");
    for k: 0..min(A.N,A.M)-1 {
        ii := k;
        jj := k;
        rho := abs(get(a,k,k));

        // println("    setting ii, jj");
        for j: k..A.M-1 {
            for i: k..A.N-1 {
                tmp := abs(get(a,i,j));
                if tmp > rho {
                    ii = i;
                    jj = j;
                    rho = tmp;
                }
            }
        }

        assert(rho != 0.0, "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        // instead of swapping, use the permutation arrays rows and cols to access the 
        // entries of a. Should be faster?
        // println("    permute rows");
        if ii > k {
            // println("        swap_rows");
            swap_rows(*a, k, ii);
            // println("        swap row array entries");
            rows[k], rows[ii] = swap(rows[k], rows[ii]);
            row_permutations += 1;
        }
        // println("    permute columns");
        if jj > k {
            swap_columns(*a, k, jj);
            cols[k], cols[jj] = swap(cols[k], cols[jj]);
            col_permutations += 1;
        }

        // if k == min(A.N, A.M)-1 then break;

        // println("    updating matrix");
        for i: k+1..A.N-1 {
            set(*a,i,k, get(a,i,k)/get(a,k,k));
        }
        // println("    updating matrix 2");
        for j: k+1..A.M-1 {
            for i: k+1..A.N-1 {
                set(*a,i,j, get(a,i,j) - get(a,i,k) * get(a,k,j));
            }
        }
    }

    // println("    returning results");
    return a, rows, cols, row_permutations, col_permutations;
}








// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.4 Solve Linear System via Gaussian Factorization with Full Pivoting, page 206
slgf_fp :: solve_linear_gaussian_factorization_full_pivot;
solve_LRfp :: solve_linear_gaussian_factorization_full_pivot;
solve_LUfp :: solve_linear_gaussian_factorization_full_pivot;
solve_linear_gaussian_factorization_full_pivot :: (a: $A/MatrixType, b: $B/VectorType) -> VectorHelper(A,B,A.N).T {
    M, rows, cols := gaussian_factorization_full_pivot(a);
    L, R := decompose_LR(M);

    y := solve_linear_left_trapezoidal(L, permute(b, ..rows));
    w := solve_linear_right_trapezoidal(R, y);
    return permute(w, ..cols);
}










// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7. Gaussian Factorization with Partial Pivoting, page 210
gaussian_factorization_row_pivot :: (m: $A/MatrixType) -> A #must, rows: [A.N]int, permutations: int {
    a := copy(m);

    rows : [A.N]int = ---;
    permutations := 0;

    for 0..A.N-1 {
        rows[it] = it;
    }

    for k: 0..#run min(A.N,A.M)-1 {
        ii := k;
        rho := abs(get(a,k,k));

        for i: k..A.N-1 {
            tmp := abs(get(a,i,k));
            if tmp > rho {
                ii = i;
                rho = tmp;
            }
        }

        assert(rho != 0.0, "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        // instead of swapping, use the permutation arrays rows and cols to access the 
        // entries of a. Should be faster?
        if ii > k {
            swap_rows(*a, k, ii);
            rows[k], rows[ii] = swap(rows[k], rows[ii]);
            permutations += 1;
        }

        for i: k+1..A.N-1 {
            set(*a,i,k, get(a,i,k)/get(a,k,k));
        }
        for j: k+1..A.M-1 {
            for i: k+1..A.N-1 {
                set(*a,i,j, get(a,i,j) - get(a,i,k) * get(a,k,j));
            }
        }
    }

    return a, rows, permutations;
}



// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.7 Solve Linear System via Gaussian Factorization with Partial Pivoting, page 215
slgf_rp :: solve_linear_gaussian_factorization_row_pivot;
solve_LRrp :: solve_linear_gaussian_factorization_row_pivot;
solve_LUrp :: solve_linear_gaussian_factorization_row_pivot;
solve_linear_gaussian_factorization_row_pivot :: (a: $A/MatrixType, b: $B/VectorType) -> VectorHelper(A,B,A.N).T {
    M, rows := gaussian_factorization_row_pivot(a);
    L, R := decompose_LR(M);

    y := solve_linear_left_trapezoidal(L, permute(b, ..rows));
    return solve_linear_right_trapezoidal(R, y);
}






// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.10.1 Matrix Inverse, page 261
inverse :: (a: $A/MatrixType) -> A {
    M, rows := gaussian_factorization_row_pivot(a);
    L, R := decompose_LR(M);

    id : A;
    mat_init_unit(*id);

    Y := solve_linear_left_trapezoidal(L, id);
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
determinant :: (a: $A/MatrixType) -> A.data_type {
    #if CHECKS {
        #run assert(A.N == A.M, "matrix not quadratic!\n (%,%)", A.N, A.M);
    }
    // TODO implement optimized formula for 3x3 matrix
    #if A.N == 2 {
        return get(a,0,0)*get(a,1,1) - get(a,0,1)*get(a,1,0);
    } else {    
        M, rows, permutations := gaussian_factorization_row_pivot(a);
        L, R := decompose_LR(M);
    
        p := one(A.data_type);
        for 0..R.N-1 {
            p *= get(R,it,it);
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
    println_push("Linear Algebra", color = .FG_WHITE);

    {
        println_push("test 1", color = .FG_GREEN);

        A := dmat(2,2, R, R, 2.0 - I, I);
        println("mat_init: %", pstr(A));

        b := dvec(2, I, 2*R+4*I);
        println("vec_init: %", pstr(b));

        res := solve_linear_2x2(A, b);
        println("result: %", pstr(res));

        sol := dvec(2,  // ? this is the correct answer
            -0.25 + 1.75*I,
            0.25 - 0.75*I
        );
        assert(res == sol);
    }

    {
        println_push("test 2", color = .FG_GREEN);

        L : DenseMatrix(C64, 3, 3) = ---;
        mat_init(*L,
            R, C_ZERO, C_ZERO,
            3.0*R, R, C_ZERO,
            I, R+I, R
        );
        b : DenseVector(C64, 3) = ---;
        vec_init(*b, 3.0*R, -I, R-2.0*I);
        res := solve_linear_left_triangular(L, b);
        println("% = % ?", pstr(L * res), pstr(b));
        assert(L * res == b);
    }

    {
        println_push("test 3", color = .FG_GREEN);

        M := dmat(3,3,
            2.0*R, I, R-2.0*I,
            C_ZERO, R, 3.0*I,
            C_ZERO, C_ZERO, -3.0*R
        );
        b := dvec(3, 3.0*R, -I, R-2.0*I);
        res := solve_linear_right_triangular(M, b);
        println("res %", pstr(res));

        test := norm_2(M * res - b);
        println("%", pstr(M * res - b));
        println("test = %", test); 
        assert(test < 1.0e-8, "res = %, M * res = % =? b = %", res, M*res, b);
    }

    {
        println_push("test 4", color = .FG_GREEN);

        m : DenseMatrix(int, 4, 4) = ---;
        mat_init(*m,
            1, 3, 4, 5,
            0, 2, 3, 4,
            0, 1, 2, 3,
            0, 0, 0, 2
        );
        cm := mat_cast(C64, m);
        
        b : DenseVector(float64, 4) = ---;
        vec_init(*b, 2,3,4,5);
        cb := vec_cast(C64, b);

        println("cm %", pstr(cm));
        println("cb %", pstr(cb));

        res := solve_linear_right_quasi_triangular(cm,cb);
        println("res = %", pstr(res));

        test := norm_2(cm * res - cb);
        assert(test < 1.0e-8);
    }

    {
        println_push("test 5", color = .FG_GREEN);

        M : DenseMatrix(C64, 3, 2) = ---;
        mat_init(*M, 
            R, C_ZERO, 
            -R, R, 
            -R, -R
        );
        b : DenseVector(C64, 3) = ---;
        vec_init(*b, R, 2.0*R, -4.0*R);

        res := solve_linear_left_trapezoidal(M, b);
        println("sllt: %", pstr(res));
        assert( approx(M*res, b) );
    }

    {
        println_push("test 6", color = .FG_GREEN);

        M : DenseMatrix(C64, 2, 3) = ---;
        mat_init(*M, 
            R, R, C_ZERO, 
            C_ZERO, R, C_ZERO 
        );
        b : DenseVector(C64, 2) = ---;
        vec_init(*b,4*R, 6*R);

        res := solve_linear_right_trapezoidal(M, b);
        println("slrt: %", pstr(res));
        assert( approx(M*res, b) );
    }
    {
        println_push("test 7", color = .FG_GREEN);

        M   : DenseMatrix(C64, 2, 2) = ---;
        res : DenseMatrix(C64, 2, 2) = ---;
        mat_init(*M,
            2*R, 3*R,
            -2*R, R
        );
        mat_init(*res,
            2*R, 3*R,
            -1*R, 4*R
        );
        A := gaussian_factorization_no_pivot(M);
        println("gaussian factorization: %", pstr(M));
        println("-> %", pstr(A));
        assert( A == res , "%\n != %\n!", str(A), str(res));
    }

    {
        println_push("test 8", color = .FG_GREEN);

        M := dmat(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := dvec(2, 4*R, 6*I);

        res := solve_LRnp(M, b);
        println("solve linear system via Gaussian factorization (no pivoting): %", pstr(res));
        assert(res == dvec(2, 2.5-4.5*I, 4.5-1.5*I), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        println_push("test 9", color = .FG_GREEN);

        M : DenseMatrix(C64, 3, 3) = ---;
        mat_init(*M,
            C_ZERO, 4*R, -3*R,
            2*R, C_ZERO, 4*R,
            3*R, -5*R, 3*R
        );
        println("M %", pstr(M));

        A, rows, cols := gaussian_factorization_full_pivot(M);
        println("A = %", pstr(A));
    }

    {
        println_push("test 10", color = .FG_GREEN);

        M : DenseMatrix(C64, 3, 3) = ---;
        mat_init(*M,
            C_ZERO, 4*R, -3*R,
            2*R, C_ZERO, 4*R,
            3*R, -5*R, 3*R
        );

        A, rows, cols := gaussian_factorization_full_pivot(M);
        L, R := decompose_LR(A);

        println("L = %", pstr(L));
        println("R = %", pstr(R));
    }

    {
        println_push("test 11", color = .FG_GREEN);

        M : DenseMatrix(C64, 2, 2) = ---;
        mat_init(*M,
            R, I,
            -3*R, 2.0 - I
        );
        b : DenseVector(C64, 2) = ---;
        vec_init(*b, 4*R, 6*I);

        res := solve_LRfp(M, b);
        println("solve linear system via Gaussian factorization (full pivoting): %", pstr(res));
        sol : DenseVector(C64, 2) = ---;
        vec_init(*sol, 2.5-4.5*I, 4.5-1.5*I);
        assert(approx(res, sol), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        println_push("test 12", color = .FG_GREEN);

        M := dmat(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := dvec(2,   4*R, 6*I);

        res := solve_LRrp(M, b);
        println("solve linear system via Gaussian factorization (row pivoting): %", pstr(res));
        assert(approx(res, dvec(2, 2.5-4.5*I, 4.5-1.5*I)), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        println_push("test 13", color = .FG_GREEN);

        L := dmat(3,3,
            R, C_ZERO, C_ZERO,
            I, R, C_ZERO,
            2*R, 3+I, R
        );
        one : DenseMatrix(C64, 3, 3);
        mat_init_unit(*one);
        res := solve_linear_left_trapezoidal(L,one);
        println("K = %", pstr(res));
        assert(approx(res, dmat(3,3,
            R, C_ZERO, C_ZERO,
            -I, R, C_ZERO, 
            -3+3*I, -3-I, R
        )));
    }

    {
        println_push("test 14", color = .FG_GREEN);

        M := dmat(3,3,
            2*R, I, -2*R,
            C_ZERO, -I, C_ZERO,
            C_ZERO, C_ZERO, 1+I
        );
        one : DenseMatrix(C64, 3, 3);
        mat_init_unit(*one);
        res := solve_linear_right_trapezoidal(M,one);
        println("K = %", pstr(res));
        assert(approx(res, dmat(3,3,
            0.5*R, 0.5*R, 0.5-0.5*I,
            C_ZERO, I, C_ZERO,
            C_ZERO, C_ZERO, 0.5-0.5*I
        )));
    }

    {
        println_push("test 15", color = .FG_GREEN);

        M := dmat(3, 3,
            2*R, I, -2*R,
            C_ZERO, -I, 2*I,
            R, C_ZERO, R+I
        );
        iM := inverse(M);
        println("iM = %", pstr(iM));
        assert(approx(iM, dmat(3,3,
            0.25*R + 0.25*I, 0.25*R+0.25*I, 0.5*R-0.5*I,
            -0.5*R, -0.5*R+I, R,
            -0.25*R, -0.25*R, 0.5*R
        )));
    }

    {
        println_push("test 16", color = .FG_GREEN);

        M := dmat(3, 3,
            2*R, I, -2*R,
            C_ZERO, -I, 2*I,
            R, C_ZERO, R+I
        );
        d := determinant(M);
        println("det(M) = %", str(d));
        assert(d == -4*I);
    }
}