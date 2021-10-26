
//  ######  ##        #######
// ##    ## ##       ##     ##
// ##       ##              ##
//  ######  ##        #######
//       ## ##       ##
// ##    ## ##       ##
//  ######  ######## #########

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
        assert(rows(a) == 2, "dimension mismatch: 2 != A.N = %\n", rows(a));
        assert(cols(a) == 2, "dimension mismatch: 2 != A.M = %\n", cols(a));
        assert(dim(b) == 2, "dimension mismatch: 2 != B.N = %\n", dim(b));
        assert(dim(x) == 2, "dimension mismatch: 2 != C.N = %\n", dim(x));
    }

    // TODO remove unnecesary allocations!

    sigma := max(abs(get(a,0,0)), abs(get(a,0,1)), abs(get(a,1,0)), abs(get(a,1,1)));
    mu : DenseMatrix(A.data_type, 2, 2) = ---;
    copy(*mu, a);
    div(*mu, sigma);
    det := get(mu,0,0) * get(mu,1,1) - get(mu,0,1) * get(mu,1,0);

    g : DenseVector(B.data_type, 2) = ---;
    copy(*g, b);
    div(*g, sigma*det);
    m := dmat(2,2, get(mu,1,1), -get(mu,0,1), -get(mu,1,0), get(mu,0,0));
    mat_mul_vec(m, g, x);
}
#scope_export





//  ######  ##       ##       ########    ###
// ##    ## ##       ##          ##      ## ##
// ##       ##       ##          ##     ##   ##
//  ######  ##       ##          ##    ##     ##
//       ## ##       ##          ##    #########
// ##    ## ##       ##          ##    ##     ##
//  ######  ######## ########    ##    ##     ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.2 Forward Linear Recurrence, page 155
// ? solve the system Lx = b where L is a unit left-triangluar matrix
// ? that means, all diagonals are 1 and all L_ij = 0 for j > i.

// TODO clarify, does this work with sllta(L, x, *x) ?
// TODO it seems, that the inner loop only traverses already updated data, so it should be possible
// ? The diagonals L_ii = 1, but they are never used since j < i at all times!

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

        // ? not necessary, since these entries are not used anyways
        // assert(is_diagonal_unit(L), "L is not unit: Diagonal entries are not 1\n");
        // assert(is_left_triangular(L), "L is not left-triangular\n");
    }

    // println_push("solve_linear_left_triangular_default");
    
    for i: 0..dim(b)-1 {
        sum := zero(C.data_type);
        for j: 0..i-1 {
            sum += get(L,i,j) * get(x,j);
        }
        set(x,i, get(b,i) - sum);
    }
}
#scope_export



//  ######  ##       ########  ########    ###
// ##    ## ##       ##     ##    ##      ## ##
// ##       ##       ##     ##    ##     ##   ##
//  ######  ##       ########     ##    ##     ##
//       ## ##       ##   ##      ##    #########
// ##    ## ##       ##    ##     ##    ##     ##
//  ######  ######## ##     ##    ##    ##     ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.4.4 Backward Linear Recurrence, page 158
// ? solve the system Rx = b where R is a right-triangluar matrix
// ? that means all R_ij = 0 for j < i.

// TODO clarify, does this work with sllta(L, x, *x) ?
// TODO It seems, that the inner loop only traverses already updated data, so it should be possible.
// TODO Note the reversed order!
// ? The diagonals R_ii = 0, but they are never used since j > i at all times!

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

        // ? not necessary, since these entries are not used anyways
        // assert(is_right_triangular(R), "R is not right-triangular\n");
    }

    i : int;
    d := dim(b);
    for 0..d-1 {
        i = d-1 - it; // ? reverse order!
        sum := zero(x.data_type);
        for j: i+1..d-1 {
            sum += get(R,i,j)*get(x,j);
        }
        set(x,i, (get(b,i) - sum)/get(R,i,i));
    }
}
#scope_export




//  ######  ##       ########   #######  ########    ###
// ##    ## ##       ##     ## ##     ##    ##      ## ##
// ##       ##       ##     ## ##     ##    ##     ##   ##
//  ######  ##       ########  ##     ##    ##    ##     ##
//       ## ##       ##   ##   ##  ## ##    ##    #########
// ##    ## ##       ##    ##  ##    ##     ##    ##     ##
//  ######  ######## ##     ##  ##### ##    ##    ##     ##

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
            tmp : DenseVector(U, 2);
            solve_linear_2x2(a, c, *tmp);
            set(x, i-1, get(tmp,0));
            set(x, i,   get(tmp,1));
            i -= 2;
        }
    }
}
#scope_export




//  ######  ##        #######  ########
// ##    ## ##       ##     ## ##     ##
// ##       ##       ##     ## ##     ##
//  ######  ##       ##     ## ########
//       ## ##       ##     ## ##
// ##    ## ##       ##     ## ##
//  ######  ########  #######  ##

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



//  ######  ##        ######   #######  ########
// ##    ## ##       ##    ## ##     ## ##     ##
// ##       ##       ##       ##     ## ##     ##
//  ######  ##        ######  ##     ## ########
//       ## ##             ## ##     ## ##
// ##    ## ##       ##    ## ##     ## ##
//  ######  ########  ######   #######  ##

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



//  ######  ##       ##       ######## ########
// ##    ## ##       ##          ##         ##
// ##       ##       ##          ##        ##
//  ######  ##       ##          ##       ##
//       ## ##       ##          ##      ##
// ##    ## ##       ##          ##     ##
//  ######  ######## ########    ##    ########

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 3.4.5.1 Left Trapezoidal, page 159
slltz :: solve_linear_left_trapezoidal;

solve_linear_left_trapezoidal :: inline (L: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_left_trapezoidal_default(L,b,x);
}

#scope_file
solve_linear_left_trapezoidal_default :: (L: $A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    #if CHECKS {
        assert(rows(L) == dim(b), "dimensions mismatch (%,%), (%)\n", rows(L), cols(L), dim(b));
        assert(cols(L) == dim(x), "dimensions mismatch (%,%), (%)\n", rows(L), cols(L), dim(x));

        // ? this test should be disabled since `solve_linear_left_triangular` works with any matrix and prevents unnecessary allocations
        // assert(is_left_trapezoidal(L)); // more rows than columns N >= M
    }

    // println_push("solve_linear_left_trapezoidal_default");

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
        assert(rows(l) == rows(m), "dimensions mismatch l(%,%), m(%,%)\n", rows(l), cols(l), rows(m), cols(m));
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


//  ######  ##       ########  ######## ########
// ##    ## ##       ##     ##    ##         ##
// ##       ##       ##     ##    ##        ##
//  ######  ##       ########     ##       ##
//       ## ##       ##   ##      ##      ##
// ##    ## ##       ##    ##     ##     ##
//  ######  ######## ##     ##    ##    ########

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

        // ? this test should be disabled since `solve_linear_right_triangular` works with any matrix and prevents unnecessary allocations
        // assert(is_right_trapezoidal(R), "Matrix is not right trapezoidal: M(%) >= N(%)?\n", cols(R), rows(R));
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
decompose_LR :: inline (m: $A/MatrixType, L: *$B/MatrixType, R: *$C/MatrixType) {
    decompose_LR_default(m,L,R);
}

#scope_file
decompose_LR_default :: (m: $A/MatrixType, L: *$B/MatrixType, R: *$C/MatrixType) {
    #if CHECKS {
        assert(rows(m) == rows(L) && cols(m) == cols(L), "dimension mismatch: m(%,%), L(%,%)\n", rows(m), cols(m), rows(L), cols(L));
        assert(rows(m) == rows(R) && cols(m) == cols(R), "dimension mismatch: m(%,%), R(%,%)\n", rows(m), cols(m), rows(R), cols(R));
    }

    mat_init_unit(*L);

    for m {
        if it_column >= it_row {
            set(R, it_row, it_column, it);
        } else {
            set(L, it_row, it_column, it);
        }
    }
}
#scope_export




//  ######   ########         ##    ## ########
// ##    ##  ##               ###   ## ##     ##
// ##        ##               ####  ## ##     ##
// ##   #### ######           ## ## ## ########
// ##    ##  ##               ##  #### ##
// ##    ##  ##               ##   ### ##
//  ######   ##       ####### ##    ## ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.1 Gaussian Factorization - No Pivoting, page 199
gf_np :: gaussian_factorization_no_pivot;
gaussian_factorization_no_pivot :: inline (m: *$A/MatrixType) {
    gaussian_factorization_no_pivot_default(m);
}

#scope_file
gaussian_factorization_no_pivot_default :: (m: *$A/MatrixType) {

    r := rows(m);
    c := cols(m);

    for k: 0..min(r,c)-1 {
        assert(abs(get(m,k,k)) != zero(A.data_type), "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        for i: k+1..r-1 {
            set(m,i,k, get(m,i,k)/get(m,k,k));
        }
        for j: k+1..c-1 {
            for i: k+1..r-1 {
                set(m,i,j, get(m,i,j) - get(m,i,k) * get(m,k,j));
            }
        }
    }
}
#scope_export




//  ######  ##        ######   ########         ##    ## ########
// ##    ## ##       ##    ##  ##               ###   ## ##     ##
// ##       ##       ##        ##               ####  ## ##     ##
//  ######  ##       ##   #### ######           ## ## ## ########
//       ## ##       ##    ##  ##               ##  #### ##
// ##    ## ##       ##    ##  ##               ##   ### ##
//  ######  ########  ######   ##       ####### ##    ## ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.2 Solve Linear System via Gaussian Factorization - No Pivoting, page 199
slgf_np :: solve_linear_gaussian_factorization_no_pivot;
solve_LRnp :: solve_linear_gaussian_factorization_no_pivot;
solve_LUnp :: solve_linear_gaussian_factorization_no_pivot;
solve_linear_gaussian_factorization_no_pivot :: inline (a: *$A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    solve_linear_gaussian_factorization_no_pivot_default(a,b,x);
}

#scope_file
solve_linear_gaussian_factorization_no_pivot_default :: (a: *$A/MatrixType, b: $B/VectorType, x: *$C/VectorType) {
    gaussian_factorization_no_pivot(a);
    solve_linear_left_trapezoidal( <<a, b, x); // ? this should be possible, see function implementations!
    solve_linear_right_trapezoidal(<<a, <<x, x); // ? this should be possible, see function implementations!
}
#scope_export




//  ######   ########         ######## ########
// ##    ##  ##               ##       ##     ##
// ##        ##               ##       ##     ##
// ##   #### ######           ######   ########
// ##    ##  ##               ##       ##
// ##    ##  ##               ##       ##
//  ######   ##       ####### ##       ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.3 Gaussian Factorization with Full Pivoting, page 204
gf_fp :: gaussian_factorization_full_pivot;
gaussian_factorization_full_pivot :: inline (m: *$A/MatrixType, rows_perm: *$B/VectorType, cols_perm: *$C/VectorType) -> row_permutations: int, col_permutations: int {
    rp, cp := gaussian_factorization_full_pivot_default(m,rows_perm,cols_perm);
    return rp, cp;
}

#scope_file
gaussian_factorization_full_pivot_default :: (m: *$A/MatrixType, rows_perm: *$B/VectorType, cols_perm: *$C/VectorType) -> row_permutations: int, col_permutations: int {
// gaussian_factorization_full_pivot :: (m: $A/MatrixType)
//  -> A #must, 
//     rows_perm: [A.N]int #must, 
//     cols_perm: [A.M]int #must, 
//     row_permutations: int, 
//     col_permutations: int 
// {
    // println("gaussian_factorization_full_pivot");

    #if CHECKS {
        assert(rows(m) == dim(rows_perm), "dimensions mismatch m(%,%), rows(%)\n", rows(m), cols(m), dim(rows_perm));
        assert(cols(m) == dim(cols_perm), "dimensions mismatch m(%,%), cols(%)\n", rows(m), cols(m), dim(cols_perm));
    }

    mr := rows(m);
    mc := cols(m);

    row_permutations := 0;
    col_permutations := 0;

    for rows_perm {
        set(rows_perm, it_index, ncast(rows_perm.data_type, it_index));
    }
    for cols_perm {
        set(cols_perm, it_index, ncast(cols_perm.data_type, it_index));
    }

    // println("    entering for k");
    for k: 0..min(mr, mc)-1 {
        ii := k;
        jj := k;
        rho := abs(get(m,k,k));

        // println("    setting ii, jj");
        for j: k..mc-1 {
            for i: k..mr-1 {
                tmp := abs(get(m,i,j));
                if tmp > rho {
                    ii = i;
                    jj = j;
                    rho = tmp;
                }
            }
        }

        assert(rho != 0.0, "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        // instead of swapping, use the permutation arrays rows and cols to access the 
        // entries of m. Should be faster?
        // println("    permute rows");
        if ii > k {
            // println("        swap_rows");
            swap_rows(m, k, ii);
            // println("        swap row array entries");
            swap(rows_perm, k, ii);
            row_permutations += 1;
        }
        // println("    permute columns");
        if jj > k {
            swap_cols(m, k, jj);
            swap(cols_perm, k, jj);
            col_permutations += 1;
        }

        // if k == min(A.N, A.M)-1 then break;

        // println("    updating matrix");
        for i: k+1..mr-1 {
            set(m,i,k, get(m,i,k)/get(m,k,k));
        }
        // println("    updating matrix 2");
        for j: k+1..mc-1 {
            for i: k+1..mr-1 {
                set(m,i,j, get(m,i,j) - get(m,i,k) * get(m,k,j));
            }
        }
    }

    // println("    returning results");
    return row_permutations, col_permutations;
}
#scope_export



//  ######  ##        ######   ########         ######## ########
// ##    ## ##       ##    ##  ##               ##       ##     ##
// ##       ##       ##        ##               ##       ##     ##
//  ######  ##       ##   #### ######           ######   ########
//       ## ##       ##    ##  ##               ##       ##
// ##    ## ##       ##    ##  ##               ##       ##
//  ######  ########  ######   ##       ####### ##       ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.4 Solve Linear System via Gaussian Factorization with Full Pivoting, page 206
slgf_fp :: solve_linear_gaussian_factorization_full_pivot;
solve_LRfp :: solve_linear_gaussian_factorization_full_pivot;
solve_LUfp :: solve_linear_gaussian_factorization_full_pivot;

solve_linear_gaussian_factorization_full_pivot :: inline (
    a: *$A/MatrixType, 
    b: *$B/VectorType, 
    x: *$C/VectorType, 
    rows_perm: *$D/VectorType, 
    cols_perm: *$E/VectorType
) -> num_rows_perm: int, num_cols_perm: int {
    rp, cp := solve_linear_gaussian_factorization_full_pivot_default(a,b,x,rows_perm,cols_perm);
    return rp, cp;
}

#scope_file
solve_linear_gaussian_factorization_full_pivot_default :: (
    a: *$A/MatrixType, 
    b: *$B/VectorType,
    x: *$C/VectorType, 
    rows_perm: *$D/VectorType, 
    cols_perm: *$E/VectorType
) -> num_rows_perm: int, num_cols_perm: int {
    #if CHECKS {
        assert(rows(a) == dim(rows_perm));
        assert(cols(a) == dim(cols_perm));
    }

    num_rows_perm, num_cols_perm := gaussian_factorization_full_pivot(a, rows_perm, cols_perm);
    permute(b, <<rows_perm);
    solve_linear_left_trapezoidal( <<a, <<b, x);
    solve_linear_right_trapezoidal(<<a, <<x, x);
    permute(x, <<cols_perm);

    return num_rows_perm, num_cols_perm;
}
#scope_export






//  ######   ########         ########  ########
// ##    ##  ##               ##     ## ##     ##
// ##        ##               ##     ## ##     ##
// ##   #### ######           ########  ########
// ##    ##  ##               ##   ##   ##
// ##    ##  ##               ##    ##  ##
//  ######   ##       ####### ##     ## ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7. Gaussian Factorization with Partial Pivoting, page 210
gf_rp :: gaussian_factorization_row_pivot;
gaussian_factorization_row_pivot :: inline (m: *$A/MatrixType, row_perms: *$B/VectorType) -> permutations: int {
    return gaussian_factorization_row_pivot_default(m,row_perms);
}

#scope_file
gaussian_factorization_row_pivot_default :: (m: *$A/MatrixType, row_perms: *$B/VectorType) -> permutations: int {
    #if CHECKS {
        assert(rows(m) == dim(row_perms));
    }
    permutations := 0;

    for row_perms {
        set(row_perms, it_index, ncast(B.data_type, it_index));
    }

    mr := rows(m);
    mc := cols(m);

    for k: 0..min(mr,mc)-1 {
        ii := k;
        rho := abs(get(m,k,k));

        for i: k..mr-1 {
            tmp := abs(get(m,i,k));
            if tmp > rho {
                ii = i;
                rho = tmp;
            }
        }

        assert(rho != 0.0, "\nCannot perform Gaussain factorization due to 0 in diagonal.");

        if ii > k {
            swap_rows(m, k, ii);
            swap(row_perms, k, ii);
            permutations += 1;
        }

        for i: k+1..mr-1 {
            set(m,i,k, get(m,i,k)/get(m,k,k));
        }
        for j: k+1..mc-1 {
            for i: k+1..mr-1 {
                set(m,i,j, get(m,i,j) - get(m,i,k) * get(m,k,j));
            }
        }
    }

    return permutations;
}
#scope_export


//  ######  ##        ######   ########         ########  ########
// ##    ## ##       ##    ##  ##               ##     ## ##     ##
// ##       ##       ##        ##               ##     ## ##     ##
//  ######  ##       ##   #### ######           ########  ########
//       ## ##       ##    ##  ##               ##   ##   ##
// ##    ## ##       ##    ##  ##               ##    ##  ##
//  ######  ########  ######   ##       ####### ##     ## ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.7.7 Solve Linear System via Gaussian Factorization with Partial Pivoting, page 215
slgf_rp :: solve_linear_gaussian_factorization_row_pivot;
solve_LRrp :: solve_linear_gaussian_factorization_row_pivot;
solve_LUrp :: solve_linear_gaussian_factorization_row_pivot;
solve_linear_gaussian_factorization_row_pivot :: inline (
    a: *$A/MatrixType, 
    b: *$B/VectorType, 
    x: *$C/VectorType, 
    rows_perm: *$D/VectorType
) -> num_rows_perm: int {
    return solve_linear_gaussian_factorization_row_pivot_default(a,b,x,rows_perm);
}

#scope_file
solve_linear_gaussian_factorization_row_pivot_default :: (
    a: *$A/MatrixType, 
    b: *$B/VectorType, 
    x: *$C/VectorType, 
    rows_perm: *$D/VectorType
) -> num_rows_perm: int {
    #if CHECKS {
        assert(rows(a) == dim(rows_perm));
    }

    num_rows_perm := gaussian_factorization_row_pivot(a, rows_perm);
    permute(b, <<rows_perm);
    solve_linear_left_trapezoidal( <<a, <<b, x);
    solve_linear_right_trapezoidal(<<a, <<x, x);

    return num_rows_perm;
}
#scope_export




// #### ##    ## ##     ## ######## ########   ######  ########
//  ##  ###   ## ##     ## ##       ##     ## ##    ## ##
//  ##  ####  ## ##     ## ##       ##     ## ##       ##
//  ##  ## ## ## ##     ## ######   ########   ######  ######
//  ##  ##  ####  ##   ##  ##       ##   ##         ## ##
//  ##  ##   ###   ## ##   ##       ##    ##  ##    ## ##
// #### ##    ##    ###    ######## ##     ##  ######  ########

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? Algorithm 3.10.1 Matrix Inverse, page 261
inv :: inverse;
inverse :: inline (a: $A/MatrixType, ainv: *A) {
    inverse_default(a,ainv);
}

#scope_file
inverse_default :: (a: $A/MatrixType, ainv: *A) {
    row_perms := dhvec(u32, rows(a), false); // ! alloc, u32 limitation?
    defer free(row_perms);

    gaussian_factorization_row_pivot(*a, *row_perms);

    mat_init_unit(ainv);

    solve_linear_left_trapezoidal( a, <<ainv, ainv);
    solve_linear_right_trapezoidal(a, <<ainv, ainv);

    permute_cols(ainv, row_perms);
}
#scope_export



// ########  ######## ######## ######## ########  ##     ## #### ##    ##    ###    ##    ## ########
// ##     ## ##          ##    ##       ##     ## ###   ###  ##  ###   ##   ## ##   ###   ##    ##
// ##     ## ##          ##    ##       ##     ## #### ####  ##  ####  ##  ##   ##  ####  ##    ##
// ##     ## ######      ##    ######   ########  ## ### ##  ##  ## ## ## ##     ## ## ## ##    ##
// ##     ## ##          ##    ##       ##   ##   ##     ##  ##  ##  #### ######### ##  ####    ##
// ##     ## ##          ##    ##       ##    ##  ##     ##  ##  ##   ### ##     ## ##   ###    ##
// ########  ########    ##    ######## ##     ## ##     ## #### ##    ## ##     ## ##    ##    ##

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 3.11 Determinants, page 263
// ?
// ? it is not recommended to calculate the determinant of large matrices
det :: determinant;
determinant :: inline (a: *$A/MatrixType) -> A.data_type {
    return determinant_default(a);
}

#scope_file
determinant_default :: (a: *$A/MatrixType) -> A.data_type {
    #if CHECKS {
        assert(rows(a) == cols(a), "matrix not quadratic!\n (%,%)", rows(a), cols(a));
    }
    // TODO implement optimized formula for 3x3 matrix
    if rows(a) == 2 {
        return get(a,0,0)*get(a,1,1) - get(a,0,1)*get(a,1,0);
    } else {    
        row_perms := dhvec(A.data_type, rows(a), false); // ! alloc
        defer free(row_perms);

        permutations := gaussian_factorization_row_pivot(a, *row_perms);
    
        p := one(A.data_type);
        for 0..rows(a)-1 {
            p *= get(a,it,it);
        }

        // +1 for even, -1 for odd num of permutations
        s : float64 = 2.0 * cast(float64)(permutations % 2 == 0) - 1.0; // ! this will give the wrong return type at some point
    
        return s * p;
    }
}
#scope_export































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

        x : DenseVector(C64, 2);

        solve_linear_2x2(A, b, *x);
        println("result: %", pstr(x));

        sol := dvec(2,  // ? this is the correct answer
            -0.25 + 1.75*I,
            0.25 - 0.75*I
        );
        assert(x == sol);
    }

    {
        println_push("test 2", color = .FG_GREEN);

        L := dmat(C64, 3, 3,
            R, C_ZERO, C_ZERO,
            3.0*R, R, C_ZERO,
            I, R+I, R
        );
        b := dvec(C64, 3, 3.0*R, -I, R-2.0*I);
        res : DenseVector(C64, 3) = ---;

        solve_linear_left_triangular(L, b, *res);

        sol : DenseVector(C64, 3);
        mat_mul_vec(L, res, *sol);

        println("% = % ?", pstr(sol), pstr(b));
        assert(sol == b);
    }

    {
        println_push("test 3", color = .FG_GREEN);

        M := dmat(3,3,
            2.0*R, I, R-2.0*I,
            C_ZERO, R, 3.0*I,
            C_ZERO, C_ZERO, -3.0*R
        );
        b := dvec(3, 3.0*R, -I, R-2.0*I);
        res : DenseVector(C64, 3) = ---;
        solve_linear_right_triangular(M, b, *res);
        println("res %", pstr(res));

        sol : DenseVector(C64, 3);
        tmp : DenseVector(C64, 3) = ---;
        mat_mul_vec(M, res, *sol);
        copy(*tmp, sol);
        sub(*sol, b);
        sol_norm := norm_2(sol);
        println("%", pstr(sol));
        println("sol_norm = %", sol_norm); 
        assert(sol_norm < 1.0e-8, "res = %, M * res = % =? b = %", res, tmp, b);
    }

    {
        println_push("test 4", color = .FG_GREEN);

        m := dmat(int, 4, 4,
            1, 3, 4, 5,
            0, 2, 3, 4,
            0, 1, 2, 3,
            0, 0, 0, 2
        );
        cm : DenseMatrix(C64, 4, 4) = ---;
        cast_to(*cm, m);
        
        b := dvec(float64, 4, 2,3,4,5);
        cb : DenseVector(C64, 4) = ---;
        cast_to(*cb, b);

        println("cm %", pstr(cm));
        println("cb %", pstr(cb));

        res : DenseVector(C64, 4) = ---;
        solve_linear_right_quasi_triangular(cm,cb, *res);
        println("res = %", pstr(res));

        // * check result
        sol : DenseVector(C64, 4);
        mat_mul_vec(cm, res, *sol);
        sub(*sol, cb);
        assert(norm_2(sol) < 1.0e-8);
    }

    {
        println_push("test 5", color = .FG_GREEN);

        M := dmat(C64, 3, 2, 
            R, C_ZERO, 
            -R, R, 
            -R, -R
        );
        b := dvec(C64, 3, R, 2.0*R, -4.0*R);
        res : DenseVector(C64, 2);

        println("calling solve_linear_left_trapezoidal...");
        solve_linear_left_trapezoidal(M, b, *res);
        println("sllt: %", pstr(res));

        // * check solution
        sol : DenseVector(C64, 3);
        mat_mul_vec(M, res, *sol);
        println("sol = %", sol);
        println("b = %", b);
        assert( approx(sol, b) );
    }

    {
        println_push("test 6", color = .FG_GREEN);

        M := dmat(C64, 2, 3, 
            R, R, C_ZERO, 
            C_ZERO, R, C_ZERO 
        );
        b := dvec(C64, 2, 4*R, 6*R);
        res : DenseVector(C64, 3);

        solve_linear_right_trapezoidal(M, b, *res);
        
        // * check solution
        sol : DenseVector(C64, 2) = ---;
        // copy(*sol, res);
        mat_mul_vec(M, res, *sol);
        
        println("slrt: %", pstr(res));
        assert( approx(sol, b) );
    }
    {
        println_push("test 7", color = .FG_GREEN);

        M := dmat(C64, 2, 2,
            2*R, 3*R,
            -2*R, R
        );
        println("M = %", pstr(M));
        res := dmat(C64, 2, 2,
            2*R, 3*R,
            -1*R, 4*R
        );
        gaussian_factorization_no_pivot(*M);
        println("gaussian factorization: %", pstr(M));
        assert( M == res , "%\n != %\n!", str(M), str(res));
    }

    {
        println_push("test 8", color = .FG_GREEN);

        M := dmat(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := dvec(2, 4*R, 6*I);
        x : DenseVector(C64, 2) = ---;

        slgf_np(*M, b, *x);
        println("solve linear system via Gaussian factorization (no pivoting): %", pstr(x));
        assert(x == dvec(2, 2.5-4.5*I, 4.5-1.5*I), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        println_push("test 9", color = .FG_GREEN);

        M := dmat(C64, 3, 3,
            C_ZERO, 4*R, -3*R,
            2*R, C_ZERO, 4*R,
            3*R, -5*R, 3*R
        );
        println("M %", pstr(M));

        rows : DenseVector(u8, 3) = ---;
        cols : DenseVector(u8, 3) = ---;

        rps, cps := gaussian_factorization_full_pivot(*M, *rows, *cols);
        println("M = %", pstr(M));
        println("rows = %", pstr(rows));
        println("cols = %", pstr(cols));
        println("rps %, cps %", rps, cps);
    }

    {
        println_push("test 10", color = .FG_GREEN);

        M := dmat(C64, 3, 3,
            C_ZERO, 4*R, -3*R,
            2*R, C_ZERO, 4*R,
            3*R, -5*R, 3*R
        );
        rows : DenseVector(u8, 3) = ---;
        cols : DenseVector(u8, 3) = ---;

        rp, cp := gaussian_factorization_full_pivot(*M, *rows, *cols);

        println("M = %", pstr(M));
    }

    {
        println_push("test 11", color = .FG_GREEN);

        M := dmat(2, 2,
            R, I,
            -3*R, 2.0 - I
        );
        b := dvec(2, 4*R, 6*I);
        res : DenseVector(C64, 2) = ---;

        rperm : DenseVector(u8, 2) = ---;
        cperm : DenseVector(u8, 2) = ---;

        println("solving lin sys via gf_fp...");
        rps, cps := slgf_fp(*M, *b, *res, *rperm, *cperm);
        println("solution: %", pstr(res));
        sol := dvec(2, 2.5-4.5*I, 4.5-1.5*I);
        assert(approx(res, sol), "\nCould not solve linear system via Gaussain factorization!");
    }

    {
        println_push("test 12", color = .FG_GREEN);

        M := dmat(2,2,
            R, I,
            -3*R, 2.0 - I
        );
        b := dvec(2,   4*R, 6*I);
        res : DenseVector(C64, 2);
        rowp : DenseVector(u8, 2) = ---;

        nrowp := slgf_rp(*M, *b, *res, *rowp);
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

        res : DenseMatrix(C64, 3, 3) = ---;

        solve_linear_left_trapezoidal(L,one, *res);
        
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

        res : DenseMatrix(C64, 3, 3) = ---;

        solve_linear_right_trapezoidal(M,one,*res);
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
        iM : DenseMatrix(C64, 3, 3) = ---;
        inverse(M, *iM);
        
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
        d := determinant(*M);
        println("det(M) = %", str(d));
        assert(d == -4*I);
    }
}