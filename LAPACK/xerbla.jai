


xerbla :: (msg: string, info: Any) {
    assert(false, sprint("LAPACK ERROR: %, %\n", msg, info));
}