// ? LAPACK *hla_transtype()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chla_transtype
// *HLA_TRANSTYPE translates a BLAST transpose constant to a character code.

chla_transtype :: ($$trans: int) -> u8 {
    BLAS_NO_TRANS   :: 111;
    BLAS_TRANS      :: 112;
    BLAS_CONJ_TRANS :: 113;

    if trans == BLAS_NO_TRANS {
        return #char "N";
    } else if trans == BLAS_TRANS {
        return #char "T";
    } else if trans == BLAS_CONJ_TRANS {
        return #char "C";
    } else {
        return #char "X";
    }
}
