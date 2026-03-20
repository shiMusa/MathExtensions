// ? LAPACK *hb2st_kernels()
// ? http://www.netlib.org/lapack/explore-html/
//
// ? original comment chb2st_kernels/zhb2st_kernels
//
// *HB2ST_KERNELS is an internal routine used by *HETRD_HB2ST.
//
chb2st_kernels :: #bake_arguments _hb2st_kernels(RT=float32, MT=C32);
zhb2st_kernels :: #bake_arguments _hb2st_kernels(RT=float64, MT=C64);

#scope_file
_hb2st_kernels :: (
    $$uplo: u8,
    $$wantz: bool,
    $$ttype: int,
    $$st: int,
    $$ed: int,
    $$sweep: int,
    $$n: int,
    $$nb: int,
    $$ib: int,
    a: *$A/ArrayType1D,
    $$lda: int,
    v: *$V/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    $$ldvt: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           V.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    // Local aliases/values.
    zero_t :: #run zero(MT);
    one_t  :: #run one(MT);

    upper : bool = lsame(uplo, #char "U");

    // Index layout in the workspace buffer (as per Netlib).
    ajeter : int = ib + ldvt; // Not used in this reduced port, kept for parity.

    dpos, ofdpos : int = ---;
    if upper {
        dpos   = 2*nb + 1;
        ofdpos = 2*nb;
    } else {
        dpos   = 1;
        ofdpos = 2;
    }

    // Treat the A workspace as a 2D matrix view with leading dimension LDA.
    a_mat := array2d_view(a, lda, row_major=false);

    // For this routine, WANTZ does not change the computed VPOS/TAUPOS.
    vpos   : int = mod(sweep - 1, 2) * n + st;
    taupos : int = mod(sweep - 1, 2) * n + st;

    if upper {
        if ttype == 1 {
            lm : int = ed - st + 1;

            fset(v, vpos, one_t);
            for i:1..lm-1 {
                z := conjugate(fget(a_mat, ofdpos - i, st + i));
                fset(v, vpos + i, z);
                fset(a_mat, ofdpos - i, st + i, zero_t);
            }

            // CLARFG needs ALPHA by reference; we store the conjugated alpha into A first.
            alpha_ptr := a + (st - 1)*lda + (ofdpos - 1);
            fset(a_mat, ofdpos, st, conjugate(fget(a_mat, ofdpos, st)));

            #if MT == C32 {
                clarfg(lm, alpha_ptr, v + vpos, 1, tau + (taupos - 1));
            } else {
                zlarfg(lm, alpha_ptr, v + vpos, 1, tau + (taupos - 1));
            }

            tau_conj := conjugate(fget(tau, taupos));
            c_raw := a + (st - 1)*lda + (dpos - 1);
            c_view := array2d_view(c_raw, lda - 1, row_major=false);

            // Apply CLARFY to the Hermitian block.
            #if MT == C32 {
                clarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            } else {
                zlarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            }
        }

        if ttype == 3 {
            lm : int = ed - st + 1;

            tau_conj := conjugate(fget(tau, taupos));
            c_raw := a + (st - 1)*lda + (dpos - 1);
            c_view := array2d_view(c_raw, lda - 1, row_major=false);

            #if MT == C32 {
                clarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            } else {
                zlarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            }
        }

        if ttype == 2 {
            j1 : int = ed + 1;
            j2 : int = min(ed + nb, n);
            ln : int = ed - st + 1;
            lm : int = j2 - j1 + 1;

            if lm > 0 {
                // Apply reflector from the left: C := H * C
                c1_raw := a + (j1 - 1)*lda + ((dpos - nb) - 1);
                c1_view := array2d_view(c1_raw, lda - 1, row_major=false);
                tau_conj := conjugate(fget(tau, taupos));

                #if MT == C32 {
                    clarfx(#char "L", ln, lm, v + (vpos - 1), tau_conj, *c1_view, lda - 1, work);
                } else {
                    zlarfx(#char "L", ln, lm, v + (vpos - 1), tau_conj, *c1_view, lda - 1, work);
                }

                // Update V/TAU positions for the newly created bulge.
                vpos   = mod(sweep - 1, 2) * n + j1;
                taupos = mod(sweep - 1, 2) * n + j1;

                fset(v, vpos, one_t);
                for i:1..lm-1 {
                    z := conjugate(fget(a_mat, dpos - nb - i, j1 + i));
                    fset(v, vpos + i, z);
                    fset(a_mat, dpos - nb - i, j1 + i, zero_t);
                }

                alpha_ptr2 := a + (j1 - 1)*lda + ((dpos - nb) - 1);
                fset(a_mat, dpos - nb, j1, conjugate(fget(a_mat, dpos - nb, j1)));

                #if MT == C32 {
                    clarfg(lm, alpha_ptr2, v + vpos, 1, tau + (taupos - 1));
                } else {
                    zlarfg(lm, alpha_ptr2, v + vpos, 1, tau + (taupos - 1));
                }

                // Apply reflector from the right.
                c2_raw := a + (j1 - 1)*lda + ((dpos - nb + 1) - 1);
                c2_view := array2d_view(c2_raw, lda - 1, row_major=false);

                tau_val := fget(tau, taupos);
                #if MT == C32 {
                    clarfx(#char "R", ln - 1, lm, v + (vpos - 1), tau_val, *c2_view, lda - 1, work);
                } else {
                    zlarfx(#char "R", ln - 1, lm, v + (vpos - 1), tau_val, *c2_view, lda - 1, work);
                }
            }
        }
    } else {
        // Lower triangle case.
        if ttype == 1 {
            lm : int = ed - st + 1;

            fset(v, vpos, one_t);
            for i:1..lm-1 {
                z := fget(a_mat, ofdpos + i, st - 1);
                fset(v, vpos + i, z);
                fset(a_mat, ofdpos + i, st - 1, zero_t);
            }

            alpha_ptr := a + ((st - 1) - 1)*lda + (ofdpos - 1);
            #if MT == C32 {
                clarfg(lm, alpha_ptr, v + vpos, 1, tau + (taupos - 1));
            } else {
                zlarfg(lm, alpha_ptr, v + vpos, 1, tau + (taupos - 1));
            }

            tau_conj := conjugate(fget(tau, taupos));
            c_raw := a + (st - 1)*lda + (dpos - 1);
            c_view := array2d_view(c_raw, lda - 1, row_major=false);

            #if MT == C32 {
                clarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            } else {
                zlarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            }
        }

        if ttype == 3 {
            lm : int = ed - st + 1;
            tau_conj := conjugate(fget(tau, taupos));
            c_raw := a + (st - 1)*lda + (dpos - 1);
            c_view := array2d_view(c_raw, lda - 1, row_major=false);

            #if MT == C32 {
                clarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            } else {
                zlarfy(uplo, lm, v + (vpos - 1), 1, tau_conj, *c_view, lda - 1, work);
            }
        }

        if ttype == 2 {
            j1 : int = ed + 1;
            j2 : int = min(ed + nb, n);
            ln : int = ed - st + 1;
            lm : int = j2 - j1 + 1;

            if lm > 0 {
                // Apply reflector from the right.
                c1_raw := a + (st - 1)*lda + ((dpos + nb) - 1);
                c1_view := array2d_view(c1_raw, lda - 1, row_major=false);
                tau_val := fget(tau, taupos);

                #if MT == C32 {
                    clarfx(#char "R", lm, ln, v + (vpos - 1), tau_val, *c1_view, lda - 1, work);
                } else {
                    zlarfx(#char "R", lm, ln, v + (vpos - 1), tau_val, *c1_view, lda - 1, work);
                }

                // Update V/TAU positions for the newly created bulge.
                vpos   = mod(sweep - 1, 2) * n + j1;
                taupos = mod(sweep - 1, 2) * n + j1;

                fset(v, vpos, one_t);
                for i:1..lm-1 {
                    z := fget(a_mat, dpos + nb + i, st);
                    fset(v, vpos + i, z);
                    fset(a_mat, dpos + nb + i, st, zero_t);
                }

                alpha_ptr2 := a + (st - 1)*lda + ((dpos + nb) - 1);
                #if MT == C32 {
                    clarfg(lm, alpha_ptr2, v + vpos, 1, tau + (taupos - 1));
                } else {
                    zlarfg(lm, alpha_ptr2, v + vpos, 1, tau + (taupos - 1));
                }

                // Apply reflector from the left: C := H * C
                c2_raw := a + (st)*lda + ((dpos + nb - 1) - 1);
                c2_view := array2d_view(c2_raw, lda - 1, row_major=false);
                tau_conj := conjugate(fget(tau, taupos));

                #if MT == C32 {
                    clarfx(#char "L", lm, ln - 1, v + (vpos - 1), tau_conj, *c2_view, lda - 1, work);
                } else {
                    zlarfx(#char "L", lm, ln - 1, v + (vpos - 1), tau_conj, *c2_view, lda - 1, work);
                }
            }
        }
    }

    // Keep ldvt and wantz parameters "used" (for warnings / style parity).
    _ = ajeter;
    _ = wantz;
    _ = ldvt;
}

