// ? LAPACK *larrv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarrv/dlarrv/clarrv/zlarrv
// *LARRV computes eigenvectors of a tridiagonal matrix from its representation tree.

slarrv :: #bake_arguments _larrv(RT=float32, MT=float32);
dlarrv :: #bake_arguments _larrv(RT=float64, MT=float64);
clarrv :: #bake_arguments _larrv(RT=float32, MT=C32);
zlarrv :: #bake_arguments _larrv(RT=float64, MT=C64);

_larrv :: (
    $$n: int, $$vl: $RT, $$vu: RT, d: $D/ArrayType1D, l: $L/ArrayType1D, $$pivmin: RT,
    isplit: $ISPLIT/ArrayType1D, $$m: int, $$dol: int, $$dou: int, $$minrgp: RT,
    $$rtol1: RT, $$rtol2: RT, w: $W/ArrayType1D, werr: $WERR/ArrayType1D, wgap: $WGAP/ArrayType1D,
    iblock: $IBLOCK/ArrayType1D, indexw: $INDEXW/ArrayType1D, gers: $GERS/ArrayType2D, z: *$Z/ArrayType2D, $$ldz: int,
    isuppz: *$ISUPPZ/ArrayType1D, work: *$WORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           D.data_type == RT && L.data_type == RT && ISPLIT.data_type == int &&
           W.data_type == RT && WERR.data_type == RT && WGAP.data_type == RT &&
           IBLOCK.data_type == int && INDEXW.data_type == int && GERS.data_type == RT &&
           Z.data_type == RT && ISUPPZ.data_type == int && WORK.data_type == RT && IWORK.data_type == int;
} {
    info.* = 0;
    _ = n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, MT;
    // Staged port note: the full *larrv body is still pending.
}
