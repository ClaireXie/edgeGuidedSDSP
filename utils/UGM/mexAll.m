
mex -outdir minFunc_2009 minFunc_2009/lbfgsC.c
mex -outdir minFunc_2009 minFunc_2009/mcholC.c

mex -outdir KPM KPM/max_mult.c

mex -Imex -outdir compiled mex/UGM_makeEdgeVEC.c
mex -Imex -outdir compiled mex/UGM_Decode_ExactC.c
mex -Imex -outdir compiled mex/UGM_Infer_ExactC.c
mex -Imex -outdir compiled mex/UGM_Infer_ChainC.c
mex -Imex -outdir compiled mex/UGM_makeClampedPotentialsC.c
mex -Imex -outdir compiled mex/UGM_Decode_ICMC.c
mex -Imex -outdir compiled mex/UGM_Decode_GraphCutC.c
mex -Imex -outdir compiled mex/UGM_Sample_GibbsC.c
mex -Imex -outdir compiled mex/UGM_Infer_MFC.c
mex -Imex -outdir compiled mex/UGM_Infer_LBPC.c
mex -Imex -outdir compiled mex/UGM_Decode_LBPC.c
mex -Imex -outdir compiled mex/UGM_Infer_TRBPC.c
mex -Imex -outdir compiled mex/UGM_Decode_TRBPC.c
mex -Imex -outdir compiled mex/UGM_CRF_makePotentialsC.c
mex -Imex -outdir compiled mex/UGM_CRF_PseudoNLLC.c
mex -Imex -outdir compiled mex/UGM_LogConfigurationPotentialC.c
mex -Imex -outdir compiled mex/UGM_Decode_AlphaExpansionC.c
mex -Imex -outdir compiled mex/UGM_Decode_AlphaExpansionBetaShrinkC.c
mex -Imex -outdir compiled mex/UGM_CRF_NLLC.c
