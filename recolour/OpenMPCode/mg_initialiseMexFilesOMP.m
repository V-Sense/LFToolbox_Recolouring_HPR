mex -g 'recolour/OpenMPCode/mex_mgRecolourParallel_1.cpp' COMPFLAGS="/openmp $COMPFLAGS" -lgomp
mex -g 'recolour/OpenMPCode/mex_mgRecolourParallel_Mask.cpp' COMPFLAGS="/openmp $COMPFLAGS" -lgomp
mex -g 'recolour/OpenMPCode/mex_mgRecolourParallelTPS.cpp' COMPFLAGS="/openmp $COMPFLAGS" -lgomp