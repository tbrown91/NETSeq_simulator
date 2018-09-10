#!/bin/bash

./test_parameterSet WT_NETseqDoubleStallParams.txt WT_NETseqDoubleStallParams &
./test_parameterSet WT_NETseqHalfStallParams.txt WT_NETseqHalfStallParams &
./test_parameterSet WT_NETseqDoubleBacktrackParams.txt WT_NETseqDoubleBacktrackParams &
./test_parameterSet WT_NETseqHalfBacktrackParams.txt WT_NETseqHalfBacktrackParams &

./test_parameterSet WT_NETseqTripleStallParams.txt WT_NETseqTripleStallParams &
./test_parameterSet WT_NETseqThirdStallParams.txt WT_NETseqThirdStallParams &
./test_parameterSet WT_NETseqTripleBacktrackParams.txt WT_NETseqTripleBacktrackParams &
./test_parameterSet WT_NETseqThirdBacktrackParams.txt WT_NETseqThirdBacktrackParams &

./test_parameterSet WT_NETseqQuadrupleStallParams.txt WT_NETseqQuadrupleStallParams &
./test_parameterSet WT_NETseqQuarterStallParams.txt WT_NETseqQuarterStallParams &
./test_parameterSet WT_NETseqQuadrupleBacktrackParams.txt WT_NETseqQuadrupleBacktrackParams &
./test_parameterSet WT_NETseqQuarterBacktrackParams.txt WT_NETseqQuarterBacktrackParams &

#./test_parameterSet WT_NETseqNoBacktrackingParams.txt WT_NETseqNoBacktrackingParams &
#./test_parameterSet WT_NETseqNoBacktrackReleaseParams.txt WT_NETseqNoBacktrackReleaseParams &
