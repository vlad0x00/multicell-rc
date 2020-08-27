/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

#ifndef __MY_DEFINE_H__
#define __MY_DEFINE_H__

#include "biocellion.h"

/* define constants to be used inside model functions */

#define SYSTEM_DIMENSION 3

/* MODEL START */

/* A Uniform Random Number Generator. One should not use c++ function, rand(), as it is not thread-safe
 Instead, users should use Biocellion's built-in RNG */
typedef enum _model_rng_type_e {
  MODEL_RNG_UNIFORM,
  MODEL_RNG_GAUSSIAN,
  NUM_MODEL_RNGS
} model_rng_type_e;

const BOOL ENABLE_SUMMARY = false;

/* A baseline time step is the largest time step used in Biocellion
 Users can split a baseline time step into one or more state-and-grid time steps */
const REAL BASELINE_TIME_STEP_DURATION = 1.0;
const S32 NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE = 1;

const S32 SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL = 2;
const S32 NUM_AMR_LEVELS = 1;
const REAL GRID_PHI_NORM_THRESHOLD = 1e-10;

const REAL SIDE_BOUNDARY_INPUT_WALL_FRACTION = 0.3;

extern S32 gNumGenes;
extern S32 gNumCells;
extern S32 gNumCellTypes;
extern S32* gNv;
extern S32* gVarfOffsets;
extern S32* gTtOffsets;
extern S32* gVarf;
extern S32* gTt;
extern S32* gGeneInitialStates;
extern S32* gInputSignal;
extern REAL gAlphaInput;
extern REAL gBetaInput;
extern REAL gAlphaCytokines;
extern REAL gBetaCytokines;
extern S32 gNumCytokines;
extern REAL gSecretionLow;
extern REAL gSecretionHigh;
extern REAL gCytokineThreshold;
extern REAL gCellRadius;
extern REAL gIfGridSpacing;
extern REAL gDirichletBoundary;
extern REAL gInputThreshold;
extern S32 gZLayers;
extern S32 gYLayers;
extern S32 gXLayers;

typedef struct {
  S64 numGenes;
  S64 numCells;
  S64 numCellTypes;
  S64 nv;
	S64 varfOffsets;
	S64 ttOffsets;
  S64 varf;
  S64 tt;
  S64 geneInitialStates;
  S64 inputSignal;
	S64 alphaInput;
	S64 betaInput;
	S64 alphaCytokines;
	S64 betaCytokines;
	S64 numCytokines;
	S64 secretionLow;
	S64 secretionHigh;
	S64 cytokineThreshold;
	S64 cellRadius;
	S64 ifGridSpacing;
	S64 dirichletBoundary;
	S64 inputThreshold;
	S64 zLayers;
	S64 yLayers;
	S64 xLayers;
} GlobalDataFormat;

static inline S32* getNv(S32 cellType) {
	return gNv + cellType * gNumGenes;
}

static inline S32* getVarf(S32 cellType, S32 gene) {
	return gVarf + gVarfOffsets[cellType * gNumGenes + gene];
}

static inline S32* getTt(S32 cellType, S32 gene) {
	return gTt + gTtOffsets[cellType * gNumGenes + gene];
}

static inline S32* getGeneInitialStates(S32 cell) {
	return gGeneInitialStates + cell * gNumGenes;
}

/* MODEL END */

#endif/* #ifndef __MY_DEFINE_H__ */

