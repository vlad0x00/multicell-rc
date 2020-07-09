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

const REAL CELL_RADIUS = 2.0;

/* IF_GRID_SPACING is the unit length of each voxel in the Simulation Domain
 The Simulation Domain size is set in the model XML file
 The Grid spacing can not be less than maximum cell agent diameter */
const REAL IF_GRID_SPACING = 10.0; // Micrometers

/* A baseline time step is the largest time step used in Biocellion
 Users can split a baseline time step into one or more state-and-grid time steps */
const REAL BASELINE_TIME_STEP_DURATION = 1.0;
const S32 NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE = 1;

const S32 SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL = 2;
const S32 NUM_AMR_LEVELS = 1;
const S32 NUM_PDE_TIME_STEPS_PER_STATE_AND_GRID_TIME_STEP = 1;
const REAL GRID_PHI_NORM_THRESHOLD = 1e-10;

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
	S64 alpha;
	S64 beta;
	S64 numCytokines;
	S64 secretionLow;
	S64 secretionHigh;
	S64 cytokineThreshold;
} GlobalDataFormat;

static inline GlobalDataFormat getGlobalDataFormat() {
	const auto& g = Info::getGlobalDataRef();
	GlobalDataFormat format = *((GlobalDataFormat*)(&(g[0])));
	return format;
}

static inline S32 getNumGenes() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((S32*)(&(g[f.numGenes])));
}

static inline S32 getNumCells() { 
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();	
	return *((S32*)(&(g[f.numCells])));
}

static inline S32 getNumCellTypes() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((S32*)(&(g[f.numCellTypes])));
}

static inline S32* getNv(S32 cellType) {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return (S32*)(&(g[f.nv])) + cellType * getNumGenes();
}

static inline S32* getVarf(S32 cellType, S32 gene) {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return (S32*)(&(g[f.varf])) + ((S32*)(&(g[f.varfOffsets])))[cellType * getNumGenes() + gene];
}

static inline S32* getTt(S32 cellType, S32 gene) {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();	
	return (S32*)(&(g[f.tt])) + ((S32*)(&(g[f.ttOffsets])))[cellType * getNumGenes() + gene];
}

static inline S32* getGeneInitialStates(S32 cellType) {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();	
	return (S32*)(&(g[f.geneInitialStates + cellType * getNumGenes()]));
}

static inline S32* getInputSignal() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();	
	return (S32*)(&(g[f.inputSignal]));
}

static inline REAL getAlpha() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((REAL*)(&(g[f.alpha])));
}

static inline REAL getBeta() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((REAL*)(&(g[f.beta])));
}

static inline S32 getNumCytokines() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((S32*)(&(g[f.numCytokines])));
}

static inline REAL getSecretionLow() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((REAL*)(&(g[f.secretionLow])));
}

static inline REAL getSecretionHigh() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((REAL*)(&(g[f.secretionHigh])));
}

static inline REAL getCytokineThreshold() {
	const auto& g = Info::getGlobalDataRef();
	const auto f = getGlobalDataFormat();
	return *((REAL*)(&(g[f.cytokineThreshold])));
}

/* MODEL END */

#endif/* #ifndef __MY_DEFINE_H__ */

