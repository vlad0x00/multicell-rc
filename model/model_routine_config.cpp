/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

#include <fstream>
#include <sstream>
#include <cmath>

S32 gNumGenes;
S32 gNumCells;
S32 gNumCellTypes;
S32* gNv;
S32* gVarfOffsets;
S32* gTtOffsets;
S32* gVarf;
S32* gTt;
S32* gGeneInitialStates;
S32* gInputSignal;
REAL gAlphaInput;
REAL gBetaInput;
REAL gAlphaCytokines;
REAL gBetaCytokines;
S32 gNumCytokines;
REAL gSecretionLow;
REAL gSecretionHigh;
REAL gCytokineThreshold;
REAL gCellRadius;
REAL gIfGridSpacing;
REAL gDirichletBoundary;
REAL gInputThreshold;
S32 gZLayers;
S32 gYLayers;
S32 gXLayers;
BOOL gKappa;

// Parameters provided to Biocellion and their indices in the input string array
const S32 NUM_GENES_PARAM = 0;
const S32 NUM_CELLS_PARAM = 1;
const S32 NUM_CELL_TYPES_PARAM = 2;
const S32 NV_FILE_PARAM = 3;
const S32 VARF_FILE_PARAM = 4;
const S32 TT_FILE_PARAM = 5;
const S32 GENE_INITIAL_STATES_FILE_PARAM = 6;
const S32 INPUT_SIGNAL_FILE_PARAM = 7;
const S32 ALPHA_CYTOKINES_PARAM = 8;
const S32 BETA_CYTOKINES_PARAM = 9;
const S32 NUM_CYTOKINES_PARAM = 10;
const S32 SECRETION_LOW_PARAM = 11;
const S32 SECRETION_HIGH_PARAM = 12;
const S32 CYTOKINE_THRESHOLD_PARAM = 13;
const S32 CELL_RADIUS_PARAM = 14;
const S32 IF_GRID_SPACING_PARAM = 15;
const S32 DIRICHLET_BOUNDARY_PARAM = 16;
const S32 INPUT_THRESHOLD_PARAM = 17;
const S32 ALPHA_INPUT_PARAM = 18;
const S32 BETA_INPUT_PARAM = 19;
const S32 Z_LAYERS_PARAM = 20;
const S32 Y_LAYERS_PARAM = 21;
const S32 X_LAYERS_PARAM = 22;
const S32 KAPPA_PARAM = 23;

// Returns a vector of strings of passed arguments
static Vector<std::string> readXMLParameters() {
  Vector<string> params;

  std::string paramsString = Info::getModelParam();
  std::string delimiter = " ";

  size_t pos = 0;
  std::string token;
  while ((pos = paramsString.find(delimiter)) != std::string::npos) {
    token = paramsString.substr(0, pos);
    if (token == "") { continue; }
    params.push_back(token);
    paramsString.erase(0, pos + delimiter.length());
  }
  if (!paramsString.empty()) {
    params.push_back(paramsString);
  }

  return params;
}

void ModelRoutine::updateIfGridSpacing( REAL& ifGridSpacing ) {
  /* MODEL START */

  const auto params = readXMLParameters();
  ifGridSpacing = std::stod(params[IF_GRID_SPACING_PARAM]);

  /* MODEL END */

  return;
}

void ModelRoutine::updateOptModelRoutineCallInfo( OptModelRoutineCallInfo& callInfo ) {
  /* MODEL START */

  callInfo.numComputeMechIntrctIters = 0;
  callInfo.numUpdateIfGridVarPreStateAndGridStepIters = 1; // Single ifGridUpdate call
  callInfo.numUpdateIfGridVarPostStateAndGridStepIters = 0;

  /* MODEL END */

  return;
}

void ModelRoutine::updateDomainBdryType( domain_bdry_type_e a_domainBdryType[DIMENSION] ) {
  /* MODEL START */

  /* We set the Boundary Types of the Simulation Domain
     There are only two Boundary Types: DOMAIN_BDRY_TYPE_PERIODIC, DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL */
  a_domainBdryType[0] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL; //+-x direction walls
  a_domainBdryType[1] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL; //+-y direction walls
  a_domainBdryType[2] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL; //+-z direction walls

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferBdryType( pde_buffer_bdry_type_e& pdeBufferBdryType ) {
  /* MODEL START */

  /* Set the Boundary Type at between Buffer space, and Simulation space
    We don't use buffers in this model, so set as PDE_BUFFER_BDRY_TYPE_HARD_WALL */
  pdeBufferBdryType = PDE_BUFFER_BDRY_TYPE_HARD_WALL;

  /* MODEL END */

  return;
}

void ModelRoutine::updateTimeStepInfo( TimeStepInfo& timeStepInfo ) {
  /* MODEL START */

  /* There are two different TimeStep sizes in Biocellion:
    Baseline time step and the finer, num state and grid time step */
  timeStepInfo.durationBaselineTimeStep = BASELINE_TIME_STEP_DURATION;
  timeStepInfo.numStateAndGridTimeStepsPerBaseline = NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;

  /* MODEL END */

  return;
}

void ModelRoutine::updateSyncMethod( sync_method_e& mechIntrctSyncMethod, sync_method_e& updateIfGridVarSyncMethod/* dummy if both callUpdateIfGridVarPreStateAndGridStep and callUpdateIfGridVarPostStateAndGridStep are set to false in ModelRoutine::updateOptModelRoutineCallInfo */ ) {
  /* MODEL START */

  /* Update the sync method in model_routine_mech_intrct.cpp and model_routine_grid.cpp
   Since mechanical interactions and PDEs are not used in this model, these are dummy */
  mechIntrctSyncMethod = SYNC_METHOD_PER_ATTR;
  updateIfGridVarSyncMethod = SYNC_METHOD_PER_ATTR;

  /* MODEL END */

  return;
}

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentInfo( Vector<SpAgentInfo>& v_spAgentInfo ) {/* set the mechanical interaction range & the numbers of model specific variables */
  /* MODEL START */

  const auto params = readXMLParameters();
  const S32 numGenes = std::stoi(params[NUM_GENES_PARAM]);
  const S32 numCellTypes = std::stoi(params[NUM_CELL_TYPES_PARAM]);
  const REAL ifGridSpacing = std::stod(params[IF_GRID_SPACING_PARAM]);

  /* Provide information about the discrete agent types in the user model */
  SpAgentInfo info;

  for (S32 cellType = 0; cellType < numCellTypes; cellType++) {
    info.dMax = ifGridSpacing;
    info.numBoolVars = numGenes;
    info.numStateModelReals = 0;
    info.numStateModelInts = 1; // For cell id
    info.v_mechIntrctModelRealInfo.clear();
    info.v_mechIntrctModelIntInfo.clear();
    info.v_odeNetInfo.clear();

    v_spAgentInfo.push_back(info);
  }

  /* MODEL END */

  return;
}
#endif

void ModelRoutine::updateJunctionEndInfo( Vector<JunctionEndInfo>& v_junctionEndInfo ) {/* set the numbers of model specific variables */
  /* MODEL START */

  /* No junctions in this model */
  v_junctionEndInfo.clear();

  /* MODEL END */

  return;
}

void ModelRoutine::updatePhiPDEInfo( Vector<PDEInfo>& v_phiPDEInfo ) {
  /* MODEL START */

  const auto params = readXMLParameters();
  const S32 numCytokines = std::stoi(params[NUM_CYTOKINES_PARAM]);

  PDEInfo pdeInfo;
  GridPhiInfo gridPhiInfo;

  // Input signal
  pdeInfo.pdeIdx = 0;
  pdeInfo.pdeType = PDE_TYPE_REACTION_DIFFUSION_STEADY_STATE_LINEAR;
  pdeInfo.numLevels = NUM_AMR_LEVELS;
  pdeInfo.numTimeSteps = 0;
  pdeInfo.v_tagExpansionSize.assign(NUM_AMR_LEVELS, 0); /* v_tagExpansionSize[0] should be always 0 */
  pdeInfo.ifLevel = NUM_AMR_LEVELS - 1;
  pdeInfo.mgSolveInfo.numPre = 3;
  pdeInfo.mgSolveInfo.numPost = 3;
  pdeInfo.mgSolveInfo.numBottom = 3;
  pdeInfo.mgSolveInfo.vCycle = true;
  pdeInfo.mgSolveInfo.maxIters = 30;
  pdeInfo.mgSolveInfo.epsilon = 1e-8;
  pdeInfo.mgSolveInfo.hang = 1e-6;
  pdeInfo.mgSolveInfo.normThreshold = 1e-10;

  pdeInfo.advectionInfo.courantNumber = 0.5; /* dummy */

  pdeInfo.splittingInfo.v_diffusionTimeSteps.assign(1, 1); /* dummy */
  pdeInfo.splittingInfo.odeStiff = ODE_STIFF_NORMAL; /* dummy */
  pdeInfo.splittingInfo.odeH = 0.5; /* dummy */
  pdeInfo.splittingInfo.odeHm = 0.1; /* dummy */
  pdeInfo.splittingInfo.odeEpsilon = 1e-6; /* dummy */
  pdeInfo.splittingInfo.odeThreshold = 1e-19; /* dummy */

  pdeInfo.callAdjustRHSTimeDependentLinear = false;

  gridPhiInfo.elemIdx = 0;
  gridPhiInfo.name = "input_substance";
  gridPhiInfo.syncMethod = VAR_SYNC_METHOD_DELTA;
  gridPhiInfo.aa_bcType[0][0] = BC_TYPE_DIRICHLET_MODEL;
  gridPhiInfo.aa_bcVal[0][0] = 0.0;
  gridPhiInfo.aa_bcType[0][1] = BC_TYPE_DIRICHLET_MODEL;
  gridPhiInfo.aa_bcVal[0][1] = 0.0;
  gridPhiInfo.aa_bcType[1][0] = BC_TYPE_DIRICHLET_MODEL;
  gridPhiInfo.aa_bcVal[1][0] = 0.0;
  gridPhiInfo.aa_bcType[1][1] = BC_TYPE_DIRICHLET_MODEL;
  gridPhiInfo.aa_bcVal[1][1] = 0.0;
  gridPhiInfo.aa_bcType[2][0] = BC_TYPE_DIRICHLET_MODEL;
  gridPhiInfo.aa_bcVal[2][0] = 0.0;
  gridPhiInfo.aa_bcType[2][1] = BC_TYPE_DIRICHLET_MODEL;
  gridPhiInfo.aa_bcVal[2][1] = 0.0;
  gridPhiInfo.errorThresholdVal = GRID_PHI_NORM_THRESHOLD * -1.0;
  gridPhiInfo.warningThresholdVal = GRID_PHI_NORM_THRESHOLD * -1.0;
  gridPhiInfo.setNegToZero = true;
  pdeInfo.v_gridPhiInfo.clear();
  pdeInfo.v_gridPhiInfo.push_back(gridPhiInfo);
  v_phiPDEInfo.push_back(pdeInfo);

  // Initialize cytokine PDEs
  for (S32 cytokine = 0; cytokine < numCytokines; cytokine++) {
    pdeInfo.pdeIdx = 1 + cytokine;
    pdeInfo.pdeType = PDE_TYPE_REACTION_DIFFUSION_STEADY_STATE_LINEAR;
    pdeInfo.numLevels = NUM_AMR_LEVELS;
    pdeInfo.numTimeSteps = 0;
    pdeInfo.v_tagExpansionSize.assign(NUM_AMR_LEVELS, 0); /* v_tagExpansionSize[0] should be always 0 */
    pdeInfo.ifLevel = NUM_AMR_LEVELS - 1;
    pdeInfo.mgSolveInfo.numPre = 3;
    pdeInfo.mgSolveInfo.numPost = 3;
    pdeInfo.mgSolveInfo.numBottom = 3;
    pdeInfo.mgSolveInfo.vCycle = true;
    pdeInfo.mgSolveInfo.maxIters = 30;
    pdeInfo.mgSolveInfo.epsilon = 1e-8;
    pdeInfo.mgSolveInfo.hang = 1e-6;
    pdeInfo.mgSolveInfo.normThreshold = 1e-10;

    pdeInfo.advectionInfo.courantNumber = 0.5; /* dummy */

    pdeInfo.splittingInfo.v_diffusionTimeSteps.assign(1, 1); /* dummy */
    pdeInfo.splittingInfo.odeStiff = ODE_STIFF_NORMAL; /* dummy */
    pdeInfo.splittingInfo.odeH = 0.5; /* dummy */
    pdeInfo.splittingInfo.odeHm = 0.1; /* dummy */
    pdeInfo.splittingInfo.odeEpsilon = 1e-6; /* dummy */
    pdeInfo.splittingInfo.odeThreshold = 1e-19; /* dummy */

    pdeInfo.callAdjustRHSTimeDependentLinear = false;

    gridPhiInfo.elemIdx = 1 + cytokine;
    gridPhiInfo.name = "cytokine" + std::to_string(cytokine);
    gridPhiInfo.syncMethod = VAR_SYNC_METHOD_DELTA;
    gridPhiInfo.aa_bcType[0][0] = BC_TYPE_NEUMANN_CONST; /* dummy */
    gridPhiInfo.aa_bcVal[0][0] = 0.0;                    /* dummy */
    gridPhiInfo.aa_bcType[0][1] = BC_TYPE_NEUMANN_CONST; /* dummy */
    gridPhiInfo.aa_bcVal[0][1] = 0.0;                    /* dummy */
    gridPhiInfo.aa_bcType[1][0] = BC_TYPE_NEUMANN_CONST; /* dummy */
    gridPhiInfo.aa_bcVal[1][0] = 0.0;                    /* dummy */
    gridPhiInfo.aa_bcType[1][1] = BC_TYPE_NEUMANN_CONST; /* dummy */
    gridPhiInfo.aa_bcVal[1][1] = 0.0;                    /* dummy */
    gridPhiInfo.aa_bcType[2][0] = BC_TYPE_NEUMANN_CONST; /* dummy */
    gridPhiInfo.aa_bcVal[2][0] = 0.0;                    /* dummy */
    gridPhiInfo.aa_bcType[2][1] = BC_TYPE_NEUMANN_CONST; /* dummy */
    gridPhiInfo.aa_bcVal[2][1] = 0.0;                    /* dummy */
    gridPhiInfo.errorThresholdVal = GRID_PHI_NORM_THRESHOLD * -1.0;
    gridPhiInfo.warningThresholdVal = GRID_PHI_NORM_THRESHOLD * -1.0;
    gridPhiInfo.setNegToZero = true;
    pdeInfo.v_gridPhiInfo.clear();
    pdeInfo.v_gridPhiInfo.push_back(gridPhiInfo);
    v_phiPDEInfo.push_back(pdeInfo);
  }

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfGridModelVarInfo( Vector<IfGridModelVarInfo>& v_ifGridModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridModelIntInfo ) {
  /* MODEL START */

  const auto params = readXMLParameters();
  const S32 numCytokines = std::stoi(params[NUM_CYTOKINES_PARAM]);

  // Initialize RHS for input and cytokines
  v_ifGridModelRealInfo.clear();
  IfGridModelVarInfo info;

  info.name = "input_substance_rhs";
  info.syncMethod = VAR_SYNC_METHOD_DELTA;
  v_ifGridModelRealInfo.push_back(info);

  for (S32 cytokine = 0; cytokine < numCytokines; cytokine++) {
    info.name = std::to_string(cytokine) + "_rhs";
    info.syncMethod = VAR_SYNC_METHOD_DELTA;
    v_ifGridModelRealInfo.push_back(info);
  }

  // Voxels have the occupied volume property in order to calculate kappa
  info.name = "occupied_volume";
  info.syncMethod = VAR_SYNC_METHOD_OVERWRITE;
  v_ifGridModelRealInfo.push_back(info);

  v_ifGridModelIntInfo.clear();

  /* MODEL END */

  return;
}

void ModelRoutine::updateRNGInfo( Vector<RNGInfo>& v_rngInfo ) {
  /* MODEL START */

  /* We use two different RGN's in this model */
  CHECK(NUM_MODEL_RNGS == 2);

  v_rngInfo.resize(NUM_MODEL_RNGS);

  RNGInfo rngInfo;
  /* Uniform distribution (min=0 and max=1) */
  rngInfo.type = RNG_TYPE_UNIFORM;
  rngInfo.param0 = 0.0;
  rngInfo.param1 = 1.0;
  rngInfo.param2 = 0.0;/* dummy */
  v_rngInfo[MODEL_RNG_UNIFORM] = rngInfo;

  /* Gaussian distribution (mean=0 and std=1) */
  rngInfo.type = RNG_TYPE_GAUSSIAN;
  rngInfo.param0 = 0.0;
  rngInfo.param1 = 1.0;
  rngInfo.param2 = 0.0;/* dummy */
  v_rngInfo[MODEL_RNG_GAUSSIAN] = rngInfo;

  /* MODEL END */

  return;
}

void ModelRoutine::updateFileOutputInfo( FileOutputInfo& fileOutputInfo ) {
  /* MODEL START */

  const auto params = readXMLParameters();
  const auto numGenes = std::stoi(params[NUM_GENES_PARAM]);
  const auto numCytokines = std::stoi(params[NUM_CYTOKINES_PARAM]);

  /* FileOutputInfo class holds the information related to file output of simulation results. */
  fileOutputInfo.particleOutput = true;
  fileOutputInfo.v_particleExtraOutputScalarVarName.clear();
  // Gene values are outputted every timestep and are packed into bits of a float64
  // (Biocellion doesn't seem to support custom int outputs)
  for (S32 i = 0; i < std::ceil(double(numGenes) / (sizeof(REAL) * 8)); i++) {
    fileOutputInfo.v_particleExtraOutputScalarVarName.push_back("genebits_" + std::to_string(i));
  }
  fileOutputInfo.v_particleExtraOutputVectorVarName.clear();
  // Don't output phi for every step
  fileOutputInfo.v_gridPhiOutput.assign(1 + numCytokines, false);
  fileOutputInfo.v_gridPhiOutputDivideByKappa.assign(1 + numCytokines, false);

  /* MODEL END */

  return;
}

void ModelRoutine::updateSummaryOutputInfo( Vector<SummaryOutputInfo>& v_summaryOutputRealInfo, Vector<SummaryOutputInfo>& v_summaryOutputIntInfo ) {
  /* MODEL START */

  /* Declare information you want to output when you run the simulation
    Biocellion supports a summary (reduction) mechanism for interface grid
    variables. Users add a fixed number of grid state variables to every unit box in the interface
    grid (or assign a fixed number of attributes to the interface grid), and for each attribute, users
    can ask Biocellion to visit every unit box to reduce the variables for the attribute to a single
    value. summary type e is used to set the reduction method. Choose one of these summary types:
    {SUMMARY_TYPE_SUM, SUMMARY_TYPE_AVG, SUMMARY_TYPE_MIN, SUMMARY_TYPE_MAX} */

  const auto params = readXMLParameters();
  const auto numCytokines = std::stoi(params[NUM_CYTOKINES_PARAM]);

  SummaryOutputInfo info;
  v_summaryOutputIntInfo.clear();
  v_summaryOutputRealInfo.clear();

  if (ENABLE_SUMMARY) {
    info.name = "input_signal_min";
    info.type = SUMMARY_TYPE_MIN;
    v_summaryOutputRealInfo.push_back(info);
    info.name = "input_signal_max";
    info.type = SUMMARY_TYPE_MAX;
    v_summaryOutputRealInfo.push_back(info);
    info.name = "input_signal_avg";
    info.type = SUMMARY_TYPE_AVG;
    v_summaryOutputRealInfo.push_back(info);

    for (S32 cytokine = 0; cytokine < numCytokines; cytokine++) {
      info.name = "cytokine" + std::to_string(cytokine) + "_min";
      info.type = SUMMARY_TYPE_MIN;
      v_summaryOutputRealInfo.push_back(info);
      info.name = "cytokine" + std::to_string(cytokine) + "_max";
      info.type = SUMMARY_TYPE_MAX;
      v_summaryOutputRealInfo.push_back(info);
      info.name = "cytokine" + std::to_string(cytokine) + "_avg";
      info.type = SUMMARY_TYPE_AVG;
      v_summaryOutputRealInfo.push_back(info);
    }
  }

  /* MODEL END */

  return;
}

void ModelRoutine::initGlobal( Vector<U8>& v_globalData ) {
  /* MODEL START */

  // Read all the arguments passed to Biocellion
  const auto params = readXMLParameters();
  const S32 numGenes = std::stoi(params[NUM_GENES_PARAM]);
  const S32 numCells = std::stoi(params[NUM_CELLS_PARAM]);
  const S32 numCellTypes = std::stoi(params[NUM_CELL_TYPES_PARAM]);
  auto nvFile = std::ifstream(params[NV_FILE_PARAM]);
  auto varfFile = std::ifstream(params[VARF_FILE_PARAM]);
  auto ttFile = std::ifstream(params[TT_FILE_PARAM]);
  auto geneInitialStatesFile = std::ifstream(params[GENE_INITIAL_STATES_FILE_PARAM]);
  auto inputSignalFile = std::ifstream(params[INPUT_SIGNAL_FILE_PARAM]);
  const REAL alphaCytokines = std::stod(params[ALPHA_CYTOKINES_PARAM]);
  const REAL betaCytokines = std::stod(params[BETA_CYTOKINES_PARAM]);
  const S32 numCytokines = std::stoi(params[NUM_CYTOKINES_PARAM]);
  const REAL secretionLow = std::stod(params[SECRETION_LOW_PARAM]);
  const REAL secretionHigh = std::stod(params[SECRETION_HIGH_PARAM]);
  const REAL cytokineThreshold = std::stod(params[CYTOKINE_THRESHOLD_PARAM]);
  const REAL cellRadius = std::stod(params[CELL_RADIUS_PARAM]);
  const REAL ifGridSpacing = std::stod(params[IF_GRID_SPACING_PARAM]);
  const REAL dirichletBoundary = std::stod(params[DIRICHLET_BOUNDARY_PARAM]);
  const REAL inputThreshold = std::stod(params[INPUT_THRESHOLD_PARAM]);
  const REAL alphaInput = std::stod(params[ALPHA_INPUT_PARAM]);
  const REAL betaInput = std::stod(params[BETA_INPUT_PARAM]);
  const S32 zLayers = std::stoi(params[Z_LAYERS_PARAM]);
  const S32 yLayers = std::stoi(params[Y_LAYERS_PARAM]);
  const S32 xLayers = std::stoi(params[X_LAYERS_PARAM]);
  const BOOL kappa = (params[KAPPA_PARAM] == "True" ? true : false);

  Vector<S32> nv;
  Vector<S32> varfOffsets;
  Vector<S32> ttOffsets;
  Vector<S32> varf;
  Vector<S32> tt;
  Vector<S32> geneInitialStates;
  Vector<S32> inputSignal;

  S32 i;
  std::string line;

  while (nvFile >> i) {
    nv.push_back(i);
  }

  S32 currentVarfOffset = 0;
  for (const auto v : nv) {
    varfOffsets.push_back(currentVarfOffset);
    currentVarfOffset += v;
  }

  S32 currentTtOffset = 0;
  for (const auto v : nv) {
    ttOffsets.push_back(currentTtOffset);
    if (v > 0) {
      currentTtOffset += std::pow(2, v);
    }
  }

  while (varfFile >> i) {
    varf.push_back(i);
  }

  while (ttFile >> i) {
    tt.push_back(i);
  }

  while (geneInitialStatesFile >> i) {
    geneInitialStates.push_back(i);
  }

  while (inputSignalFile >> i) {
    inputSignal.push_back(i);
  }

  // Setup offsets within global data
  GlobalDataFormat format;
  S64 size = sizeof(format);

  format.numGenes = size;
  size += sizeof(numGenes);

  format.numCells = size;
  size += sizeof(numCells);

  format.numCellTypes = size;
  size += sizeof(numCellTypes);

  format.nv = size;
  size += nv.size() * sizeof(nv[0]);

  format.varfOffsets = size;
  size += varfOffsets.size() * sizeof(varfOffsets[0]);

  format.ttOffsets = size;
  size += ttOffsets.size() * sizeof(ttOffsets[0]);

  format.varf = size;
  size += varf.size() * sizeof(varf[0]);

  format.tt = size;
  size += tt.size() * sizeof(tt[0]);

  format.geneInitialStates = size;
  size += geneInitialStates.size() * sizeof(geneInitialStates[0]);

  format.inputSignal = size;
  size += inputSignal.size() * sizeof(inputSignal[0]);

  format.alphaInput = size;
  size += sizeof(alphaInput);

  format.betaInput = size;
  size += sizeof(betaInput);

  format.alphaCytokines = size;
  size += sizeof(alphaCytokines);

  format.betaCytokines = size;
  size += sizeof(betaCytokines);

  format.numCytokines = size;
  size += sizeof(numCytokines);

  format.secretionLow = size;
  size += sizeof(secretionLow);

  format.secretionHigh = size;
  size += sizeof(secretionHigh);

  format.cytokineThreshold = size;
  size += sizeof(cytokineThreshold);

  format.cellRadius = size;
  size += sizeof(cellRadius);

  format.ifGridSpacing = size;
  size += sizeof(ifGridSpacing);

  format.dirichletBoundary = size;
  size += sizeof(dirichletBoundary);

  format.inputThreshold = size;
  size += sizeof(inputThreshold);

  format.zLayers = size;
  size += sizeof(zLayers);

  format.yLayers = size;
  size += sizeof(yLayers);

  format.xLayers = size;
  size += sizeof(xLayers);

  format.kappa = size;
  size += sizeof(kappa);

  // Copy passed arguments into global data so it can be accessed in other processes
  v_globalData.resize(size);
  memcpy(&(v_globalData[0]), &format, sizeof(format));
  memcpy(&(v_globalData[format.numGenes]), &numGenes, sizeof(numGenes));
  memcpy(&(v_globalData[format.numCells]), &numCells, sizeof(numCells));
  memcpy(&(v_globalData[format.numCellTypes]), &numCellTypes, sizeof(numCellTypes));
  memcpy(&(v_globalData[format.nv]), &(nv[0]), nv.size() * sizeof(nv[0]));
  memcpy(&(v_globalData[format.varfOffsets]), &(varfOffsets[0]), varfOffsets.size() * sizeof(varfOffsets[0]));
  memcpy(&(v_globalData[format.ttOffsets]), &(ttOffsets[0]), ttOffsets.size() * sizeof(ttOffsets[0]));
  memcpy(&(v_globalData[format.varf]), &(varf[0]), varf.size() * sizeof(varf[0]));
  memcpy(&(v_globalData[format.tt]), &(tt[0]), tt.size() * sizeof(tt[0]));
  memcpy(&(v_globalData[format.geneInitialStates]), &(geneInitialStates[0]), geneInitialStates.size() * sizeof(geneInitialStates[0]));
  memcpy(&(v_globalData[format.inputSignal]), &(inputSignal[0]), inputSignal.size() * sizeof(inputSignal[0]));
  memcpy(&(v_globalData[format.alphaInput]), &alphaInput, sizeof(alphaInput));
  memcpy(&(v_globalData[format.betaInput]), &betaInput, sizeof(betaInput));
  memcpy(&(v_globalData[format.alphaCytokines]), &alphaCytokines, sizeof(alphaCytokines));
  memcpy(&(v_globalData[format.betaCytokines]), &betaCytokines, sizeof(betaCytokines));
  memcpy(&(v_globalData[format.numCytokines]), &numCytokines, sizeof(numCytokines));
  memcpy(&(v_globalData[format.secretionLow]), &secretionLow, sizeof(secretionLow));
  memcpy(&(v_globalData[format.secretionHigh]), &secretionHigh, sizeof(secretionHigh));
  memcpy(&(v_globalData[format.cytokineThreshold]), &cytokineThreshold, sizeof(cytokineThreshold));
  memcpy(&(v_globalData[format.cellRadius]), &cellRadius, sizeof(cellRadius));
  memcpy(&(v_globalData[format.ifGridSpacing]), &ifGridSpacing, sizeof(ifGridSpacing));
  memcpy(&(v_globalData[format.dirichletBoundary]), &dirichletBoundary, sizeof(dirichletBoundary));
  memcpy(&(v_globalData[format.inputThreshold]), &inputThreshold, sizeof(inputThreshold));
  memcpy(&(v_globalData[format.zLayers]), &zLayers, sizeof(zLayers));
  memcpy(&(v_globalData[format.yLayers]), &yLayers, sizeof(yLayers));
  memcpy(&(v_globalData[format.xLayers]), &xLayers, sizeof(xLayers));
  memcpy(&(v_globalData[format.kappa]), &kappa, sizeof(kappa));

  /* MODEL END */

  return;
}

void ModelRoutine::init( void ) {
  /* MODEL START */

  // Load arguments from global data into global variables for efficient access
  const auto& g = Info::getGlobalDataRef();
  const GlobalDataFormat f = *((GlobalDataFormat*)(&(g[0])));

  gNumGenes = *((S32*)(&(g[f.numGenes])));
  gNumCells = *((S32*)(&(g[f.numCells])));
  gNumCellTypes = *((S32*)(&(g[f.numCellTypes])));
  gNv = (S32*)(&(g[f.nv]));
  gVarfOffsets = (S32*)(&(g[f.varfOffsets]));
  gTtOffsets = (S32*)(&(g[f.ttOffsets]));
  gVarf = (S32*)(&(g[f.varf]));
  gTt = (S32*)(&(g[f.tt]));
  gGeneInitialStates = (S32*)(&(g[f.geneInitialStates]));
  gInputSignal = (S32*)(&(g[f.inputSignal]));
  gAlphaInput = *((REAL*)(&(g[f.alphaInput])));
  gBetaInput = *((REAL*)(&(g[f.betaInput])));
  gAlphaCytokines = *((REAL*)(&(g[f.alphaCytokines])));
  gBetaCytokines = *((REAL*)(&(g[f.betaCytokines])));
  gNumCytokines = *((S32*)(&(g[f.numCytokines])));
  gSecretionLow = *((REAL*)(&(g[f.secretionLow])));
  gSecretionHigh = *((REAL*)(&(g[f.secretionHigh])));
  gCytokineThreshold = *((REAL*)(&(g[f.cytokineThreshold])));
  gCellRadius = *((REAL*)(&(g[f.cellRadius])));
  gIfGridSpacing = *((REAL*)(&(g[f.ifGridSpacing])));
  gDirichletBoundary = *((REAL*)(&(g[f.dirichletBoundary])));
  gInputThreshold = *((REAL*)(&(g[f.inputThreshold])));
  gZLayers = *((S32*)(&(g[f.zLayers])));
  gYLayers = *((S32*)(&(g[f.yLayers])));
  gXLayers = *((S32*)(&(g[f.xLayers])));
  gKappa = *((BOOL*)(&(g[f.kappa])));

  /* MODEL END */

  return;
}

void ModelRoutine::term( void ) {
  /* MODEL START */

  // Empty

  /* MODEL END */

  return;
}

void ModelRoutine::setPDEBuffer( const VIdx& startVIdx, const VIdx& regionSize, BOOL& isPDEBuffer ) {
  /* MODEL START */

  isPDEBuffer = false;

  /* MODEL END */

  return;
}

void ModelRoutine::setHabitable( const VIdx& vIdx, BOOL& isHabitable ) {
  /* MODEL START */

  /* Used to check if the initialized/simulated cells are placed in habitable space */
  isHabitable = true;

  /* MODEL END */

  return;
}