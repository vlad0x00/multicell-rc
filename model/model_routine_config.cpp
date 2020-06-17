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

static Vector<string> readXMLParameters() {
	Vector<string> v_modelParam; // This will be returned

	std::string paramsString = Info::getModelParam();
	std::string delimiter = " ";

	size_t pos = 0;
	std::string token;
	while ((pos = paramsString.find(delimiter)) != std::string::npos) {
		token = paramsString.substr(0, pos);
		if (token == "") { continue; }
		v_modelParam.push_back(token);
		paramsString.erase(0, pos + delimiter.length());
	}
	if (!paramsString.empty()) {
		v_modelParam.push_back(paramsString);
	}

	return v_modelParam;
}

void ModelRoutine::updateIfGridSpacing( REAL& ifGridSpacing ) {
	/* MODEL START */

  /* Initialize Interface Grid Spacing (which was declared in model_define.h) */
  ifGridSpacing = IF_GRID_SPACING;

	/* MODEL END */

	return;
}

void ModelRoutine::updateOptModelRoutineCallInfo( OptModelRoutineCallInfo& callInfo ) {
	/* MODEL START */

	// Empty

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

	auto params = readXMLParameters();

	/* Provide information about the discrete agent types in the user model */
	SpAgentInfo info;

	info.dMax = IF_GRID_SPACING;
	info.numBoolVars = std::stoi(params[0]);
	info.numStateModelReals = 0;
	info.numStateModelInts = 0;
	info.v_mechIntrctModelRealInfo.clear();
	info.v_mechIntrctModelIntInfo.clear();
	info.v_odeNetInfo.clear();

	v_spAgentInfo.push_back(info);

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

	/* No PDE's in this model */
	v_phiPDEInfo.clear();

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridModelVarInfo( Vector<IfGridModelVarInfo>& v_ifGridModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridModelIntInfo ) {
	/* MODEL START */

	/* No model_routine_grid.cpp in this model */
  v_ifGridModelRealInfo.clear();
  v_ifGridModelIntInfo.clear();

	/* MODEL END */

	return;
}

void ModelRoutine::updateRNGInfo( Vector<RNGInfo>& v_rngInfo ) {
	/* MODEL START */

	/* We use two different RGN's in this model */
	CHECK( NUM_MODEL_RNGS == 2 );

	v_rngInfo.resize( NUM_MODEL_RNGS );

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

	/* FileOutputInfo class holds the information related to file output of simulation results. */
	fileOutputInfo.particleOutput = true;                          
	fileOutputInfo.v_particleExtraOutputScalarVarName.clear();
	fileOutputInfo.v_particleExtraOutputVectorVarName.clear();
	fileOutputInfo.v_gridPhiOutput.clear();
	fileOutputInfo.v_gridPhiOutputDivideByKappa.clear();

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

	Vector<string> v_modelParam = readXMLParameters();
	S32 numGenes = std::stoi(v_modelParam[0]);

	SummaryOutputInfo info;
	v_summaryOutputIntInfo.resize(2 + numGenes);
	v_summaryOutputRealInfo.clear();

	info.name = "Timestep";
	info.type = SUMMARY_TYPE_SUM;
	v_summaryOutputIntInfo[0] = info;

	info.name = "Input";
	info.type = SUMMARY_TYPE_SUM;
	v_summaryOutputIntInfo[1] = info;

	for (S32 i = 0; i < numGenes; i++) {
		info.name = "Gene " + std::to_string(i);
		info.type = SUMMARY_TYPE_SUM;
		v_summaryOutputIntInfo[i + 2] = info;
	}
	
	/* MODEL END */

	return;
}

void ModelRoutine::initGlobal( Vector<U8>& v_globalData ) {
	/* MODEL START */

	Vector<string> v_modelParam = readXMLParameters();

	S32 param = 0;
	const S32 numGenes = std::stoi(v_modelParam[param++]);
	const S32 numCells = std::stoi(v_modelParam[param++]);
  auto nvFile = std::ifstream(v_modelParam[param++]);
	auto varfFile = std::ifstream(v_modelParam[param++]);
	auto ttFile = std::ifstream(v_modelParam[param++]);
	auto bnInitialStateFile = std::ifstream(v_modelParam[param++]);
	auto inputArrayFile = std::ifstream(v_modelParam[param++]);

	Vector<S32> nv;
	Vector<S32> varfOffsets;
	Vector<S32> ttOffsets;
	Vector<S32> varf;
	Vector<S32> tt;
	Vector<S32> bnInitialState;
	Vector<S32> inputArray;

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

	while (bnInitialStateFile >> i) {
		bnInitialState.push_back(i);
	}

	while (inputArrayFile >> i) {
		inputArray.push_back(i);
	}

	GlobalDataFormat format;
	S64 size = sizeof(format);

	format.numGenes = size;
	size += sizeof(numGenes);

	format.numCells = size;
	size += sizeof(numCells);

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

	format.bnInitialState = size;
	size += bnInitialState.size() * sizeof(bnInitialState[0]);

	format.inputArray = size;
	size += inputArray.size() * sizeof(inputArray[0]);

	v_globalData.resize(size);
	memcpy(&(v_globalData[0]), &format, sizeof(format));
	memcpy(&(v_globalData[format.numGenes]), &numGenes, sizeof(numGenes));
	memcpy(&(v_globalData[format.numCells]), &numCells, sizeof(numCells));
	memcpy(&(v_globalData[format.nv]), &(nv[0]), nv.size() * sizeof(nv[0]));
	memcpy(&(v_globalData[format.varfOffsets]), &(varfOffsets[0]), varfOffsets.size() * sizeof(varfOffsets[0]));
	memcpy(&(v_globalData[format.ttOffsets]), &(ttOffsets[0]), ttOffsets.size() * sizeof(ttOffsets[0]));
	memcpy(&(v_globalData[format.varf]), &(varf[0]), varf.size() * sizeof(varf[0]));
	memcpy(&(v_globalData[format.tt]), &(tt[0]), tt.size() * sizeof(tt[0]));
	memcpy(&(v_globalData[format.bnInitialState]), &(bnInitialState[0]), bnInitialState.size() * sizeof(bnInitialState[0]));
	memcpy(&(v_globalData[format.inputArray]), &(inputArray[0]), inputArray.size() * sizeof(inputArray[0]));

	/* MODEL END */

	return;
}

void ModelRoutine::init( void ) {
	/* MODEL START */

	// Empty

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