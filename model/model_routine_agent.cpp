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

#if HAS_SPAGENT
void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentOffset ) {/* initialization */
	/* MODEL START */

	if (init == true) {
		VIdx vIdx;      // The location of each cell defined by the unit box index, and position offset.
		VReal vOffset;  // Poisition offset is the vector distance from the center of the unit box.
		SpAgentState state;

		for (S32 cell = 0; cell < gNumCells; cell++) {
			S32 cellType = cell % gNumCellTypes;

			vIdx[0] = regionSize[0] * (0.5 * Util::getModelRand(MODEL_RNG_UNIFORM) + 0.250);
			vIdx[1] = regionSize[1] * (0.5 * Util::getModelRand(MODEL_RNG_UNIFORM) + 0.250);
			vIdx[2] = regionSize[2] * (0.5 * Util::getModelRand(MODEL_RNG_UNIFORM) + 0.250);

			vOffset[0] = IF_GRID_SPACING * (Util::getModelRand(MODEL_RNG_UNIFORM) - 0.5);
			vOffset[1] = IF_GRID_SPACING * (Util::getModelRand(MODEL_RNG_UNIFORM) - 0.5);
			vOffset[2] = IF_GRID_SPACING * (Util::getModelRand(MODEL_RNG_UNIFORM) - 0.5);

			/* Initialize state */
			state.setType(cellType);
			state.setRadius(CELL_RADIUS);

			for (S32 gene = 0; gene < gNumGenes; gene++) {
				state.setBoolVal(gene, getGeneInitialStates(cell)[gene]);
			}
			state.setModelInt(0, cell);

			/* Initialize return vectors {v_spAgentState, v_spAgentVIdx, v_spAgentOffset} by using .push_back() */
			CHECK(ifGridHabitableBoxData.get(vIdx) == true);
			v_spAgentVIdx.push_back(vIdx);
			v_spAgentOffset.push_back(vOffset);
			v_spAgentState.push_back(state);

		}
	}

	/* MODEL END */

	return;
}

void ModelRoutine::spAgentCRNODERHS( const S32 odeNetIdx, const VIdx& vIdx, const SpAgent& spAgent, const NbrUBEnv& nbrUBEnv, const Vector<double>& v_y, Vector<double>& v_f ) {
	/* MODEL START */

	// Empty

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */ ) {
	/* MODEL START */

	const auto cellType = state.getType();

	REAL aaaRatio[3][3][3];
	REAL* avgPhi = nullptr;
	if (gNumCytokines > 0) {
		avgPhi = new REAL[gNumCytokines];
		for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
			avgPhi[cytokine] = 0.0;
		}
		Util::computeSphereUBVolOvlpRatio(SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, vOffset, state.getRadius(), aaaRatio);
		for(S32 i = -1 ; i <= 1; i++) {
			for(S32 j = -1 ; j <= 1; j++) {
				for(S32 k = -1 ; k <= 1; k++) {
					for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
						avgPhi[cytokine] += nbrUBEnv.getPhi(i, j, k, cytokine) * aaaRatio[i + 1][j + 1][k + 1];
					}
				}
			}
		}
	}

	Vector<BOOL> newBools;

	const S32 baselineStep = Info::getCurBaselineTimeStep();
	newBools.push_back(gInputSignal[baselineStep + 1]);

	for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
		if (avgPhi[cytokine] > gCytokineThreshold) {
			newBools.push_back(1);	
		} else {
			newBools.push_back(0);
		}
	}

	for (S32 gene = 1 + gNumCytokines; gene < gNumGenes; gene++) {
		const auto nv = getNv(cellType)[gene];
		CHECK(nv >= 0);

		if (nv > 0 ) {
			const auto varf = getVarf(cellType, gene);
			const auto tt = getTt(cellType, gene);

			U32 ttEntry = 0;
			for (U32 j = 0; j < (U32)nv; j++) {
				CHECK(varf[j] >= 0 && varf[j] < gNumGenes);
				const auto var_val = state.getBoolVal(varf[j]);
				CHECK(var_val != -1);
				if (var_val) {
					ttEntry |= (1 << j);
				}
			}

			const auto tt_val = tt[ttEntry];
			CHECK(tt_val != -1);
			newBools.push_back(tt_val);
		} else {
			newBools.push_back(state.getBoolVal(gene));
		}
	}

	state.setBoolValArray(newBools);

	if (gNumCytokines > 0) {
		delete[] avgPhi;
	}

	/* MODEL END */

	return;
}

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentDisp ) {
	/* MODEL START */

	// Empty

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, BOOL& divide, BOOL& disappear ) {
	/* MODEL START */

 	divide = false;
  disappear = false;

	/* MODEL END */

	return;
}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, VReal& disp ) {/* if not dividing or disappearing */
	/* MODEL START */

	// Empty

	/* MODEL END */

	return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& motherState/* INOUT */, VReal& motherDisp, SpAgentState& daughterState, VReal& daughterDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

	// Empty

	/* MODEL END */

	return;
}
#endif

