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

#include <cmath>

#if HAS_SPAGENT
void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentOffset ) {/* initialization */
  /* MODEL START */

  if (init == true) {
    VIdx vIdx;      // The location of each cell defined by the unit box index, and position offset.
    VReal vOffset;  // Position offset is the vector distance from the center of the unit box.
    SpAgentState state;

    const S32 tissueVolume = gXLayers * gYLayers * gZLayers;
    CHECK(tissueVolume >= gNumCells);

    // Setup a grid of cell locations
    Vector<VIdx> coordsVIdx;
    Vector<VReal> coordsOffset;
    for (S32 z = 0; z < gZLayers; z += 1) {
      for (S32 y = -gYLayers / 2; y < gYLayers - gYLayers / 2; y += 1) {
        for (S32 x = -gXLayers / 2; x < gXLayers - gXLayers / 2; x += 1) {
          vIdx[0] = x + regionSize[0] / 2;
          vIdx[1] = y + regionSize[1] / 2;
          vIdx[2] = z + gSourceDist;
          vOffset[0] = 0;
          vOffset[1] = 0;
          vOffset[2] = 0;
          coordsVIdx.push_back(vIdx);
          coordsOffset.push_back(vOffset);
        }
      }
    }

    CHECK(coordsVIdx.size() == (U64)(tissueVolume));
    CHECK(coordsOffset.size() == (U64)(tissueVolume));

    // Place all the cells and assign them type
    const S32 skipCellsPerLayer = gXLayers * gYLayers - gCellsPerLayer;
    CHECK(gNumCellTypes > 0);
    for (S32 cell = 0, cellLocation = 0; cell < gNumCells; cell++, cellLocation++) {
      S32 cellType = gNumCellTypes * Util::getModelRand(MODEL_RNG_UNIFORM);
      if (cellType == gNumCellTypes) { cellType--; }
      CHECK(cellType >= 0 && cellType < gNumCellTypes);
      std::cout << ("Cell: " + std::to_string(cell) + ", Type: " + std::to_string(cellType) + "\n") << std::flush;

      if (cellLocation % (gXLayers * gYLayers) >= gCellsPerLayer) {
        cellLocation += skipCellsPerLayer;
      }
      CHECK((U64)(cellLocation) < coordsVIdx.size() && (U64)(cellLocation) < coordsOffset.size());
      vIdx = coordsVIdx[cellLocation];
      vOffset = coordsOffset[cellLocation];

      /* Initialize state */
      state.setType(cellType);
      state.setRadius(gCellRadius);

      // Setup the initial state
      for (S32 gene = 0; gene < gNumGenes; gene++) {
        state.setBoolVal(gene, getGeneInitialStates(cell)[gene]);
      }
      // Set the id the of the cell
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

  // Empty

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

  const auto cellType = state.getType();

  // Determine input signal and cytokine levels in the cell's surroundings
  REAL aaaRatio[3][3][3];
  REAL* avgPhi = nullptr;
  avgPhi = new REAL[1 + gNumCytokines];
  avgPhi[0] = 0.0;
  for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
    avgPhi[1 + cytokine] = 0.0;
  }
  // Calculate the cell overlap with neighbouring voxels
  Util::computeSphereUBVolOvlpRatio(SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, vOffset, state.getRadius(), aaaRatio);
  // Iterate the surrounding voxels and sum up phi (input/cytokine levels)
  for(S32 i = -1 ; i <= 1; i++) {
    for(S32 j = -1 ; j <= 1; j++) {
      for(S32 k = -1 ; k <= 1; k++) {
        avgPhi[0] += nbrUBEnv.getPhi(i, j, k, 0) * aaaRatio[i + 1][j + 1][k + 1];
        for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
          avgPhi[1 + cytokine] += nbrUBEnv.getPhi(i, j, k, 1 + cytokine) * aaaRatio[i + 1][j + 1][k + 1];
        }
      }
    }
  }

  // Updated gene values
  Vector<BOOL> newBools;

  // Push new input signal gene value and update it immediately so that internal genes can see the change
  newBools.push_back((avgPhi[0] > gInputThreshold) ? 1 : 0);
  state.setBoolVal(0, newBools.back());

  // Push new cytokine input gene values and update them immediately for internal genes
  for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
    if (avgPhi[1 + cytokine] / gCytokineNormalization > gCytokineThreshold) {
      newBools.push_back(1);
    } else {
      newBools.push_back(0);
    }
    CHECK(S32(newBools.size() - 1) == 1 + cytokine);
    state.setBoolVal(newBools.size() - 1, newBools.back());
  }

  // Update cytokine output and internal genes
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

  delete[] avgPhi;

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

