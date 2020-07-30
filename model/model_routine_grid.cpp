/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

#include "model_define.h"

using namespace std;

void ModelRoutine::initIfGridVar( const VIdx& vIdx, const UBAgentData& ubAgentData, UBEnv& ubEnv ) {
  /* MODEL START */

  ubEnv.setPhi(0, 0.0);
  ubEnv.setModelReal(0, 0.0);
  for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
    ubEnv.setPhi(1 + cytokine, 0.0);
    ubEnv.setModelReal(1 + cytokine, 0.0);
  }

  /* MODEL END */

  return;
}

void ModelRoutine::initIfSubgridKappa( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridKappa ) {/* relevant only if v_gridPhiOutputDivideByKappa[pdeIdx] is set to true in updateFileOutputInfo() */
  /* MODEL START */

  gridKappa = 1.0;

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfGridVar( const BOOL pre, const S32 iter, const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, NbrUBEnv& nbrUBEnv/* [INOUT] */ ) {
  /* MODEL START */

  CHECK(pre == true);
  CHECK(iter == 0);

  static const REAL dt = BASELINE_TIME_STEP_DURATION / NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE;

  REAL *rhs = nullptr;
  if (gNumCytokines > 0) {
    rhs = new REAL[gNumCytokines];

    for(S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
      rhs[cytokine] = 0.0;
    }

    /* iterate over 3 * 3 * 3 boxes */
    for(S32 i = -1; i <= 1; i++) {
      const auto x = i + vIdx[0];
      for(S32 j = -1; j <= 1; j++) {
        const auto y = j + vIdx[1];
        for(S32 k = -1; k <= 1; k++) {
          const auto z = k + vIdx[2];

          if (x < 0 || x >= Info::getDomainSize(0)) { continue; }
          if (y < 0 || y >= Info::getDomainSize(1)) { continue; }
          if (z < 0 || z >= Info::getDomainSize(2)) { continue; }

          const UBAgentData& ubAgentData = *(nbrUBAgentData.getConstPtr(i, j, k));
          VIdx ubVIdxOffset;
          ubVIdxOffset[0] = i* -1;
          ubVIdxOffset[1] = j* -1;
          ubVIdxOffset[2] = k* -1;

          for (const auto& agent: ubAgentData.v_spAgent) {
            const REAL ratio = Util::computeSphereUBVolOvlpRatio(SPHERE_UB_VOL_OVLP_RATIO_MAX_LEVEL, agent.vOffset, agent.state.getRadius(), ubVIdxOffset);

            if(ratio > 0.0) {
              // Don't really need those
              //REAL radius = agent.state.getRadius();
              //REAL vol = (4.0 * MY_PI / 3.0) * radius * radius * radius;

              for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
                rhs[cytokine] += ratio * dt * (agent.state.getBoolVal(1 + gNumCytokines + cytokine) ? gSecretionHigh : gSecretionLow);
              }
            }
          }
        }
      }
    }

    for(S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
      rhs[cytokine] /= (IF_GRID_SPACING * IF_GRID_SPACING * IF_GRID_SPACING);
      nbrUBEnv.setModelReal(0, 0, 0, 1 + cytokine, rhs[cytokine]);
    }

    delete[] rhs;
  }

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfSubgridKappa( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridKappa ) {
  /* MODEL START */

  gridKappa = 1.0;

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfSubgridAlpha( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridAlpha/* decay (-) */ ) {
  /* MODEL START */

  gridAlpha = -gAlpha;

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfSubgridBetaInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& subgridVOffset0, const UBAgentData& ubAgentData0, const UBEnv& ubEnv0, const VIdx& vIdx1, const VIdx& subgridVOffset1, const UBAgentData& ubAgentData1, const UBEnv& ubEnv1, REAL& gridBeta ) {
  /* MODEL START */

  gridBeta = gBeta;

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfSubgridBetaPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridBeta ) {
  /* MODEL START */

  gridBeta = gBeta;

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfSubgridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridBeta ) {
  /* MODEL START */

  gridBeta = gBeta;

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfSubgridRHSLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
  /* MODEL START */

  gridRHS = ubEnv.getModelReal(elemIdx);

  /* MODEL END */

  return;
}

void ModelRoutine::adjustIfSubgridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnvModelVar& ubEnvModelVar, const REAL gridPhi, REAL& gridRHS/* INOUT */ ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfSubgridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& vIdx, const VIdx& subgridVOffset, const UBAgentData& ubAgentData, const UBEnvModelVar& ubEnvModelVar, const Vector<double>& v_gridPhi/* [idx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {/* for Wnt & SFRP */
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfGridAMRTags( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<S32>& v_finestLevel/* [pdeIdx] */ ) {
  /* MODEL START */

  for(S32 pdeIdx = 0; pdeIdx < gNumCytokines; pdeIdx++) {
    v_finestLevel[pdeIdx] = NUM_AMR_LEVELS - 1;
  }

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfGridDirichletBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const UBEnvModelVar a_ubEnvModelVar[3], const Vector<REAL> av_gridPhi[3]/* av_gridPhi[].size() == ratio * raito * ratio (ratio = Info::envAuxDataInfo.v_phiRatioFromIfGridToIfSubgrid[elemIdx]), use VIdx::getIdx3DTo1D() to index */, REAL& bcVal ) {
  /* MODEL START */

  CHECK(elemIdx == 0);
  bcVal = gInputSignal[Info::getCurBaselineTimeStep() + 1] * INPUT_SIGNAL_DIRICHLET_VAL;

  /* MODEL END */

  return;
}

void ModelRoutine::updateIfGridNeumannBCVal( const S32 elemIdx, const VReal& pos, const S32 dim, const BOOL lowSide, const UBEnvModelVar a_ubEnvModelVar[3], const Vector<REAL> av_gridPhi[3]/* av_gridPhi[].size() == ratio * raito * ratio (ratio = Info::envAuxDataInfo.v_phiRatioFromIfGridToIfSubgrid[elemIdx]), use VIdx::getIdx3DTo1D() to index */, REAL& bcVal ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::initPDEBufferPhi( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, Vector<REAL>& v_gridPhi/* [idx] */ ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::initPDEBufferKappa( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridKappa ) {/* relevant only if v_gridPhiOutputDivideByKappa[pdeIdx] is set to true in updateFileOutputInfo() */
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferKappa( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridKappa ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferAlpha( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridAlpha/* decay (-) */ ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferBetaInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, REAL& gridBeta ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferRHSLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* uptake(-) and secretion (+) */ ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::adjustPDEBufferRHSTimeDependentLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const REAL gridPhi, REAL& gridRHS/* INOUT */ ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxSize, const Vector<double>& v_gridPhi/* [idx] */, Vector<double>& v_gridRHS/* [idx], uptake(-) and secretion (+) */ ) {
  /* MODEL START */

  ERROR( "unimplemented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferDirichletBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferFaceSize, const S32 dim, const BOOL lowSide, REAL& bcVal ) {
  /* MODEL START */

  ERROR( "unimplmented." );

  /* MODEL END */

  return;
}

void ModelRoutine::updatePDEBufferNeumannBCVal( const S32 elemIdx, const VReal& startPos, const VReal& pdeBufferFaceSize, const S32 dim, const BOOL lowSide, REAL& bcVal ) {
  /* MODEL START */

  ERROR( "unimplmented." );

  /* MODEL END */

  return;
}

