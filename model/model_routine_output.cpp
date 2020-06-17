/*

Copyright © 2013 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* UESR START */

#include "model_define.h"

/* UESR END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentOutput( const VIdx& vIdx, const SpAgent& spAgent, REAL& color, Vector<REAL>& v_extraScalar, Vector<VReal>& v_extraVector ) {
	/* MODEL START */

	color = spAgent.state.getType();
	CHECK(v_extraScalar.size() == 0);
	CHECK(v_extraVector.size() == 0);

	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updateSummaryVar( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<REAL>& v_realVal/* [elemIdx] */, Vector<S32>& v_intVal/* [elemIdx] */ ) {
	/* MODEL START */

	S32 numGenes = getNumGenes();

	CHECK(v_realVal.size() == 0);
	CHECK(v_intVal.size() == 2 + (U32)numGenes);

	const UBAgentData& ubAgentData = *(nbrUBAgentData.getConstPtr(0, 0, 0));

	if (ubAgentData.v_spAgent.size() > 0) {
		const S32 baselineStep = Info::getCurBaselineTimeStep();
		v_intVal[0] = baselineStep;
		v_intVal[1] = getInputArray()[baselineStep];
		for (S32 i = 0; i < numGenes; i++) {
			for (auto agent : ubAgentData.v_spAgent) {
				v_intVal[i + 2] = agent.state.getBoolVal(i);
			}
		}
	} else {
		v_intVal[0] = 0;
		v_intVal[1] = 0;
		for (S32 i = 0; i < numGenes; i++) {
			v_intVal[i + 2] = 0;
		}
	}
  
	/* MODEL END */

	return;
}

