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

	const auto numGenebits = std::ceil(double(gNumGenes) / (sizeof(REAL) * 8));

	color = spAgent.state.getType();
	CHECK(v_extraScalar.size() == size_t(numGenebits));
	CHECK(v_extraVector.size() == 0);

	CHECK(sizeof(REAL) == sizeof(U64));
	for (S32 genebits = 0; genebits < numGenebits; genebits++) {
		*((U64*)(&(v_extraScalar[genebits]))) = 0;
	}

	for (S32 gene = 0; gene < gNumGenes; gene++) {
		const S32 genebits = gene / (sizeof(REAL) * 8);
		const S32 genebitPos = gene % (sizeof(REAL) * 8);
		if (spAgent.state.getBoolVal(gene)) {
			*((U64*)(&(v_extraScalar[genebits]))) |= (1 << genebitPos);
		}
	}

	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updateSummaryVar( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<REAL>& v_realVal/* [elemIdx] */, Vector<S32>& v_intVal/* [elemIdx] */ ) {
	/* MODEL START */

	CHECK(v_realVal.size() == (U32)(gNumCytokines * 3));
	CHECK(v_intVal.size() == 0);

	for (S32 cytokine = 0; cytokine < gNumCytokines; cytokine++) {
		for (S32 i = 0; i < 3; i++) {
			v_realVal[cytokine * 3 + i] = nbrUBEnv.getPhi(0, 0, 0, cytokine);
		}
	}

	/* MODEL END */

	return;
}

