/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "CSRMatrix.h"
#include "SkylineMatrix.h"
#include "SparseMatrix.h"

//!	Base class for a solver
/*	New solver should be derived from this base class, and match the storage scheme
    of the global stiffness matrix employed in Domain class. */
class CSolver
{
public:
    CSolver(SparseMatrix<double>& _K){};
};

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver : public CSolver
{
protected:
    CSkylineMatrix<double>& K;

public:
    //!	Constructor
    CLDLTSolver(CSkylineMatrix<double>& _K) : CSolver(_K), K(_K){};

    //!	Perform L*D*L(T) factorization of the stiffness matrix
    void LDLT();

    //!	Reduce right-hand-side load vector and back substitute
    void BackSubstitution(double* Force);
#ifdef _VIB_
	void Multiple(double* acc,double* Force,unsigned int numeq,unsigned int vib_m);
#endif
};

class CSRSolver : public CSolver
{
protected:
    CSRMatrix<double>& K;

public:
    CSRSolver(CSRMatrix<double>& _K) : CSolver(_K), K(_K){};

    void solve(double* Force, unsigned NLCase);
};
