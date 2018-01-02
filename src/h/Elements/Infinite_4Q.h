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

#include "Element.h"

using namespace std;

//! infinite element class
class CInfEle : public CElement
{
public:

//!	Constructor
	CInfEle();

//!	Desconstructor
	~CInfEle();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream 
	virtual void Write(COutputter& output, unsigned int Ele);

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* Q4stress, double* Displacement);

//! Calculate coordinate of gauss point
	virtual void ElementGauss(double* Coordinate);

//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();

#ifdef _VIB_
	//!	Calculate element mass matrix (Upper triangular matrix, stored as an array column by colum)
	virtual void ElementMass(double* mass);
#endif

};
