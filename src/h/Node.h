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

#include "Outputter.h"

#include <iostream>
#include <fstream>

using namespace std;

//!	Node class
class CNode
{
public:

//!	Maximum number of degrees of freedom per node
/*! the nodes of structure elements need 6 degrees of freedom, and for nodes of solid elements, value of the last 3 bcodes are 1 */
	const static unsigned int NDF= 6;

//!	Node number
	unsigned int NodeNumber;

//!	x, y and z coordinates of the node
	double XYZ[3];

//!	Boundary code of each degree of freedom of the node
/*!		0: The corresponding degree of freedom is active (defined in the global system) */
/*!		1: The corresponding degree of freedom in nonactive (not defined) */
/*!	After call Domain::CalculateEquationNumber(), bcode stores the global equation number */
/*!	corresponding to each degree of freedom of the node */
	unsigned int bcode[NDF];

//! Boundary code flag marking whether the last 3 bcodes are given while inputting 
/*!     1: Given       !*/
/*!     0: Not given   !*/
	bool RotationDOFManuallyInputFlag;

//!	Constructor
	CNode(double X = 0, double Y = 0, double Z = 0);

//!	Read nodal point data from stream Input
	bool Read(ifstream& Input, unsigned int np);

//!	Output nodal point data to stream
	void Write(COutputter& output, unsigned int np);

//!	Output equation numbers of nodal point to stream OutputFile
	void WriteEquationNo(COutputter& OutputFile, unsigned int np);

//!	Write nodal displacement
	void WriteNodalDisplacement(COutputter& OutputFile, unsigned int np, double* Displacement);

#ifdef _VIB_
	void WriteVibrationDisplacement(COutputter& OutputFile, unsigned int np, double* Displacement);

#endif
};
