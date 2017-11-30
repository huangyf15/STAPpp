/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/TimoshenkoSRINT.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CTimoshenkoSRINT::CTimoshenkoSRINT()
{
	NEN = 2;	// Each element has 2 nodes
	nodes = new CNode*[NEN];

	ND = 12;
	LocationMatrix = new unsigned int[ND];

	ElementMaterial = nullptr;
}

//	Desconstructor
CTimoshenkoSRINT::~CTimoshenkoSRINT()
{
}

//	Read element data from stream Input
bool CTimoshenkoSRINT::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "    Expected element : " << Ele + 1 << endl
			<< "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
	ElementMaterial = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CTimoshenkoSRINT::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes[0]->NodeNumber
		<< setw(9) << nodes[1]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CTimoshenkoSRINT::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN; N++)
		for (unsigned int D = 0; D < 6; D++)
			LocationMatrix[i++] = nodes[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node Timoshenko beam element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int CTimoshenkoSRINT::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CTimoshenkoSRINT::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

}

//	Calculate element stress 
void CTimoshenkoSRINT::ElementStress(double stress[3], double force[12], double* Displacement)
{
	CBarMaterial* material = dynamic_cast<CBarMaterial*>(ElementMaterial);	// Pointer to material of the element
}
