/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/Frustum.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CFrustum::CFrustum()
{
	NEN = 2;	// Each element has 2 nodes
	nodes = new CNode*[NEN];

	ND = 6;
	LocationMatrix = new unsigned int[ND];

	ElementMaterial = NULL;
}

//	Desconstructor
CFrustum::~CFrustum()
{
}

//	Read element data from stream Input
bool CFrustum::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	ElementMaterial = static_cast<CShellMaterial*>(MaterialSets) + MSet - 1;
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CFrustum::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes[0]->NodeNumber
		<< setw(9) << nodes[1]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CFrustum::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN; N++)
		for (unsigned int D = 0; D < ND; D++)
			LocationMatrix[i++] = nodes[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node frustum shell element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int CFrustum::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CFrustum::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//	Only X and Z matter
	double DX = nodes[1]->XYZ[0] - nodes[0]->XYZ[0];
	double DZ = nodes[1]->XYZ[2] - nodes[0]->XYZ[2];
	double l = sqrt(DX * DX + DZ * DZ);
	double S = DZ / l;
	double C = DX / l;
	double r = 0.5*(nodes[0]->XYZ[2] + nodes[1]->XYZ[2]);

	//	Calculate element stiffness matrix
	CShellMaterial* material = static_cast<CShellMaterial*>(ElementMaterial);	// Pointer to material of the element

	const double nu = material->nu;
	const double k = material->E / (1 - nu*nu);
	const double h = material->h;

	Matrix[0] = S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*3.5E1 + C*S*(l*l*l)*3.5E1)*(2.0 / 3.5E1)) / r) + C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) - S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	
	Matrix[3] = -C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*1.0 / (l*l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*3.5E1 + C*S*(l*l*l)*3.5E1)*(2.0 / 3.5E1)) / r) - S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) - S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[5] = -C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*3.5E1 + C*S*(l*l*l)*3.5E1)*(2.0 / 3.5E1)) / r) + S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) - S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	
	Matrix[10] = (h*k*3.141592653589793*((C*C)*(l*l*l*l) + (S*S)*(l*l)*1.4E1 + (h*h)*(r*r)*3.5E1)*(2.0 / 1.05E2)) / (l*r);
	Matrix[12] = (C*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.1E1 + (S*S)*(l*l)*2.1E1 + (h*h)*(r*r)*1.05E2)*(1.0 / 1.05E2)) / r - (S*h*k*l*3.141592653589793*((S*S)*5.0 - C*S*l*3.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[14] = (C*h*k*l*3.141592653589793*((S*S)*5.0 - C*S*l*3.0 + C*nu*r*5.0)*(-1.0 / 3.0E1)) / r - (S*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.1E1 + (S*S)*(l*l)*2.1E1 + (h*h)*(r*r)*1.05E2)*(1.0 / 1.05E2)) / r;
	
	Matrix[21] = S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*3.141592653589793*((S*S)*4.2E1 + (C*C)*(l*l)*1.3E1 + (h*h)*(r*r)*3.5E1 - C*S*l*3.5E1)*(2.0 / 3.5E1)) / (l*r)) + C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) + S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[23] = (S*h*k*3.141592653589793*((S*S)*l*4.2E1 - (C*C)*(l*l*l)*1.3E1 + (h*h)*(r*r)*2.1E2)*(1.0 / 2.1E2)) / (l*r) + (C*h*k*l*3.141592653589793*((S*S)*5.0 + C*S*l*2.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[25] = -S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) + C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	Matrix[27] = -S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) - C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	
	Matrix[36] = -C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*3.141592653589793*((S*S)*4.2E1 + (C*C)*(l*l)*1.3E1 + (h*h)*(r*r)*3.5E1 - C*S*l*3.5E1)*(2.0 / 3.5E1)) / (l*r)) - S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) + S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[38] = -C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*3.141592653589793*((S*S)*4.2E1 + (C*C)*(l*l)*1.3E1 + (h*h)*(r*r)*3.5E1 - C*S*l*3.5E1)*(2.0 / 3.5E1)) / (l*r)) + S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) + S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[40] = (C*h*k*3.141592653589793*((S*S)*l*4.2E1 - (C*C)*(l*l*l)*1.3E1 + (h*h)*(r*r)*2.1E2)*(-1.0 / 2.1E2)) / (l*r) + (S*h*k*l*3.141592653589793*((S*S)*5.0 + C*S*l*2.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[42] = C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) + S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	Matrix[44] = C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) - S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	
	Matrix[55] = (h*k*3.141592653589793*((C*C)*(l*l*l*l) + (S*S)*(l*l)*1.4E1 + (h*h)*(r*r)*3.5E1)*(2.0 / 1.05E2)) / (l*r);
	Matrix[57] = (C*h*k*3.141592653589793*((S*S)*l*2.1E1 + (C*C)*(l*l*l)*1.1E1 + (h*h)*(r*r)*1.05E2)*(-1.0 / 1.05E2)) / (l*r) - (S*h*k*l*3.141592653589793*((S*S)*5.0 + C*S*l*3.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[59] = (S*h*k*3.141592653589793*((S*S)*l*2.1E1 + (C*C)*(l*l*l)*1.1E1 + (h*h)*(r*r)*1.05E2)*(1.0 / 1.05E2)) / (l*r) - (C*h*k*l*3.141592653589793*((S*S)*5.0 + C*S*l*3.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[61] = (h*k*3.141592653589793*((C*C)*(l*l*l*l)*3.0 + (S*S)*(l*l)*1.4E1 - (h*h)*(r*r)*7.0E1)*(-1.0 / 2.1E2)) / (l*r);
	Matrix[63] = (C*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*-1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*2.1E2)*(1.0 / 2.1E2)) / r + (S*h*k*l*3.141592653589793*((S*S)*5.0 - C*S*l*2.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[65] = (C*h*k*l*3.141592653589793*((S*S)*5.0 - C*S*l*2.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r - (S*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*-1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*2.1E2)*(1.0 / 2.1E2)) / r;
}

//	Calculate element stress 
void CFrustum::ElementStress(double* stress, double* Displacement)
{

}
