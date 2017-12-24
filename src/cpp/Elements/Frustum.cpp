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
		for (unsigned int D = 0; D < 3; D++)
			LocationMatrix[i++] = nodes[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node frustum shell element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CFrustum::SizeOfStiffnessMatrix() { return 21; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CFrustum::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//	Only X and Z matter
	double DX[3];		//	dx = x2-x1, dy = 0, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];

	double L2 = DX[0] * DX[0] + DX[2] * DX[2];
	double l = sqrt(L2);
	double S = DX[2] / l;
	double C = DX[0] / l;
	double r = DX[2];
	//	Calculate element stiffness matrix

	CShellMaterial* material = static_cast<CShellMaterial*>(ElementMaterial);	// Pointer to material of the element

	const double nu = material->nu;
	const double k = material->E / (1 - nu*nu);
	const double h = material->h;

	Matrix[0] = S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*3.5E1 + C*S*(l*l*l)*3.5E1)*(2.0 / 3.5E1)) / r) + C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) - S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[1] = -C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*1.0 / (l*l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*3.5E1 + C*S*(l*l*l)*3.5E1)*(2.0 / 3.5E1)) / r) - S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) - S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[2] = -C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l*l)*3.141592653589793*((C*C)*(l*l*l*l)*1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*3.5E1 + C*S*(l*l*l)*3.5E1)*(2.0 / 3.5E1)) / r) + S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) - S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[3] = (h*k*3.141592653589793*((C*C)*(l*l*l*l)*2.0 + (S*S)*(l*l)*6.7E1 + (h*h)*(r*r)*7.0E1 + C*S*(l*l*l)*9.0)*(1.0 / 1.05E2)) / (l*r);
	Matrix[4] = -S*(h*k*3.141592653589793*(S*nu*1.5E1 + C*l*nu*5.0)*(1.0 / 3.0E1) - (h*k*3.141592653589793*((S*S)*l + C*S*(l*l)*3.0)*(1.0 / 3.0E1)) / r) + (C*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*2.2E1 + (S*S)*(l*l)*1.68E2 + (h*h)*(r*r)*2.1E2 + C*S*(l*l*l)*3.9E1)*(1.0 / 2.1E2)) / r;
	Matrix[5] = -C*(h*k*3.141592653589793*(S*nu*1.5E1 + C*l*nu*5.0)*(1.0 / 3.0E1) - (h*k*3.141592653589793*((S*S)*l + C*S*(l*l)*3.0)*(1.0 / 3.0E1)) / r) - (S*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*2.2E1 + (S*S)*(l*l)*1.68E2 + (h*h)*(r*r)*2.1E2 + C*S*(l*l*l)*3.9E1)*(1.0 / 2.1E2)) / r;
	Matrix[6] = S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*3.141592653589793*((S*S)*4.2E1 + (C*C)*(l*l)*1.3E1 + (h*h)*(r*r)*3.5E1 - C*S*l*3.5E1)*(2.0 / 3.5E1)) / (l*r)) + C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) + S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[7] = C*(h*k*3.141592653589793*(S*nu*1.5E1 + C*l*nu*5.0)*(1.0 / 3.0E1) + (h*k*3.141592653589793*((S*S)*l*1.4E1 + C*S*(l*l)*2.0)*(1.0 / 3.0E1)) / r) + (S*h*k*3.141592653589793*((S*S)*l*1.68E2 - (C*C)*(l*l*l)*1.3E1 + (h*h)*(r*r)*2.1E2 - C*S*(l*l)*6.6E1)*(1.0 / 2.1E2)) / (l*r);
	Matrix[8] = -S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) + C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	Matrix[9] = -S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) - C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	Matrix[10] = -C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*3.141592653589793*((S*S)*4.2E1 + (C*C)*(l*l)*1.3E1 + (h*h)*(r*r)*3.5E1 - C*S*l*3.5E1)*(2.0 / 3.5E1)) / (l*r)) - S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) + S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[11] = -C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*3.141592653589793*((S*S)*4.2E1 + (C*C)*(l*l)*1.3E1 + (h*h)*(r*r)*3.5E1 - C*S*l*3.5E1)*(2.0 / 3.5E1)) / (l*r)) + S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*7.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*3.0 + (S*S)*(l*l) + S*l*nu*r*3.0)*(2.0 / 3.0)) / (l*r));
	Matrix[12] = S*(h*k*3.141592653589793*(S*nu*1.5E1 + C*l*nu*5.0)*(1.0 / 3.0E1) + (h*k*3.141592653589793*((S*S)*l*1.4E1 + C*S*(l*l)*2.0)*(1.0 / 3.0E1)) / r) - (C*h*k*3.141592653589793*((S*S)*l*1.68E2 - (C*C)*(l*l*l)*1.3E1 + (h*h)*(r*r)*2.1E2 - C*S*(l*l)*6.6E1)*(1.0 / 2.1E2)) / (l*r);
	Matrix[13] = C*(S*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (C*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) + S*(C*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) - (S*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	Matrix[14] = C*(C*((h*k*3.141592653589793*(S*nu*2.0E1 - C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l - (h*k*3.141592653589793*((S*S)*l*1.0E1 - C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (S*h*k*1.0 / (l*l)*3.141592653589793*((S*S)*l*8.4E1 - (C*C)*(l*l*l)*9.0 + (h*h)*(r*r)*7.0E1)*(1.0 / 3.5E1)) / r) - S*(S*((h*k*3.141592653589793*(S*nu*2.0E1 + C*l*nu*1.0E1)*(1.0 / 1.0E1)) / l + (h*k*3.141592653589793*((S*S)*l*1.0E1 + C*S*(l*l)*3.0)*(1.0 / 1.0E1)) / (l*r)) + (C*h*k*3.141592653589793*((r*r)*6.0 - (S*S)*(l*l))*(1.0 / 3.0)) / (l*r));
	Matrix[15] = (h*k*3.141592653589793*((C*C)*(l*l*l*l) + (S*S)*(l*l)*1.4E1 + (h*h)*(r*r)*3.5E1)*(2.0 / 1.05E2)) / (l*r);
	Matrix[16] = (C*h*k*3.141592653589793*((S*S)*l*2.1E1 + (C*C)*(l*l*l)*1.1E1 + (h*h)*(r*r)*1.05E2)*(-1.0 / 1.05E2)) / (l*r) - (S*h*k*l*3.141592653589793*((S*S)*5.0 + C*S*l*3.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[17] = (S*h*k*3.141592653589793*((S*S)*l*2.1E1 + (C*C)*(l*l*l)*1.1E1 + (h*h)*(r*r)*1.05E2)*(1.0 / 1.05E2)) / (l*r) - (C*h*k*l*3.141592653589793*((S*S)*5.0 + C*S*l*3.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[18] = (h*k*3.141592653589793*((C*C)*(l*l*l*l)*3.0 + (S*S)*(l*l)*1.4E1 - (h*h)*(r*r)*7.0E1 + C*S*(l*l*l)*1.2E1)*(-1.0 / 2.1E2)) / (l*r);
	Matrix[19] = (C*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*-1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*2.1E2)*(1.0 / 2.1E2)) / r + (S*h*k*l*3.141592653589793*((S*S)*5.0 - C*S*l*2.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r;
	Matrix[20] = (C*h*k*l*3.141592653589793*((S*S)*5.0 - C*S*l*2.0 + C*nu*r*5.0)*(1.0 / 3.0E1)) / r - (S*h*k*1.0 / (l*l)*3.141592653589793*((C*C)*(l*l*l*l)*-1.3E1 + (S*S)*(l*l)*4.2E1 + (h*h)*(r*r)*2.1E2)*(1.0 / 2.1E2)) / r;
}

//	Calculate element stress 
void CFrustum::ElementStress(double* stress, double* Displacement)
{

}
