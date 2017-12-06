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

	ElementMaterial = NULL;
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
	ElementMaterial = dynamic_cast<CTimoshenkoMaterial*>(MaterialSets) + MSet - 1;
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

	// Calculate beam length
	// DX[3] = [x2-x1, y2-y1, z2-z1]
	double DX[3];		
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
	double len = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

	// Get the Material/Section property
	// Pointer to material of the element
	CTimoshenkoMaterial* material = dynamic_cast<CTimoshenkoMaterial*>(ElementMaterial);
	
	// Calculate the transform matrix
	// Q = [Q11,Q21,Q31,Q12,Q22,Q32,Q13,Q23,Q33],
	//		in which, QIJ = eI \cdot eJ'
	double Q[9];		
	Q[0] = DX[0] / len;
	Q[1] = DX[1] / len;
	Q[2] = DX[2] / len;
	Q[3] = material->Thetay1;
	Q[4] = material->Thetay2;
	Q[5] = material->Thetay3;
	Q[6] = Q[1] * Q[5] - Q[2] * Q[4];
	Q[7] = Q[2] * Q[3] - Q[0] * Q[5];
	Q[8] = Q[0] * Q[4] - Q[1] * Q[3];

	// Calculate the elements of element stiffness matrix in local coordinates
	/*! by & bz are ONLY VALID FOR beams with rectangular section !*/
	const double kMod = 6.0 / 5.0;
	const double G = material->E / 2 / (1 + material->nu);	// Shear Modulus
	const double J = material->Iyy + material->Izz; // Polar momentum of inertia
	const double Tens = material->E * material->Area / len;
	const double Tors = G * J / len;
	const double Stfy = material->E * material->Iyy;
	const double Stfz = material->E * material->Izz;
	const double Stfs = G * material->Area / kMod;
#ifndef _TIMOSHENKO_ACCURATE_INTEGRATION_
	const double Bendy1 = Stfs / len;
	const double Bendy2 = Stfs * 0.5;
	const double Bendy3 = Stfs / 4 * len + Stfy / len;
	const double Bendy4 = Stfs / 4 * len - Stfy / len;
	const double Bendz1 = Stfs / len;
	const double Bendz2 = Stfs * 0.5;
	const double Bendz3 = Stfs / 4 * len + Stfz / len;
	const double Bendz4 = Stfs / 4 * len - Stfz / len;
#else
	const double Bendy1 = Stfs / len;
	const double Bendy2 = Stfs * 0.5;
	const double Bendy3 = Stfs / 3 * len + Stfy / len;
	const double Bendy4 = Stfs / 6 * len - Stfy / len;
	const double Bendz1 = Stfs / len;
	const double Bendz2 = Stfs * 0.5;
	const double Bendz3 = Stfs / 3 * len + Stfz / len;
	const double Bendz4 = Stfs / 6 * len - Stfz / len;
#endif
	
	// Calculate the elements of element stiffness matrix in background coordinates
  Matrix[0] = Bendz1*(Q[3]*Q[3])+Bendy1*(Q[6]*Q[6])+(Q[0]*Q[0])*Tens;
  Matrix[1] = Bendz1*(Q[4]*Q[4])+Bendy1*(Q[7]*Q[7])+(Q[1]*Q[1])*Tens;
  Matrix[2] = Bendz1*Q[3]*Q[4]+Bendy1*Q[6]*Q[7]+Q[0]*Q[1]*Tens;
  Matrix[3] = Bendz1*(Q[5]*Q[5])+Bendy1*(Q[8]*Q[8])+(Q[2]*Q[2])*Tens;
  Matrix[4] = Bendz1*Q[4]*Q[5]+Bendy1*Q[7]*Q[8]+Q[1]*Q[2]*Tens;
  Matrix[5] = Bendz1*Q[3]*Q[5]+Bendy1*Q[6]*Q[8]+Q[0]*Q[2]*Tens;
  Matrix[6] = Bendy3*(Q[3]*Q[3])+Bendz3*(Q[6]*Q[6])+(Q[0]*Q[0])*Tors;
  Matrix[7] = -Bendy2*Q[3]*Q[8]+Bendz2*Q[5]*Q[6];
  Matrix[8] = -Bendy2*Q[3]*Q[7]+Bendz2*Q[4]*Q[6];
  Matrix[9] = -Bendy2*Q[3]*Q[6]+Bendz2*Q[3]*Q[6];
  Matrix[10] = Bendy3*(Q[4]*Q[4])+Bendz3*(Q[7]*Q[7])+(Q[1]*Q[1])*Tors;
  Matrix[11] = Bendy3*Q[3]*Q[4]+Bendz3*Q[6]*Q[7]+Q[0]*Q[1]*Tors;
  Matrix[12] = -Bendy2*Q[4]*Q[8]+Bendz2*Q[5]*Q[7];
  Matrix[13] = -Bendy2*Q[4]*Q[7]+Bendz2*Q[4]*Q[7];
  Matrix[14] = -Bendy2*Q[4]*Q[6]+Bendz2*Q[3]*Q[7];
  Matrix[15] = Bendy3*(Q[5]*Q[5])+Bendz3*(Q[8]*Q[8])+(Q[2]*Q[2])*Tors;
  Matrix[16] = Bendy3*Q[4]*Q[5]+Bendz3*Q[7]*Q[8]+Q[1]*Q[2]*Tors;
  Matrix[17] = Bendy3*Q[3]*Q[5]+Bendz3*Q[6]*Q[8]+Q[0]*Q[2]*Tors;
  Matrix[18] = -Bendy2*Q[5]*Q[8]+Bendz2*Q[5]*Q[8];
  Matrix[19] = -Bendy2*Q[5]*Q[7]+Bendz2*Q[4]*Q[8];
  Matrix[20] = -Bendy2*Q[5]*Q[6]+Bendz2*Q[3]*Q[8];
  Matrix[21] = Bendz1*(Q[3]*Q[3])+Bendy1*(Q[6]*Q[6])+(Q[0]*Q[0])*Tens;
  Matrix[22] = Bendy2*Q[5]*Q[6]-Bendz2*Q[3]*Q[8];
  Matrix[23] = Bendy2*Q[4]*Q[6]-Bendz2*Q[3]*Q[7];
  Matrix[24] = Bendy2*Q[3]*Q[6]-Bendz2*Q[3]*Q[6];
  Matrix[25] = -Bendz1*Q[3]*Q[5]-Bendy1*Q[6]*Q[8]-Q[0]*Q[2]*Tens;
  Matrix[26] = -Bendz1*Q[3]*Q[4]-Bendy1*Q[6]*Q[7]-Q[0]*Q[1]*Tens;
  Matrix[27] = -Bendz1*(Q[3]*Q[3])-Bendy1*(Q[6]*Q[6])-(Q[0]*Q[0])*Tens;
  Matrix[28] = Bendz1*(Q[4]*Q[4])+Bendy1*(Q[7]*Q[7])+(Q[1]*Q[1])*Tens;
  Matrix[29] = Bendz1*Q[3]*Q[4]+Bendy1*Q[6]*Q[7]+Q[0]*Q[1]*Tens;
  Matrix[30] = Bendy2*Q[5]*Q[7]-Bendz2*Q[4]*Q[8];
  Matrix[31] = Bendy2*Q[4]*Q[7]-Bendz2*Q[4]*Q[7];
  Matrix[32] = Bendy2*Q[3]*Q[7]-Bendz2*Q[4]*Q[6];
  Matrix[33] = -Bendz1*Q[4]*Q[5]-Bendy1*Q[7]*Q[8]-Q[1]*Q[2]*Tens;
  Matrix[34] = -Bendz1*(Q[4]*Q[4])-Bendy1*(Q[7]*Q[7])-(Q[1]*Q[1])*Tens;
  Matrix[35] = -Bendz1*Q[3]*Q[4]-Bendy1*Q[6]*Q[7]-Q[0]*Q[1]*Tens;
  Matrix[36] = Bendz1*(Q[5]*Q[5])+Bendy1*(Q[8]*Q[8])+(Q[2]*Q[2])*Tens;
  Matrix[37] = Bendz1*Q[4]*Q[5]+Bendy1*Q[7]*Q[8]+Q[1]*Q[2]*Tens;
  Matrix[38] = Bendz1*Q[3]*Q[5]+Bendy1*Q[6]*Q[8]+Q[0]*Q[2]*Tens;
  Matrix[39] = Bendy2*Q[5]*Q[8]-Bendz2*Q[5]*Q[8];
  Matrix[40] = Bendy2*Q[4]*Q[8]-Bendz2*Q[5]*Q[7];
  Matrix[41] = Bendy2*Q[3]*Q[8]-Bendz2*Q[5]*Q[6];
  Matrix[42] = -Bendz1*(Q[5]*Q[5])-Bendy1*(Q[8]*Q[8])-(Q[2]*Q[2])*Tens;
  Matrix[43] = -Bendz1*Q[4]*Q[5]-Bendy1*Q[7]*Q[8]-Q[1]*Q[2]*Tens;
  Matrix[44] = -Bendz1*Q[3]*Q[5]-Bendy1*Q[6]*Q[8]-Q[0]*Q[2]*Tens;
  Matrix[45] = Bendy3*(Q[3]*Q[3])+Bendz3*(Q[6]*Q[6])+(Q[0]*Q[0])*Tors;
  Matrix[46] = Bendy2*Q[3]*Q[8]-Bendz2*Q[5]*Q[6];
  Matrix[47] = Bendy2*Q[3]*Q[7]-Bendz2*Q[4]*Q[6];
  Matrix[48] = Bendy2*Q[3]*Q[6]-Bendz2*Q[3]*Q[6];
  Matrix[49] = Bendy4*Q[3]*Q[5]+Bendz4*Q[6]*Q[8]-Q[0]*Q[2]*Tors;
  Matrix[50] = Bendy4*Q[3]*Q[4]+Bendz4*Q[6]*Q[7]-Q[0]*Q[1]*Tors;
  Matrix[51] = Bendy4*(Q[3]*Q[3])+Bendz4*(Q[6]*Q[6])-(Q[0]*Q[0])*Tors;
  Matrix[52] = -Bendy2*Q[3]*Q[8]+Bendz2*Q[5]*Q[6];
  Matrix[53] = -Bendy2*Q[3]*Q[7]+Bendz2*Q[4]*Q[6];
  Matrix[54] = -Bendy2*Q[3]*Q[6]+Bendz2*Q[3]*Q[6];
  Matrix[55] = Bendy3*(Q[4]*Q[4])+Bendz3*(Q[7]*Q[7])+(Q[1]*Q[1])*Tors;
  Matrix[56] = Bendy3*Q[3]*Q[4]+Bendz3*Q[6]*Q[7]+Q[0]*Q[1]*Tors;
  Matrix[57] = Bendy2*Q[4]*Q[8]-Bendz2*Q[5]*Q[7];
  Matrix[58] = Bendy2*Q[4]*Q[7]-Bendz2*Q[4]*Q[7];
  Matrix[59] = Bendy2*Q[4]*Q[6]-Bendz2*Q[3]*Q[7];
  Matrix[60] = Bendy4*Q[4]*Q[5]+Bendz4*Q[7]*Q[8]-Q[1]*Q[2]*Tors;
  Matrix[61] = Bendy4*(Q[4]*Q[4])+Bendz4*(Q[7]*Q[7])-(Q[1]*Q[1])*Tors;
  Matrix[62] = Bendy4*Q[3]*Q[4]+Bendz4*Q[6]*Q[7]-Q[0]*Q[1]*Tors;
  Matrix[63] = -Bendy2*Q[4]*Q[8]+Bendz2*Q[5]*Q[7];
  Matrix[64] = -Bendy2*Q[4]*Q[7]+Bendz2*Q[4]*Q[7];
  Matrix[65] = -Bendy2*Q[4]*Q[6]+Bendz2*Q[3]*Q[7];
  Matrix[66] = Bendy3*(Q[5]*Q[5])+Bendz3*(Q[8]*Q[8])+(Q[2]*Q[2])*Tors;
  Matrix[67] = Bendy3*Q[4]*Q[5]+Bendz3*Q[7]*Q[8]+Q[1]*Q[2]*Tors;
  Matrix[68] = Bendy3*Q[3]*Q[5]+Bendz3*Q[6]*Q[8]+Q[0]*Q[2]*Tors;
  Matrix[69] = Bendy2*Q[5]*Q[8]-Bendz2*Q[5]*Q[8];
  Matrix[70] = Bendy2*Q[5]*Q[7]-Bendz2*Q[4]*Q[8];
  Matrix[71] = Bendy2*Q[5]*Q[6]-Bendz2*Q[3]*Q[8];
  Matrix[72] = Bendy4*(Q[5]*Q[5])+Bendz4*(Q[8]*Q[8])-(Q[2]*Q[2])*Tors;
  Matrix[73] = Bendy4*Q[4]*Q[5]+Bendz4*Q[7]*Q[8]-Q[1]*Q[2]*Tors;
  Matrix[74] = Bendy4*Q[3]*Q[5]+Bendz4*Q[6]*Q[8]-Q[0]*Q[2]*Tors;
  Matrix[75] = -Bendy2*Q[5]*Q[8]+Bendz2*Q[5]*Q[8];
  Matrix[76] = -Bendy2*Q[5]*Q[7]+Bendz2*Q[4]*Q[8];
  Matrix[77] = -Bendy2*Q[5]*Q[6]+Bendz2*Q[3]*Q[8];

/* ---------- JUST FOR DEBUG ---------- */
/*
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++){
			if (j < i) {
				cout << "            ";
			}
			else {
				cout << Matrix[j*(j + 1) / 2 + j - i] << "  ";
			}
		}
		cout << endl << endl;
	}
*/
}

// Calculate element stress
// stress               double[3] = {sigmaXX, sigmaXY, sigmaXZ};
// force                double[12] = {Fx1, Fx2, Fy1, Fy2, Fz1, Fz2, Mx1, Mx2, My1, My2, Mz1, Mz2}
// Displacement:        double[NEQ], represent the nodal displacement.
// Positions:           double[12], represent 3d position for 3 gauss points.
// GaussDisplacements:  double[12], represent 3d displacements for 3 gauss points.
// Weights:             double[2], represent integrate weights.
void CTimoshenkoSRINT::ElementStress(double stress[3], double force[12], double* Displacement)
{
	// Get the Material/Section property
	// Pointer to material of the element
	CTimoshenkoMaterial* material = dynamic_cast<CTimoshenkoMaterial*>(ElementMaterial);

	// Calculate beam length
	// DX[3] = {x2-x1, y2-y1, z2-z1}
	// X[12] = {x1,y1,z1,...,x2,y2,z2,...}
	double DX[3];
	for (unsigned int i = 0; i < 3; i++) {
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
	}
	double len = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

	// Calculate the transform matrix
	// Q = [Q11,Q21,Q31,Q12,Q22,Q32,Q13,Q23,Q33],
	//		in which, QIJ = eI \cdot eJ'
	double Q[9];		
	Q[0] = DX[0] / len;
	Q[1] = DX[1] / len;
	Q[2] = DX[2] / len;
	Q[3] = material->Thetay1;
	Q[4] = material->Thetay2;
	Q[5] = material->Thetay3;
	Q[6] = Q[1] * Q[5] - Q[2] * Q[4];
	Q[7] = Q[2] * Q[3] - Q[0] * Q[5];
	Q[8] = Q[0] * Q[4] - Q[1] * Q[3];

	// Tranfer ele-background disp "displacement" to ele-local disp "Eledisp"
	// Disp[12] = dispX1,dispY1,dispZ1,...,dispX2,dispY2,dispZ2,...}
	double Disp[12];
	for (unsigned index = 0; index < 12; ++index)
	{
		if (LocationMatrix[index])
		{
			Disp[index] = Displacement[LocationMatrix[index] - 1];
		}
		else
		{
			Disp[index] = 0.0;
		}
	}
	double EleDisp[12];
	EleDisp[0] = Q[0] * Disp[0] + Q[1] * Disp[1] + Q[2] * Disp[2];
	EleDisp[1] = Q[3] * Disp[0] + Q[4] * Disp[1] + Q[5] * Disp[2];
	EleDisp[2] = Q[6] * Disp[0] + Q[7] * Disp[1] + Q[8] * Disp[2];
	EleDisp[3] = Q[0] * Disp[3] + Q[1] * Disp[4] + Q[2] * Disp[5];
	EleDisp[4] = Q[3] * Disp[3] + Q[4] * Disp[4] + Q[5] * Disp[5];
	EleDisp[5] = Q[6] * Disp[3] + Q[7] * Disp[4] + Q[8] * Disp[5];
	EleDisp[6] = Q[0] * Disp[6] + Q[1] * Disp[7] + Q[2] * Disp[8];
	EleDisp[7] = Q[3] * Disp[6] + Q[4] * Disp[7] + Q[5] * Disp[8];
	EleDisp[8] = Q[6] * Disp[6] + Q[7] * Disp[7] + Q[8] * Disp[8];
	EleDisp[9] = Q[0] * Disp[9] + Q[1] * Disp[10] + Q[2] * Disp[11];
	EleDisp[10] = Q[3] * Disp[9] + Q[4] * Disp[10] + Q[5] * Disp[11];
	EleDisp[11] = Q[6] * Disp[9] + Q[7] * Disp[10] + Q[8] * Disp[11];

	// Calculate the elements of element stiffness matrix in local coordinates
	/*! by & bz are ONLY VALID FOR beams with rectangular section !*/
	const double kMod = 6.0 / 5.0;
	const double G = material->E / 2 / (1 + material->nu);	// Shear Modulus
	const double J = material->Iyy + material->Izz;			// Polar Momentum of Inertia
	const double Tens = material->E * material->Area / len;
	const double Tors = G * J / len;
	const double Stfy = material->E * material->Iyy;
	const double Stfz = material->E * material->Izz;
	const double Stfs = G * material->Area / kMod;

	// Calculate the stress
	stress[0] = material->E * (EleDisp[6] - EleDisp[0]) / len;
	stress[1] = G / kMod / len * ((EleDisp[7] - EleDisp[1])-0.5*(EleDisp[5] + EleDisp[11])*len);
	stress[2] = G / kMod / len * ((EleDisp[8] - EleDisp[2])+0.5*(EleDisp[4] + EleDisp[10])*len);
	// Calculate the internal forces on nodal points
	/* ONLY use the shear force on the midpoint of the element*/
	force[1] = Tens * (EleDisp[6] - EleDisp[0]);
	force[0] = - force[1];
	force[3] = Stfs / len * (EleDisp[7] - EleDisp[1] - 0.5 * (EleDisp[5] + EleDisp[11]) * len);
	force[2] = - force[3];
	force[5] = Stfs / len * (EleDisp[8] - EleDisp[2] + 0.5 * (EleDisp[4]  + EleDisp[10])* len);
	force[4] = - force[5];
	force[6] = Tors * (EleDisp[9] - EleDisp[3]);
	force[7] = - force[6];
	force[9] = Stfy / len * (EleDisp[10] - EleDisp[4]);
	force[8] = - force[9] - force[6] * len;
	force[11] = Stfz / len * (EleDisp[11] - EleDisp[5]);
	force[10] = - force[11] + force[4] * len;
}
