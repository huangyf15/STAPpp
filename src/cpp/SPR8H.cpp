#include "Elements/8H.h"
#include "Eigen/dense"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;
using namespace std;


void  CHex::ElementPostSPR(double* stressG, double* Displacement , double* PrePositions, double* PostPositions, double* PositionG, unsigned int* mingtianzaishuo)
{
	// get original position: preposition
	for (unsigned int i =0 ; i<3; i++)
	{
		for (unsigned int j=0;j < 8; j++)
		{
			PrePositions[i+3*j] = nodes[j]->XYZ[i];	 			
		}
	}

	// Get nodal displacements [LM can be used here]
	for (unsigned int i = 0; i < 24; i++)
	{
		if (LocationMatrix[i])
			//locatiion matrix start from 1 not 0
			PostPositions[i] = Displacement[LocationMatrix[i]-1];
		else
			PostPositions[i] = 0.0;
	}

	// Construct constitutive matrix
	CHexMaterial* material = static_cast<CHexMaterial*>(ElementMaterial);	// Pointer to material of the element
	double v = material->nu;
	double k = material->E * (1-v)/(1+v)/(1-2*v);
	double D[3];
	D[0] = k;
	D[1] = k * v / (1 - v);
	D[2] = k * (1 - 2 * v) / 2.0 / (1 - v);


	// Construct Jacobi matrix
	const double xi8[8] = { 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 };
	const double eta8[8] = { -0.577350269189626 , 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,0.577350269189626 , 0.577350269189626 ,-0.577350269189626 };
	const double zeta8[8] = { -0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 , 0.577350269189626 , 0.577350269189626 };

	double stressXYZ[6][8];	// 8 gauss points, 6 stress components
	for (unsigned p = 0; p < 8; p++)
	{
		double xi   = xi8[p];
		double eta  = eta8[p];
		double zeta = zeta8[p];

		double GN[12];
		GN[0] = (1-eta)*(1-zeta) / 8.0;
		GN[1] = (1+eta)*(1-zeta) / 8.0;
		GN[2] = (1-eta)*(1+zeta) / 8.0;
		GN[3] = (1+eta)*(1+zeta) / 8.0;
		GN[4] = (1+xi)*(1-zeta) / 8.0;
		GN[5] = (1-xi)*(1-zeta) / 8.0;
		GN[6] = (1+xi)*(1+zeta) / 8.0;
		GN[7] = (1-xi)*(1+zeta) / 8.0;
		GN[8] = (1+xi)*(1-eta) / 8.0;
		GN[9] = (1+xi)*(1+eta) / 8.0;
		GN[10] = (1-xi)*(1+eta) / 8.0;
		GN[11] = (1-xi)*(1-eta) / 8.0;

		double J[9];
		J[0] = PrePositions[0] * GN[0] + PrePositions[3] * GN[1] - PrePositions[6] * GN[1] - PrePositions[9] * GN[0] + PrePositions[12] * GN[2] + PrePositions[15] * GN[3] - PrePositions[18] * GN[3] - PrePositions[21] * GN[2];
		J[1] = -PrePositions[0] * GN[4] + PrePositions[3] * GN[4] + PrePositions[6] * GN[5] - PrePositions[9] * GN[5] - PrePositions[12] * GN[6] + PrePositions[15] * GN[6] + PrePositions[18] * GN[7] - PrePositions[21] * GN[7];
		J[2] = -PrePositions[0] * GN[8] - PrePositions[3] * GN[9] - PrePositions[6] * GN[10] - PrePositions[9] * GN[11] + PrePositions[12] * GN[8] + PrePositions[15] * GN[9] + PrePositions[18] * GN[10] + PrePositions[21] * GN[11];
		J[3] = PrePositions[1] * GN[0] + PrePositions[4] * GN[1] - PrePositions[7] * GN[1] - PrePositions[10] * GN[0] + PrePositions[13] * GN[2] + PrePositions[16] * GN[3] - PrePositions[19] * GN[3] - PrePositions[22] * GN[2];
		J[4] = -PrePositions[1] * GN[4] + PrePositions[4] * GN[4] + PrePositions[7] * GN[5] - PrePositions[10] * GN[5] - PrePositions[13] * GN[6] + PrePositions[16] * GN[6] + PrePositions[19] * GN[7] - PrePositions[22] * GN[7];
		J[5] = -PrePositions[1] * GN[8] - PrePositions[4] * GN[9] - PrePositions[7] * GN[10] - PrePositions[10] * GN[11] + PrePositions[13] * GN[8] + PrePositions[16] * GN[9] + PrePositions[19] * GN[10] + PrePositions[22] * GN[11];
		J[6] = PrePositions[2] * GN[0] + PrePositions[5] * GN[1] - PrePositions[8] * GN[1] - PrePositions[11] * GN[0] + PrePositions[14] * GN[2] + PrePositions[17] * GN[3] - PrePositions[20] * GN[3] - PrePositions[23] * GN[2];
		J[7] = -PrePositions[2] * GN[4] + PrePositions[5] * GN[4] + PrePositions[8] * GN[5] - PrePositions[11] * GN[5] - PrePositions[14] * GN[6] + PrePositions[17] * GN[6] + PrePositions[20] * GN[7] - PrePositions[23] * GN[7];
		J[8] = -PrePositions[2] * GN[8] - PrePositions[5] * GN[9] - PrePositions[8] * GN[10] - PrePositions[11] * GN[11] + PrePositions[14] * GN[8] + PrePositions[17] * GN[9] + PrePositions[20] * GN[10] + PrePositions[23] * GN[11];

		double detJ = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];

		double invJ[9];
		invJ[0] = (J[4]*J[8]-J[5]*J[7])/detJ;
		invJ[1] = -(J[1]*J[8]-J[2]*J[7])/detJ;
		invJ[2] = (J[1]*J[5]-J[2]*J[4])/detJ;
		invJ[3] = -(J[3]*J[8]-J[5]*J[6])/detJ;
		invJ[4] = (J[0]*J[8]-J[2]*J[6])/detJ;
		invJ[5] = -(J[0]*J[5]-J[2]*J[3])/detJ;
		invJ[6] = (J[3]*J[7]-J[4]*J[6])/detJ;
		invJ[7] = -(J[0]*J[7]-J[1]*J[6])/detJ;
		invJ[8] = (J[0]*J[4]-J[1]*J[3])/detJ;

		double kerB[24];
		kerB[0] = GN[0] * invJ[0] - GN[4] * invJ[3] - GN[8] * invJ[6];
		kerB[1] = GN[0] * invJ[1] - GN[4] * invJ[4] - GN[8] * invJ[7];
		kerB[2] = GN[0] * invJ[2] - GN[4] * invJ[5] - GN[8] * invJ[8];
		kerB[3] = GN[1] * invJ[0] + GN[4] * invJ[3] - GN[9] * invJ[6];
		kerB[4] = GN[1] * invJ[1] + GN[4] * invJ[4] - GN[9] * invJ[7];
		kerB[5] = GN[1] * invJ[2] + GN[4] * invJ[5] - GN[9] * invJ[8];
		kerB[6] = -GN[1] * invJ[0] + GN[5] * invJ[3] - GN[10] * invJ[6];
		kerB[7] = -GN[1] * invJ[1] + GN[5] * invJ[4] - GN[10] * invJ[7];
		kerB[8] = -GN[1] * invJ[2] + GN[5] * invJ[5] - GN[10] * invJ[8];
		kerB[9] = -GN[0] * invJ[0] - GN[5] * invJ[3] - GN[11] * invJ[6];
		kerB[10] = -GN[0] * invJ[1] - GN[5] * invJ[4] - GN[11] * invJ[7];
		kerB[11] = -GN[0] * invJ[2] - GN[5] * invJ[5] - GN[11] * invJ[8];
		kerB[12] = GN[2] * invJ[0] - GN[6] * invJ[3] + GN[8] * invJ[6];
		kerB[13] = GN[2] * invJ[1] - GN[6] * invJ[4] + GN[8] * invJ[7];
		kerB[14] = GN[2] * invJ[2] - GN[6] * invJ[5] + GN[8] * invJ[8];
		kerB[15] = GN[3] * invJ[0] + GN[6] * invJ[3] + GN[9] * invJ[6];
		kerB[16] = GN[3] * invJ[1] + GN[6] * invJ[4] + GN[9] * invJ[7];
		kerB[17] = GN[3] * invJ[2] + GN[6] * invJ[5] + GN[9] * invJ[8];
		kerB[18] = -GN[3] * invJ[0] + GN[7] * invJ[3] + GN[10] * invJ[6];
		kerB[19] = -GN[3] * invJ[1] + GN[7] * invJ[4] + GN[10] * invJ[7];
		kerB[20] = -GN[3] * invJ[2] + GN[7] * invJ[5] + GN[10] * invJ[8];
		kerB[21] = -GN[2] * invJ[0] - GN[7] * invJ[3] + GN[11] * invJ[6];
		kerB[22] = -GN[2] * invJ[1] - GN[7] * invJ[4] + GN[11] * invJ[7];
		kerB[23] = -GN[2] * invJ[2] - GN[7] * invJ[5] + GN[11] * invJ[8];

		stressXYZ[0][p] = D[0] * PostPositions[0] * kerB[0] + D[1] * PostPositions[1] * kerB[1] + D[1] * PostPositions[2] * kerB[2] + D[0] * PostPositions[3] * kerB[3] + D[1] * PostPositions[4] * kerB[4] + D[1] * PostPositions[5] * kerB[5] + D[0] * PostPositions[6] * kerB[6] + D[1] * PostPositions[7] * kerB[7] + D[1] * PostPositions[8] * kerB[8] + D[0] * PostPositions[9] * kerB[9] + D[1] * PostPositions[10] * kerB[10] + D[1] * PostPositions[11] * kerB[11] + D[0] * PostPositions[12] * kerB[12] + D[1] * PostPositions[13] * kerB[13] + D[1] * PostPositions[14] * kerB[14] + D[0] * PostPositions[15] * kerB[15] + D[1] * PostPositions[16] * kerB[16] + D[1] * PostPositions[17] * kerB[17] + D[0] * PostPositions[18] * kerB[18] + D[1] * PostPositions[19] * kerB[19] + D[1] * PostPositions[20] * kerB[20] + D[0] * PostPositions[21] * kerB[21] + D[1] * PostPositions[22] * kerB[22] + D[1] * PostPositions[23] * kerB[23];
		stressXYZ[1][p] = D[1] * PostPositions[0] * kerB[0] + D[0] * PostPositions[1] * kerB[1] + D[1] * PostPositions[2] * kerB[2] + D[1] * PostPositions[3] * kerB[3] + D[0] * PostPositions[4] * kerB[4] + D[1] * PostPositions[5] * kerB[5] + D[1] * PostPositions[6] * kerB[6] + D[0] * PostPositions[7] * kerB[7] + D[1] * PostPositions[8] * kerB[8] + D[1] * PostPositions[9] * kerB[9] + D[0] * PostPositions[10] * kerB[10] + D[1] * PostPositions[11] * kerB[11] + D[1] * PostPositions[12] * kerB[12] + D[0] * PostPositions[13] * kerB[13] + D[1] * PostPositions[14] * kerB[14] + D[1] * PostPositions[15] * kerB[15] + D[0] * PostPositions[16] * kerB[16] + D[1] * PostPositions[17] * kerB[17] + D[1] * PostPositions[18] * kerB[18] + D[0] * PostPositions[19] * kerB[19] + D[1] * PostPositions[20] * kerB[20] + D[1] * PostPositions[21] * kerB[21] + D[0] * PostPositions[22] * kerB[22] + D[1] * PostPositions[23] * kerB[23];
		stressXYZ[2][p] = D[1] * PostPositions[0] * kerB[0] + D[1] * PostPositions[1] * kerB[1] + D[0] * PostPositions[2] * kerB[2] + D[1] * PostPositions[3] * kerB[3] + D[1] * PostPositions[4] * kerB[4] + D[0] * PostPositions[5] * kerB[5] + D[1] * PostPositions[6] * kerB[6] + D[1] * PostPositions[7] * kerB[7] + D[0] * PostPositions[8] * kerB[8] + D[1] * PostPositions[9] * kerB[9] + D[1] * PostPositions[10] * kerB[10] + D[0] * PostPositions[11] * kerB[11] + D[1] * PostPositions[12] * kerB[12] + D[1] * PostPositions[13] * kerB[13] + D[0] * PostPositions[14] * kerB[14] + D[1] * PostPositions[15] * kerB[15] + D[1] * PostPositions[16] * kerB[16] + D[0] * PostPositions[17] * kerB[17] + D[1] * PostPositions[18] * kerB[18] + D[1] * PostPositions[19] * kerB[19] + D[0] * PostPositions[20] * kerB[20] + D[1] * PostPositions[21] * kerB[21] + D[1] * PostPositions[22] * kerB[22] + D[0] * PostPositions[23] * kerB[23];
		stressXYZ[3][p] = D[2] * PostPositions[0] * kerB[1] + D[2] * PostPositions[1] * kerB[0] + D[2] * PostPositions[3] * kerB[4] + D[2] * PostPositions[4] * kerB[3] + D[2] * PostPositions[6] * kerB[7] + D[2] * PostPositions[7] * kerB[6] + D[2] * PostPositions[9] * kerB[10] + D[2] * PostPositions[10] * kerB[9] + D[2] * PostPositions[12] * kerB[13] + D[2] * PostPositions[13] * kerB[12] + D[2] * PostPositions[15] * kerB[16] + D[2] * PostPositions[16] * kerB[15] + D[2] * PostPositions[18] * kerB[19] + D[2] * PostPositions[19] * kerB[18] + D[2] * PostPositions[21] * kerB[22] + D[2] * PostPositions[22] * kerB[21];
		stressXYZ[4][p] = D[2] * PostPositions[1] * kerB[2] + D[2] * PostPositions[2] * kerB[1] + D[2] * PostPositions[4] * kerB[5] + D[2] * PostPositions[5] * kerB[4] + D[2] * PostPositions[7] * kerB[8] + D[2] * PostPositions[8] * kerB[7] + D[2] * PostPositions[10] * kerB[11] + D[2] * PostPositions[11] * kerB[10] + D[2] * PostPositions[13] * kerB[14] + D[2] * PostPositions[14] * kerB[13] + D[2] * PostPositions[16] * kerB[17] + D[2] * PostPositions[17] * kerB[16] + D[2] * PostPositions[19] * kerB[20] + D[2] * PostPositions[20] * kerB[19] + D[2] * PostPositions[22] * kerB[23] + D[2] * PostPositions[23] * kerB[22];
		stressXYZ[5][p] = D[2] * PostPositions[0] * kerB[2] + D[2] * PostPositions[2] * kerB[0] + D[2] * PostPositions[3] * kerB[5] + D[2] * PostPositions[5] * kerB[3] + D[2] * PostPositions[6] * kerB[8] + D[2] * PostPositions[8] * kerB[6] + D[2] * PostPositions[9] * kerB[11] + D[2] * PostPositions[11] * kerB[9] + D[2] * PostPositions[12] * kerB[14] + D[2] * PostPositions[14] * kerB[12] + D[2] * PostPositions[15] * kerB[17] + D[2] * PostPositions[17] * kerB[15] + D[2] * PostPositions[18] * kerB[20] + D[2] * PostPositions[20] * kerB[18] + D[2] * PostPositions[21] * kerB[23] + D[2] * PostPositions[23] * kerB[21];
	}

	// stress recovery for stress on nodes
	double interpo[4] = {2.549038105676658, -0.683012701892219, 0.183012701892219, -0.049038105676658};
	double recovery[8];
	for (unsigned i = 0; i < 6; i++)
	{
		recovery[0] = interpo[0]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][3] + interpo[2]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][7] + interpo[3]*stressXYZ[i][6];
		recovery[1] = interpo[0]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][6] + interpo[3]*stressXYZ[i][7];
		recovery[2] = interpo[0]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][5] + interpo[3]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][7];
		recovery[3] = interpo[1]*stressXYZ[i][0] + interpo[0]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][4] + interpo[1]*stressXYZ[i][7] + interpo[2]*stressXYZ[i][6] + interpo[3]*stressXYZ[i][5];
		recovery[4] = interpo[1]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][1] + interpo[0]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][3] + interpo[3]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][5] + interpo[1]*stressXYZ[i][7] + interpo[2]*stressXYZ[i][6];
		recovery[5] = interpo[1]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][2] + interpo[0]*stressXYZ[i][5] + interpo[1]*stressXYZ[i][4] + interpo[3]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][7];
		recovery[6] = interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][1] + interpo[3]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][3] + interpo[0]*stressXYZ[i][6] + interpo[1]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][4] + interpo[1]*stressXYZ[i][7];
		recovery[7] = interpo[2]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][3] + interpo[2]*stressXYZ[i][2] + interpo[3]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][4] + interpo[0]*stressXYZ[i][7] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][5];
		for (unsigned j = 0; j < 8; j++)
		{
			stress[6 * j + i] = recovery[j];
		}
	}

}


