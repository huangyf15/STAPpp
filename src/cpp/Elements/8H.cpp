/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/8H.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CHex::CHex()
{
	NEN = 8;	// Each element has 2 nodes
	nodes = new CNode*[NEN];
    
    ND = 24; // number of degree of freedom of each element
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = nullptr;
}

//	Destructor
CHex::~CHex()
{
}

//	Read element data from stream Input
bool CHex::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial = dynamic_cast<CHexMaterial*>(MaterialSets) + MSet - 1; //stored in local node number
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	nodes[3] = &NodeList[N4 - 1];
	nodes[4] = &NodeList[N5 - 1];
	nodes[5] = &NodeList[N6 - 1];
	nodes[6] = &NodeList[N7 - 1];
	nodes[7] = &NodeList[N8 - 1];
	return true;
}

//	Write element data to stream
void CHex::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber
		   << setw(9) << nodes[1]->NodeNumber
		   << setw(9) << nodes[2]->NodeNumber
		   << setw(9) << nodes[3]->NodeNumber
		   << setw(9) << nodes[4]->NodeNumber
		   << setw(9) << nodes[5]->NodeNumber
		   << setw(9) << nodes[6]->NodeNumber
		   << setw(9) << nodes[7]->NodeNumber
		   << setw(12) << ElementMaterial->nset << endl;
}

//	Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CHex::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)  // dimension of FEM space
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CHex::SizeOfStiffnessMatrix() { return 300; } // reture the number df stiffness matrix

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CHex::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//Construct constitutive matrix
	CHexMaterial* material = dynamic_cast<CHexMaterial*>(ElementMaterial);	// Pointer to material of the element
	double v = material->nu;
	double k = material->E * (1-v)/(1+v)/(1-2*v);
	double D[3];

	D[0] = k;
	D[1] = k * v / (1 - v);
	D[2] = k * (1 - 2 * v) / 2 / (1 - v);

	//construct coordinate matrix
	double COORXYZ[24];
	for (unsigned int i = 0; i < 8; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			COORXYZ[3 * i + j]=nodes[i]->XYZ[j];
		}
	}

	// construct Jacobi matrix
	const double xi8[8]   = { 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 };
	const double eta8[8]  = {-0.577350269189626 , 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,0.577350269189626 , 0.577350269189626 ,-0.577350269189626 };
	const double zeta8[8] = {-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 , 0.577350269189626 , 0.577350269189626 };

	for (unsigned p = 0; p < 8; p++)
	{
		double xi   = xi8[p];
		double eta  = eta8[p];
		double zeta = zeta8[p];

		double GN[12];
		GN[0] = (1-eta)*(1-zeta);
		GN[1] = (1+eta)*(1-zeta);
		GN[2] = (1-eta)*(1+zeta);
		GN[3] = (1+eta)*(1+zeta);
		GN[4] = (1+xi)*(1-zeta);
		GN[5] = (1-xi)*(1-zeta);
		GN[6] = (1+xi)*(1+zeta);
		GN[7] = (1-xi)*(1+zeta);
		GN[8] = (1+xi)*(1-eta);
		GN[9] = (1+xi)*(1+eta);
		GN[10] = (1-xi)*(1+eta);
		GN[11] = (1-xi)*(1-eta);

		double J[9];
		J[0] = COORXYZ[0]*GN[0]+COORXYZ[3]*GN[1]-COORXYZ[6]*GN[1]-COORXYZ[9]*GN[0]+COORXYZ[12]*GN[2]+COORXYZ[15]*GN[3]-COORXYZ[18]*GN[3]-COORXYZ[21]*GN[2];
		J[1] = COORXYZ[0]*GN[4]-COORXYZ[3]*GN[4]+COORXYZ[6]*GN[5]-COORXYZ[9]*GN[5]-COORXYZ[12]*GN[6]+COORXYZ[15]*GN[6]+COORXYZ[18]*GN[7]-COORXYZ[21]*GN[7];
		J[2] = -COORXYZ[0]*GN[8]-COORXYZ[3]*GN[9]-COORXYZ[6]*GN[10]-COORXYZ[9]*GN[11]+COORXYZ[12]*GN[8]+COORXYZ[15]*GN[9]+COORXYZ[18]*GN[10]+COORXYZ[21]*GN[11];
		J[3] = COORXYZ[1]*GN[0]+COORXYZ[4]*GN[1]-COORXYZ[7]*GN[1]-COORXYZ[10]*GN[0]+COORXYZ[13]*GN[2]+COORXYZ[16]*GN[3]-COORXYZ[19]*GN[3]-COORXYZ[22]*GN[2];
		J[4] = COORXYZ[1]*GN[4]-COORXYZ[4]*GN[4]+COORXYZ[7]*GN[5]-COORXYZ[10]*GN[5]-COORXYZ[13]*GN[6]+COORXYZ[16]*GN[6]+COORXYZ[19]*GN[7]-COORXYZ[22]*GN[7];
		J[5] = -COORXYZ[1]*GN[8]-COORXYZ[4]*GN[9]-COORXYZ[7]*GN[10]-COORXYZ[10]*GN[11]+COORXYZ[13]*GN[8]+COORXYZ[16]*GN[9]+COORXYZ[19]*GN[10]+COORXYZ[22]*GN[11];
		J[6] = COORXYZ[2]*GN[0]+COORXYZ[5]*GN[1]-COORXYZ[8]*GN[1]-COORXYZ[11]*GN[0]+COORXYZ[14]*GN[2]+COORXYZ[17]*GN[3]-COORXYZ[20]*GN[3]-COORXYZ[23]*GN[2];
		J[7] = COORXYZ[2]*GN[4]-COORXYZ[5]*GN[4]+COORXYZ[8]*GN[5]-COORXYZ[11]*GN[5]-COORXYZ[14]*GN[6]+COORXYZ[17]*GN[6]+COORXYZ[20]*GN[7]-COORXYZ[23]*GN[7];
		J[8] = -COORXYZ[2]*GN[8]-COORXYZ[5]*GN[9]-COORXYZ[8]*GN[10]-COORXYZ[11]*GN[11]+COORXYZ[14]*GN[8]+COORXYZ[17]*GN[9]+COORXYZ[20]*GN[10]+COORXYZ[23]*GN[11];

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
		kerB[0] = GN[0]*invJ[0]+GN[4]*invJ[3]-GN[8]*invJ[6];
		kerB[1] = GN[0]*invJ[1]+GN[4]*invJ[4]-GN[8]*invJ[7];
		kerB[2] = GN[0]*invJ[2]+GN[4]*invJ[5]-GN[8]*invJ[8];
		kerB[3] = GN[1]*invJ[0]-GN[4]*invJ[3]-GN[9]*invJ[6];
		kerB[4] = GN[1]*invJ[1]-GN[4]*invJ[4]-GN[9]*invJ[7];
		kerB[5] = GN[1]*invJ[2]-GN[4]*invJ[5]-GN[9]*invJ[8];
		kerB[6] = -GN[1]*invJ[0]+GN[5]*invJ[3]-GN[10]*invJ[6];
		kerB[7] = -GN[1]*invJ[1]+GN[5]*invJ[4]-GN[10]*invJ[7];
		kerB[8] = -GN[1]*invJ[2]+GN[5]*invJ[5]-GN[10]*invJ[8];
		kerB[9] = -GN[0]*invJ[0]-GN[5]*invJ[3]-GN[11]*invJ[6];
		kerB[10] = -GN[0]*invJ[1]-GN[5]*invJ[4]-GN[11]*invJ[7];
		kerB[11] = -GN[0]*invJ[2]-GN[5]*invJ[5]-GN[11]*invJ[8];
		kerB[12] = GN[2]*invJ[0]-GN[6]*invJ[3]+GN[8]*invJ[6];
		kerB[13] = GN[2]*invJ[1]-GN[6]*invJ[4]+GN[8]*invJ[7];
		kerB[14] = GN[2]*invJ[2]-GN[6]*invJ[5]+GN[8]*invJ[8];
		kerB[15] = GN[3]*invJ[0]+GN[6]*invJ[3]+GN[9]*invJ[6];
		kerB[16] = GN[3]*invJ[1]+GN[6]*invJ[4]+GN[9]*invJ[7];
		kerB[17] = GN[3]*invJ[2]+GN[6]*invJ[5]+GN[9]*invJ[8];
		kerB[18] = -GN[3]*invJ[0]+GN[7]*invJ[3]+GN[10]*invJ[6];
		kerB[19] = -GN[3]*invJ[1]+GN[7]*invJ[4]+GN[10]*invJ[7];
		kerB[20] = -GN[3]*invJ[2]+GN[7]*invJ[5]+GN[10]*invJ[8];
		kerB[21] = -GN[2]*invJ[0]-GN[7]*invJ[3]+GN[11]*invJ[6];
		kerB[22] = -GN[2]*invJ[1]-GN[7]*invJ[4]+GN[11]*invJ[7];
		kerB[23] = -GN[2]*invJ[2]-GN[7]*invJ[5]+GN[11]*invJ[8];

		// construct element stiffness matrix
		Matrix[0] = detJ*(D[0]*(kerB[0]*kerB[0])+D[2]*(kerB[1]*kerB[1])+D[2]*(kerB[2]*kerB[2]));
		Matrix[1] = detJ*(D[0]*(kerB[1]*kerB[1])+D[2]*(kerB[0]*kerB[0])+D[2]*(kerB[2]*kerB[2]));
		Matrix[2] = detJ*(D[1]*kerB[0]*kerB[1]+D[2]*kerB[0]*kerB[1]);
		Matrix[3] = detJ*(D[0]*(kerB[2]*kerB[2])+D[2]*(kerB[0]*kerB[0])+D[2]*(kerB[1]*kerB[1]));
		Matrix[4] = detJ*(D[1]*kerB[1]*kerB[2]+D[2]*kerB[1]*kerB[2]);
		Matrix[5] = detJ*(D[1]*kerB[0]*kerB[2]+D[2]*kerB[0]*kerB[2]);
		Matrix[6] = detJ*(D[0]*(kerB[3]*kerB[3])+D[2]*(kerB[4]*kerB[4])+D[2]*(kerB[5]*kerB[5]));
		Matrix[7] = detJ*(D[1]*kerB[2]*kerB[3]+D[2]*kerB[0]*kerB[5]);
		Matrix[8] = detJ*(D[1]*kerB[1]*kerB[3]+D[2]*kerB[0]*kerB[4]);
		Matrix[9] = detJ*(D[0]*kerB[0]*kerB[3]+D[2]*kerB[1]*kerB[4]+D[2]*kerB[2]*kerB[5]);
		Matrix[10] = detJ*(D[0]*(kerB[4]*kerB[4])+D[2]*(kerB[3]*kerB[3])+D[2]*(kerB[5]*kerB[5]));
		Matrix[11] = detJ*(D[1]*kerB[3]*kerB[4]+D[2]*kerB[3]*kerB[4]);
		Matrix[12] = detJ*(D[1]*kerB[2]*kerB[4]+D[2]*kerB[1]*kerB[5]);
		Matrix[13] = detJ*(D[0]*kerB[1]*kerB[4]+D[2]*kerB[0]*kerB[3]+D[2]*kerB[2]*kerB[5]);
		Matrix[14] = detJ*(D[1]*kerB[0]*kerB[4]+D[2]*kerB[1]*kerB[3]);
		Matrix[15] = detJ*(D[0]*(kerB[5]*kerB[5])+D[2]*(kerB[3]*kerB[3])+D[2]*(kerB[4]*kerB[4]));
		Matrix[16] = detJ*(D[1]*kerB[4]*kerB[5]+D[2]*kerB[4]*kerB[5]);
		Matrix[17] = detJ*(D[1]*kerB[3]*kerB[5]+D[2]*kerB[3]*kerB[5]);
		Matrix[18] = detJ*(D[2]*kerB[0]*kerB[3]+D[0]*kerB[2]*kerB[5]+D[2]*kerB[1]*kerB[4]);
		Matrix[19] = detJ*(D[1]*kerB[1]*kerB[5]+D[2]*kerB[2]*kerB[4]);
		Matrix[20] = detJ*(D[1]*kerB[0]*kerB[5]+D[2]*kerB[2]*kerB[3]);
		Matrix[21] = detJ*(D[0]*(kerB[6]*kerB[6])+D[2]*(kerB[7]*kerB[7])+D[2]*(kerB[8]*kerB[8]));
		Matrix[22] = detJ*(D[1]*kerB[5]*kerB[6]+D[2]*kerB[3]*kerB[8]);
		Matrix[23] = detJ*(D[1]*kerB[4]*kerB[6]+D[2]*kerB[3]*kerB[7]);
		Matrix[24] = detJ*(D[0]*kerB[3]*kerB[6]+D[2]*kerB[4]*kerB[7]+D[2]*kerB[5]*kerB[8]);
		Matrix[25] = detJ*(D[1]*kerB[2]*kerB[6]+D[2]*kerB[0]*kerB[8]);
		Matrix[26] = detJ*(D[1]*kerB[1]*kerB[6]+D[2]*kerB[0]*kerB[7]);
		Matrix[27] = detJ*(D[0]*kerB[0]*kerB[6]+D[2]*kerB[1]*kerB[7]+D[2]*kerB[2]*kerB[8]);
		Matrix[28] = detJ*(D[0]*(kerB[7]*kerB[7])+D[2]*(kerB[6]*kerB[6])+D[2]*(kerB[8]*kerB[8]));
		Matrix[29] = detJ*(D[1]*kerB[6]*kerB[7]+D[2]*kerB[6]*kerB[7]);
		Matrix[30] = detJ*(D[1]*kerB[5]*kerB[7]+D[2]*kerB[4]*kerB[8]);
		Matrix[31] = detJ*(D[0]*kerB[4]*kerB[7]+D[2]*kerB[3]*kerB[6]+D[2]*kerB[5]*kerB[8]);
		Matrix[32] = detJ*(D[1]*kerB[3]*kerB[7]+D[2]*kerB[4]*kerB[6]);
		Matrix[33] = detJ*(D[1]*kerB[2]*kerB[7]+D[2]*kerB[1]*kerB[8]);
		Matrix[34] = detJ*(D[0]*kerB[1]*kerB[7]+D[2]*kerB[0]*kerB[6]+D[2]*kerB[2]*kerB[8]);
		Matrix[35] = detJ*(D[1]*kerB[0]*kerB[7]+D[2]*kerB[1]*kerB[6]);
		Matrix[36] = detJ*(D[0]*(kerB[8]*kerB[8])+D[2]*(kerB[6]*kerB[6])+D[2]*(kerB[7]*kerB[7]));
		Matrix[37] = detJ*(D[1]*kerB[7]*kerB[8]+D[2]*kerB[7]*kerB[8]);
		Matrix[38] = detJ*(D[1]*kerB[6]*kerB[8]+D[2]*kerB[6]*kerB[8]);
		Matrix[39] = detJ*(D[2]*kerB[3]*kerB[6]+D[0]*kerB[5]*kerB[8]+D[2]*kerB[4]*kerB[7]);
		Matrix[40] = detJ*(D[1]*kerB[4]*kerB[8]+D[2]*kerB[5]*kerB[7]);
		Matrix[41] = detJ*(D[1]*kerB[3]*kerB[8]+D[2]*kerB[5]*kerB[6]);
		Matrix[42] = detJ*(D[2]*kerB[0]*kerB[6]+D[0]*kerB[2]*kerB[8]+D[2]*kerB[1]*kerB[7]);
		Matrix[43] = detJ*(D[1]*kerB[1]*kerB[8]+D[2]*kerB[2]*kerB[7]);
		Matrix[44] = detJ*(D[1]*kerB[0]*kerB[8]+D[2]*kerB[2]*kerB[6]);
		Matrix[45] = detJ*(D[0]*(kerB[9]*kerB[9])+D[2]*(kerB[10]*kerB[10])+D[2]*(kerB[11]*kerB[11]));
		Matrix[46] = detJ*(D[1]*kerB[8]*kerB[9]+D[2]*kerB[6]*kerB[11]);
		Matrix[47] = detJ*(D[1]*kerB[7]*kerB[9]+D[2]*kerB[6]*kerB[10]);
		Matrix[48] = detJ*(D[0]*kerB[6]*kerB[9]+D[2]*kerB[7]*kerB[10]+D[2]*kerB[8]*kerB[11]);
		Matrix[49] = detJ*(D[1]*kerB[5]*kerB[9]+D[2]*kerB[3]*kerB[11]);
		Matrix[50] = detJ*(D[1]*kerB[4]*kerB[9]+D[2]*kerB[3]*kerB[10]);
		Matrix[51] = detJ*(D[0]*kerB[3]*kerB[9]+D[2]*kerB[4]*kerB[10]+D[2]*kerB[5]*kerB[11]);
		Matrix[52] = detJ*(D[1]*kerB[2]*kerB[9]+D[2]*kerB[0]*kerB[11]);
		Matrix[53] = detJ*(D[1]*kerB[1]*kerB[9]+D[2]*kerB[0]*kerB[10]);
		Matrix[54] = detJ*(D[0]*kerB[0]*kerB[9]+D[2]*kerB[1]*kerB[10]+D[2]*kerB[2]*kerB[11]);
		Matrix[55] = detJ*(D[0]*(kerB[10]*kerB[10])+D[2]*(kerB[9]*kerB[9])+D[2]*(kerB[11]*kerB[11]));
		Matrix[56] = detJ*(D[1]*kerB[9]*kerB[10]+D[2]*kerB[9]*kerB[10]);
		Matrix[57] = detJ*(D[1]*kerB[8]*kerB[10]+D[2]*kerB[7]*kerB[11]);
		Matrix[58] = detJ*(D[0]*kerB[7]*kerB[10]+D[2]*kerB[6]*kerB[9]+D[2]*kerB[8]*kerB[11]);
		Matrix[59] = detJ*(D[1]*kerB[6]*kerB[10]+D[2]*kerB[7]*kerB[9]);
		Matrix[60] = detJ*(D[1]*kerB[5]*kerB[10]+D[2]*kerB[4]*kerB[11]);
		Matrix[61] = detJ*(D[0]*kerB[4]*kerB[10]+D[2]*kerB[3]*kerB[9]+D[2]*kerB[5]*kerB[11]);
		Matrix[62] = detJ*(D[1]*kerB[3]*kerB[10]+D[2]*kerB[4]*kerB[9]);
		Matrix[63] = detJ*(D[1]*kerB[2]*kerB[10]+D[2]*kerB[1]*kerB[11]);
		Matrix[64] = detJ*(D[0]*kerB[1]*kerB[10]+D[2]*kerB[0]*kerB[9]+D[2]*kerB[2]*kerB[11]);
		Matrix[65] = detJ*(D[1]*kerB[0]*kerB[10]+D[2]*kerB[1]*kerB[9]);
		Matrix[66] = detJ*(D[0]*(kerB[11]*kerB[11])+D[2]*(kerB[9]*kerB[9])+D[2]*(kerB[10]*kerB[10]));
		Matrix[67] = detJ*(D[1]*kerB[10]*kerB[11]+D[2]*kerB[10]*kerB[11]);
		Matrix[68] = detJ*(D[1]*kerB[9]*kerB[11]+D[2]*kerB[9]*kerB[11]);
		Matrix[69] = detJ*(D[2]*kerB[6]*kerB[9]+D[0]*kerB[8]*kerB[11]+D[2]*kerB[7]*kerB[10]);
		Matrix[70] = detJ*(D[1]*kerB[7]*kerB[11]+D[2]*kerB[8]*kerB[10]);
		Matrix[71] = detJ*(D[1]*kerB[6]*kerB[11]+D[2]*kerB[8]*kerB[9]);
		Matrix[72] = detJ*(D[2]*kerB[3]*kerB[9]+D[0]*kerB[5]*kerB[11]+D[2]*kerB[4]*kerB[10]);
		Matrix[73] = detJ*(D[1]*kerB[4]*kerB[11]+D[2]*kerB[5]*kerB[10]);
		Matrix[74] = detJ*(D[1]*kerB[3]*kerB[11]+D[2]*kerB[5]*kerB[9]);
		Matrix[75] = detJ*(D[2]*kerB[0]*kerB[9]+D[0]*kerB[2]*kerB[11]+D[2]*kerB[1]*kerB[10]);
		Matrix[76] = detJ*(D[1]*kerB[1]*kerB[11]+D[2]*kerB[2]*kerB[10]);
		Matrix[77] = detJ*(D[1]*kerB[0]*kerB[11]+D[2]*kerB[2]*kerB[9]);
		Matrix[78] = detJ*(D[0]*(kerB[12]*kerB[12])+D[2]*(kerB[13]*kerB[13])+D[2]*(kerB[14]*kerB[14]));
		Matrix[79] = detJ*(D[1]*kerB[11]*kerB[12]+D[2]*kerB[9]*kerB[14]);
		Matrix[80] = detJ*(D[1]*kerB[10]*kerB[12]+D[2]*kerB[9]*kerB[13]);
		Matrix[81] = detJ*(D[0]*kerB[9]*kerB[12]+D[2]*kerB[10]*kerB[13]+D[2]*kerB[11]*kerB[14]);
		Matrix[82] = detJ*(D[1]*kerB[8]*kerB[12]+D[2]*kerB[6]*kerB[14]);
		Matrix[83] = detJ*(D[1]*kerB[7]*kerB[12]+D[2]*kerB[6]*kerB[13]);
		Matrix[84] = detJ*(D[0]*kerB[6]*kerB[12]+D[2]*kerB[7]*kerB[13]+D[2]*kerB[8]*kerB[14]);
		Matrix[85] = detJ*(D[1]*kerB[5]*kerB[12]+D[2]*kerB[3]*kerB[14]);
		Matrix[86] = detJ*(D[1]*kerB[4]*kerB[12]+D[2]*kerB[3]*kerB[13]);
		Matrix[87] = detJ*(D[0]*kerB[3]*kerB[12]+D[2]*kerB[4]*kerB[13]+D[2]*kerB[5]*kerB[14]);
		Matrix[88] = detJ*(D[1]*kerB[2]*kerB[12]+D[2]*kerB[0]*kerB[14]);
		Matrix[89] = detJ*(D[1]*kerB[1]*kerB[12]+D[2]*kerB[0]*kerB[13]);
		Matrix[90] = detJ*(D[0]*kerB[0]*kerB[12]+D[2]*kerB[1]*kerB[13]+D[2]*kerB[2]*kerB[14]);
		Matrix[91] = detJ*(D[0]*(kerB[13]*kerB[13])+D[2]*(kerB[12]*kerB[12])+D[2]*(kerB[14]*kerB[14]));
		Matrix[92] = detJ*(D[1]*kerB[12]*kerB[13]+D[2]*kerB[12]*kerB[13]);
		Matrix[93] = detJ*(D[1]*kerB[11]*kerB[13]+D[2]*kerB[10]*kerB[14]);
		Matrix[94] = detJ*(D[0]*kerB[10]*kerB[13]+D[2]*kerB[9]*kerB[12]+D[2]*kerB[11]*kerB[14]);
		Matrix[95] = detJ*(D[1]*kerB[9]*kerB[13]+D[2]*kerB[10]*kerB[12]);
		Matrix[96] = detJ*(D[1]*kerB[8]*kerB[13]+D[2]*kerB[7]*kerB[14]);
		Matrix[97] = detJ*(D[0]*kerB[7]*kerB[13]+D[2]*kerB[6]*kerB[12]+D[2]*kerB[8]*kerB[14]);
		Matrix[98] = detJ*(D[1]*kerB[6]*kerB[13]+D[2]*kerB[7]*kerB[12]);
		Matrix[99] = detJ*(D[1]*kerB[5]*kerB[13]+D[2]*kerB[4]*kerB[14]);
		Matrix[100] = detJ*(D[0]*kerB[4]*kerB[13]+D[2]*kerB[3]*kerB[12]+D[2]*kerB[5]*kerB[14]);
		Matrix[101] = detJ*(D[1]*kerB[3]*kerB[13]+D[2]*kerB[4]*kerB[12]);
		Matrix[102] = detJ*(D[1]*kerB[2]*kerB[13]+D[2]*kerB[1]*kerB[14]);
		Matrix[103] = detJ*(D[0]*kerB[1]*kerB[13]+D[2]*kerB[0]*kerB[12]+D[2]*kerB[2]*kerB[14]);
		Matrix[104] = detJ*(D[1]*kerB[0]*kerB[13]+D[2]*kerB[1]*kerB[12]);
		Matrix[105] = detJ*(D[0]*(kerB[14]*kerB[14])+D[2]*(kerB[12]*kerB[12])+D[2]*(kerB[13]*kerB[13]));
		Matrix[106] = detJ*(D[1]*kerB[13]*kerB[14]+D[2]*kerB[13]*kerB[14]);
		Matrix[107] = detJ*(D[1]*kerB[12]*kerB[14]+D[2]*kerB[12]*kerB[14]);
		Matrix[108] = detJ*(D[2]*kerB[9]*kerB[12]+D[0]*kerB[11]*kerB[14]+D[2]*kerB[10]*kerB[13]);
		Matrix[109] = detJ*(D[1]*kerB[10]*kerB[14]+D[2]*kerB[11]*kerB[13]);
		Matrix[110] = detJ*(D[1]*kerB[9]*kerB[14]+D[2]*kerB[11]*kerB[12]);
		Matrix[111] = detJ*(D[2]*kerB[6]*kerB[12]+D[0]*kerB[8]*kerB[14]+D[2]*kerB[7]*kerB[13]);
		Matrix[112] = detJ*(D[1]*kerB[7]*kerB[14]+D[2]*kerB[8]*kerB[13]);
		Matrix[113] = detJ*(D[1]*kerB[6]*kerB[14]+D[2]*kerB[8]*kerB[12]);
		Matrix[114] = detJ*(D[2]*kerB[3]*kerB[12]+D[0]*kerB[5]*kerB[14]+D[2]*kerB[4]*kerB[13]);
		Matrix[115] = detJ*(D[1]*kerB[4]*kerB[14]+D[2]*kerB[5]*kerB[13]);
		Matrix[116] = detJ*(D[1]*kerB[3]*kerB[14]+D[2]*kerB[5]*kerB[12]);
		Matrix[117] = detJ*(D[2]*kerB[0]*kerB[12]+D[0]*kerB[2]*kerB[14]+D[2]*kerB[1]*kerB[13]);
		Matrix[118] = detJ*(D[1]*kerB[1]*kerB[14]+D[2]*kerB[2]*kerB[13]);
		Matrix[119] = detJ*(D[1]*kerB[0]*kerB[14]+D[2]*kerB[2]*kerB[12]);
		Matrix[120] = detJ*(D[0]*(kerB[15]*kerB[15])+D[2]*(kerB[16]*kerB[16])+D[2]*(kerB[17]*kerB[17]));
		Matrix[121] = detJ*(D[1]*kerB[14]*kerB[15]+D[2]*kerB[12]*kerB[17]);
		Matrix[122] = detJ*(D[1]*kerB[13]*kerB[15]+D[2]*kerB[12]*kerB[16]);
		Matrix[123] = detJ*(D[0]*kerB[12]*kerB[15]+D[2]*kerB[13]*kerB[16]+D[2]*kerB[14]*kerB[17]);
		Matrix[124] = detJ*(D[1]*kerB[11]*kerB[15]+D[2]*kerB[9]*kerB[17]);
		Matrix[125] = detJ*(D[1]*kerB[10]*kerB[15]+D[2]*kerB[9]*kerB[16]);
		Matrix[126] = detJ*(D[0]*kerB[9]*kerB[15]+D[2]*kerB[10]*kerB[16]+D[2]*kerB[11]*kerB[17]);
		Matrix[127] = detJ*(D[1]*kerB[8]*kerB[15]+D[2]*kerB[6]*kerB[17]);
		Matrix[128] = detJ*(D[1]*kerB[7]*kerB[15]+D[2]*kerB[6]*kerB[16]);
		Matrix[129] = detJ*(D[0]*kerB[6]*kerB[15]+D[2]*kerB[7]*kerB[16]+D[2]*kerB[8]*kerB[17]);
		Matrix[130] = detJ*(D[1]*kerB[5]*kerB[15]+D[2]*kerB[3]*kerB[17]);
		Matrix[131] = detJ*(D[1]*kerB[4]*kerB[15]+D[2]*kerB[3]*kerB[16]);
		Matrix[132] = detJ*(D[0]*kerB[3]*kerB[15]+D[2]*kerB[4]*kerB[16]+D[2]*kerB[5]*kerB[17]);
		Matrix[133] = detJ*(D[1]*kerB[2]*kerB[15]+D[2]*kerB[0]*kerB[17]);
		Matrix[134] = detJ*(D[1]*kerB[1]*kerB[15]+D[2]*kerB[0]*kerB[16]);
		Matrix[135] = detJ*(D[0]*kerB[0]*kerB[15]+D[2]*kerB[1]*kerB[16]+D[2]*kerB[2]*kerB[17]);
		Matrix[136] = detJ*(D[0]*(kerB[16]*kerB[16])+D[2]*(kerB[15]*kerB[15])+D[2]*(kerB[17]*kerB[17]));
		Matrix[137] = detJ*(D[1]*kerB[15]*kerB[16]+D[2]*kerB[15]*kerB[16]);
		Matrix[138] = detJ*(D[1]*kerB[14]*kerB[16]+D[2]*kerB[13]*kerB[17]);
		Matrix[139] = detJ*(D[0]*kerB[13]*kerB[16]+D[2]*kerB[12]*kerB[15]+D[2]*kerB[14]*kerB[17]);
		Matrix[140] = detJ*(D[1]*kerB[12]*kerB[16]+D[2]*kerB[13]*kerB[15]);
		Matrix[141] = detJ*(D[1]*kerB[11]*kerB[16]+D[2]*kerB[10]*kerB[17]);
		Matrix[142] = detJ*(D[0]*kerB[10]*kerB[16]+D[2]*kerB[9]*kerB[15]+D[2]*kerB[11]*kerB[17]);
		Matrix[143] = detJ*(D[1]*kerB[9]*kerB[16]+D[2]*kerB[10]*kerB[15]);
		Matrix[144] = detJ*(D[1]*kerB[8]*kerB[16]+D[2]*kerB[7]*kerB[17]);
		Matrix[145] = detJ*(D[0]*kerB[7]*kerB[16]+D[2]*kerB[6]*kerB[15]+D[2]*kerB[8]*kerB[17]);
		Matrix[146] = detJ*(D[1]*kerB[6]*kerB[16]+D[2]*kerB[7]*kerB[15]);
		Matrix[147] = detJ*(D[1]*kerB[5]*kerB[16]+D[2]*kerB[4]*kerB[17]);
		Matrix[148] = detJ*(D[0]*kerB[4]*kerB[16]+D[2]*kerB[3]*kerB[15]+D[2]*kerB[5]*kerB[17]);
		Matrix[149] = detJ*(D[1]*kerB[3]*kerB[16]+D[2]*kerB[4]*kerB[15]);
		Matrix[150] = detJ*(D[1]*kerB[2]*kerB[16]+D[2]*kerB[1]*kerB[17]);
		Matrix[151] = detJ*(D[0]*kerB[1]*kerB[16]+D[2]*kerB[0]*kerB[15]+D[2]*kerB[2]*kerB[17]);
		Matrix[152] = detJ*(D[1]*kerB[0]*kerB[16]+D[2]*kerB[1]*kerB[15]);
		Matrix[153] = detJ*(D[0]*(kerB[17]*kerB[17])+D[2]*(kerB[15]*kerB[15])+D[2]*(kerB[16]*kerB[16]));
		Matrix[154] = detJ*(D[1]*kerB[16]*kerB[17]+D[2]*kerB[16]*kerB[17]);
		Matrix[155] = detJ*(D[1]*kerB[15]*kerB[17]+D[2]*kerB[15]*kerB[17]);
		Matrix[156] = detJ*(D[2]*kerB[12]*kerB[15]+D[0]*kerB[14]*kerB[17]+D[2]*kerB[13]*kerB[16]);
		Matrix[157] = detJ*(D[1]*kerB[13]*kerB[17]+D[2]*kerB[14]*kerB[16]);
		Matrix[158] = detJ*(D[1]*kerB[12]*kerB[17]+D[2]*kerB[14]*kerB[15]);
		Matrix[159] = detJ*(D[2]*kerB[9]*kerB[15]+D[0]*kerB[11]*kerB[17]+D[2]*kerB[10]*kerB[16]);
		Matrix[160] = detJ*(D[1]*kerB[10]*kerB[17]+D[2]*kerB[11]*kerB[16]);
		Matrix[161] = detJ*(D[1]*kerB[9]*kerB[17]+D[2]*kerB[11]*kerB[15]);
		Matrix[162] = detJ*(D[2]*kerB[6]*kerB[15]+D[0]*kerB[8]*kerB[17]+D[2]*kerB[7]*kerB[16]);
		Matrix[163] = detJ*(D[1]*kerB[7]*kerB[17]+D[2]*kerB[8]*kerB[16]);
		Matrix[164] = detJ*(D[1]*kerB[6]*kerB[17]+D[2]*kerB[8]*kerB[15]);
		Matrix[165] = detJ*(D[2]*kerB[3]*kerB[15]+D[0]*kerB[5]*kerB[17]+D[2]*kerB[4]*kerB[16]);
		Matrix[166] = detJ*(D[1]*kerB[4]*kerB[17]+D[2]*kerB[5]*kerB[16]);
		Matrix[167] = detJ*(D[1]*kerB[3]*kerB[17]+D[2]*kerB[5]*kerB[15]);
		Matrix[168] = detJ*(D[2]*kerB[0]*kerB[15]+D[0]*kerB[2]*kerB[17]+D[2]*kerB[1]*kerB[16]);
		Matrix[169] = detJ*(D[1]*kerB[1]*kerB[17]+D[2]*kerB[2]*kerB[16]);
		Matrix[170] = detJ*(D[1]*kerB[0]*kerB[17]+D[2]*kerB[2]*kerB[15]);
		Matrix[171] = detJ*(D[0]*(kerB[18]*kerB[18])+D[2]*(kerB[19]*kerB[19])+D[2]*(kerB[20]*kerB[20]));
		Matrix[172] = detJ*(D[1]*kerB[17]*kerB[18]+D[2]*kerB[15]*kerB[20]);
		Matrix[173] = detJ*(D[1]*kerB[16]*kerB[18]+D[2]*kerB[15]*kerB[19]);
		Matrix[174] = detJ*(D[0]*kerB[15]*kerB[18]+D[2]*kerB[16]*kerB[19]+D[2]*kerB[17]*kerB[20]);
		Matrix[175] = detJ*(D[1]*kerB[14]*kerB[18]+D[2]*kerB[12]*kerB[20]);
		Matrix[176] = detJ*(D[1]*kerB[13]*kerB[18]+D[2]*kerB[12]*kerB[19]);
		Matrix[177] = detJ*(D[0]*kerB[12]*kerB[18]+D[2]*kerB[13]*kerB[19]+D[2]*kerB[14]*kerB[20]);
		Matrix[178] = detJ*(D[1]*kerB[11]*kerB[18]+D[2]*kerB[9]*kerB[20]);
		Matrix[179] = detJ*(D[1]*kerB[10]*kerB[18]+D[2]*kerB[9]*kerB[19]);
		Matrix[180] = detJ*(D[0]*kerB[9]*kerB[18]+D[2]*kerB[10]*kerB[19]+D[2]*kerB[11]*kerB[20]);
		Matrix[181] = detJ*(D[1]*kerB[8]*kerB[18]+D[2]*kerB[6]*kerB[20]);
		Matrix[182] = detJ*(D[1]*kerB[7]*kerB[18]+D[2]*kerB[6]*kerB[19]);
		Matrix[183] = detJ*(D[0]*kerB[6]*kerB[18]+D[2]*kerB[7]*kerB[19]+D[2]*kerB[8]*kerB[20]);
		Matrix[184] = detJ*(D[1]*kerB[5]*kerB[18]+D[2]*kerB[3]*kerB[20]);
		Matrix[185] = detJ*(D[1]*kerB[4]*kerB[18]+D[2]*kerB[3]*kerB[19]);
		Matrix[186] = detJ*(D[0]*kerB[3]*kerB[18]+D[2]*kerB[4]*kerB[19]+D[2]*kerB[5]*kerB[20]);
		Matrix[187] = detJ*(D[1]*kerB[2]*kerB[18]+D[2]*kerB[0]*kerB[20]);
		Matrix[188] = detJ*(D[1]*kerB[1]*kerB[18]+D[2]*kerB[0]*kerB[19]);
		Matrix[189] = detJ*(D[0]*kerB[0]*kerB[18]+D[2]*kerB[1]*kerB[19]+D[2]*kerB[2]*kerB[20]);
		Matrix[190] = detJ*(D[0]*(kerB[19]*kerB[19])+D[2]*(kerB[20]*kerB[20])*2.0);
		Matrix[191] = detJ*(D[1]*kerB[18]*kerB[19]+D[2]*kerB[19]*kerB[20]);
		Matrix[192] = detJ*(D[1]*kerB[17]*kerB[19]+D[2]*kerB[16]*kerB[20]);
		Matrix[193] = detJ*(D[0]*kerB[16]*kerB[19]+D[2]*kerB[15]*kerB[20]+D[2]*kerB[17]*kerB[20]);
		Matrix[194] = detJ*(D[1]*kerB[15]*kerB[19]+D[2]*kerB[16]*kerB[20]);
		Matrix[195] = detJ*(D[1]*kerB[14]*kerB[19]+D[2]*kerB[13]*kerB[20]);
		Matrix[196] = detJ*(D[0]*kerB[13]*kerB[19]+D[2]*kerB[12]*kerB[20]+D[2]*kerB[14]*kerB[20]);
		Matrix[197] = detJ*(D[1]*kerB[12]*kerB[19]+D[2]*kerB[13]*kerB[20]);
		Matrix[198] = detJ*(D[1]*kerB[11]*kerB[19]+D[2]*kerB[10]*kerB[20]);
		Matrix[199] = detJ*(D[0]*kerB[10]*kerB[19]+D[2]*kerB[9]*kerB[20]+D[2]*kerB[11]*kerB[20]);
		Matrix[200] = detJ*(D[1]*kerB[9]*kerB[19]+D[2]*kerB[10]*kerB[20]);
		Matrix[201] = detJ*(D[1]*kerB[8]*kerB[19]+D[2]*kerB[7]*kerB[20]);
		Matrix[202] = detJ*(D[0]*kerB[7]*kerB[19]+D[2]*kerB[6]*kerB[20]+D[2]*kerB[8]*kerB[20]);
		Matrix[203] = detJ*(D[1]*kerB[6]*kerB[19]+D[2]*kerB[7]*kerB[20]);
		Matrix[204] = detJ*(D[1]*kerB[5]*kerB[19]+D[2]*kerB[4]*kerB[20]);
		Matrix[205] = detJ*(D[0]*kerB[4]*kerB[19]+D[2]*kerB[3]*kerB[20]+D[2]*kerB[5]*kerB[20]);
		Matrix[206] = detJ*(D[1]*kerB[3]*kerB[19]+D[2]*kerB[4]*kerB[20]);
		Matrix[207] = detJ*(D[1]*kerB[2]*kerB[19]+D[2]*kerB[1]*kerB[20]);
		Matrix[208] = detJ*(D[0]*kerB[1]*kerB[19]+D[2]*kerB[0]*kerB[20]+D[2]*kerB[2]*kerB[20]);
		Matrix[209] = detJ*(D[1]*kerB[0]*kerB[19]+D[2]*kerB[1]*kerB[20]);
		Matrix[210] = detJ*(D[0]*(kerB[20]*kerB[20])+D[2]*(kerB[18]*kerB[18])+D[2]*(kerB[19]*kerB[19]));
		Matrix[211] = detJ*(D[1]*kerB[19]*kerB[20]+D[2]*kerB[19]*kerB[20]);
		Matrix[212] = detJ*(D[1]*kerB[18]*kerB[20]+D[2]*kerB[18]*kerB[20]);
		Matrix[213] = detJ*(D[2]*kerB[15]*kerB[18]+D[0]*kerB[17]*kerB[20]+D[2]*kerB[16]*kerB[19]);
		Matrix[214] = detJ*(D[1]*kerB[16]*kerB[20]+D[2]*kerB[17]*kerB[19]);
		Matrix[215] = detJ*(D[1]*kerB[15]*kerB[20]+D[2]*kerB[17]*kerB[18]);
		Matrix[216] = detJ*(D[2]*kerB[12]*kerB[18]+D[0]*kerB[14]*kerB[20]+D[2]*kerB[13]*kerB[19]);
		Matrix[217] = detJ*(D[1]*kerB[13]*kerB[20]+D[2]*kerB[14]*kerB[19]);
		Matrix[218] = detJ*(D[1]*kerB[12]*kerB[20]+D[2]*kerB[14]*kerB[18]);
		Matrix[219] = detJ*(D[2]*kerB[9]*kerB[18]+D[0]*kerB[11]*kerB[20]+D[2]*kerB[10]*kerB[19]);
		Matrix[220] = detJ*(D[1]*kerB[10]*kerB[20]+D[2]*kerB[11]*kerB[19]);
		Matrix[221] = detJ*(D[1]*kerB[9]*kerB[20]+D[2]*kerB[11]*kerB[18]);
		Matrix[222] = detJ*(D[2]*kerB[6]*kerB[18]+D[0]*kerB[8]*kerB[20]+D[2]*kerB[7]*kerB[19]);
		Matrix[223] = detJ*(D[1]*kerB[7]*kerB[20]+D[2]*kerB[8]*kerB[19]);
		Matrix[224] = detJ*(D[1]*kerB[6]*kerB[20]+D[2]*kerB[8]*kerB[18]);
		Matrix[225] = detJ*(D[2]*kerB[3]*kerB[18]+D[0]*kerB[5]*kerB[20]+D[2]*kerB[4]*kerB[19]);
		Matrix[226] = detJ*(D[1]*kerB[4]*kerB[20]+D[2]*kerB[5]*kerB[19]);
		Matrix[227] = detJ*(D[1]*kerB[3]*kerB[20]+D[2]*kerB[5]*kerB[18]);
		Matrix[228] = detJ*(D[2]*kerB[0]*kerB[18]+D[0]*kerB[2]*kerB[20]+D[2]*kerB[1]*kerB[19]);
		Matrix[229] = detJ*(D[1]*kerB[1]*kerB[20]+D[2]*kerB[2]*kerB[19]);
		Matrix[230] = detJ*(D[1]*kerB[0]*kerB[20]+D[2]*kerB[2]*kerB[18]);
		Matrix[231] = detJ*(D[0]*(kerB[21]*kerB[21])+D[2]*(kerB[22]*kerB[22])+D[2]*(kerB[23]*kerB[23]));
		Matrix[232] = detJ*(D[1]*kerB[20]*kerB[21]+D[2]*kerB[18]*kerB[23]);
		Matrix[233] = detJ*(D[1]*kerB[19]*kerB[21]+D[2]*kerB[20]*kerB[22]);
		Matrix[234] = detJ*(D[0]*kerB[18]*kerB[21]+D[2]*kerB[19]*kerB[22]+D[2]*kerB[20]*kerB[23]);
		Matrix[235] = detJ*(D[1]*kerB[17]*kerB[21]+D[2]*kerB[15]*kerB[23]);
		Matrix[236] = detJ*(D[1]*kerB[16]*kerB[21]+D[2]*kerB[15]*kerB[22]);
		Matrix[237] = detJ*(D[0]*kerB[15]*kerB[21]+D[2]*kerB[16]*kerB[22]+D[2]*kerB[17]*kerB[23]);
		Matrix[238] = detJ*(D[1]*kerB[14]*kerB[21]+D[2]*kerB[12]*kerB[23]);
		Matrix[239] = detJ*(D[1]*kerB[13]*kerB[21]+D[2]*kerB[12]*kerB[22]);
		Matrix[240] = detJ*(D[0]*kerB[12]*kerB[21]+D[2]*kerB[13]*kerB[22]+D[2]*kerB[14]*kerB[23]);
		Matrix[241] = detJ*(D[1]*kerB[11]*kerB[21]+D[2]*kerB[9]*kerB[23]);
		Matrix[242] = detJ*(D[1]*kerB[10]*kerB[21]+D[2]*kerB[9]*kerB[22]);
		Matrix[243] = detJ*(D[0]*kerB[9]*kerB[21]+D[2]*kerB[10]*kerB[22]+D[2]*kerB[11]*kerB[23]);
		Matrix[244] = detJ*(D[1]*kerB[8]*kerB[21]+D[2]*kerB[6]*kerB[23]);
		Matrix[245] = detJ*(D[1]*kerB[7]*kerB[21]+D[2]*kerB[6]*kerB[22]);
		Matrix[246] = detJ*(D[0]*kerB[6]*kerB[21]+D[2]*kerB[7]*kerB[22]+D[2]*kerB[8]*kerB[23]);
		Matrix[247] = detJ*(D[1]*kerB[5]*kerB[21]+D[2]*kerB[3]*kerB[23]);
		Matrix[248] = detJ*(D[1]*kerB[4]*kerB[21]+D[2]*kerB[3]*kerB[22]);
		Matrix[249] = detJ*(D[0]*kerB[3]*kerB[21]+D[2]*kerB[4]*kerB[22]+D[2]*kerB[5]*kerB[23]);
		Matrix[250] = detJ*(D[1]*kerB[2]*kerB[21]+D[2]*kerB[0]*kerB[23]);
		Matrix[251] = detJ*(D[1]*kerB[1]*kerB[21]+D[2]*kerB[0]*kerB[22]);
		Matrix[252] = detJ*(D[0]*kerB[0]*kerB[21]+D[2]*kerB[1]*kerB[22]+D[2]*kerB[2]*kerB[23]);
		Matrix[253] = detJ*(D[0]*(kerB[22]*kerB[22])+D[2]*(kerB[21]*kerB[21])+D[2]*(kerB[23]*kerB[23]));
		Matrix[254] = detJ*(D[1]*kerB[21]*kerB[22]+D[2]*kerB[21]*kerB[22]);
		Matrix[255] = detJ*(D[1]*kerB[20]*kerB[22]+D[2]*kerB[19]*kerB[23]);
		Matrix[256] = detJ*(D[0]*kerB[19]*kerB[22]+D[2]*kerB[20]*kerB[21]+D[2]*kerB[20]*kerB[23]);
		Matrix[257] = detJ*(D[1]*kerB[18]*kerB[22]+D[2]*kerB[19]*kerB[21]);
		Matrix[258] = detJ*(D[1]*kerB[17]*kerB[22]+D[2]*kerB[16]*kerB[23]);
		Matrix[259] = detJ*(D[0]*kerB[16]*kerB[22]+D[2]*kerB[15]*kerB[21]+D[2]*kerB[17]*kerB[23]);
		Matrix[260] = detJ*(D[1]*kerB[15]*kerB[22]+D[2]*kerB[16]*kerB[21]);
		Matrix[261] = detJ*(D[1]*kerB[14]*kerB[22]+D[2]*kerB[13]*kerB[23]);
		Matrix[262] = detJ*(D[0]*kerB[13]*kerB[22]+D[2]*kerB[12]*kerB[21]+D[2]*kerB[14]*kerB[23]);
		Matrix[263] = detJ*(D[1]*kerB[12]*kerB[22]+D[2]*kerB[13]*kerB[21]);
		Matrix[264] = detJ*(D[1]*kerB[11]*kerB[22]+D[2]*kerB[10]*kerB[23]);
		Matrix[265] = detJ*(D[0]*kerB[10]*kerB[22]+D[2]*kerB[9]*kerB[21]+D[2]*kerB[11]*kerB[23]);
		Matrix[266] = detJ*(D[1]*kerB[9]*kerB[22]+D[2]*kerB[10]*kerB[21]);
		Matrix[267] = detJ*(D[1]*kerB[8]*kerB[22]+D[2]*kerB[7]*kerB[23]);
		Matrix[268] = detJ*(D[0]*kerB[7]*kerB[22]+D[2]*kerB[6]*kerB[21]+D[2]*kerB[8]*kerB[23]);
		Matrix[269] = detJ*(D[1]*kerB[6]*kerB[22]+D[2]*kerB[7]*kerB[21]);
		Matrix[270] = detJ*(D[1]*kerB[5]*kerB[22]+D[2]*kerB[4]*kerB[23]);
		Matrix[271] = detJ*(D[0]*kerB[4]*kerB[22]+D[2]*kerB[3]*kerB[21]+D[2]*kerB[5]*kerB[23]);
		Matrix[272] = detJ*(D[1]*kerB[3]*kerB[22]+D[2]*kerB[4]*kerB[21]);
		Matrix[273] = detJ*(D[1]*kerB[2]*kerB[22]+D[2]*kerB[1]*kerB[23]);
		Matrix[274] = detJ*(D[0]*kerB[1]*kerB[22]+D[2]*kerB[0]*kerB[21]+D[2]*kerB[2]*kerB[23]);
		Matrix[275] = detJ*(D[1]*kerB[0]*kerB[22]+D[2]*kerB[1]*kerB[21]);
		Matrix[276] = detJ*(D[0]*(kerB[23]*kerB[23])+D[2]*(kerB[21]*kerB[21])+D[2]*(kerB[22]*kerB[22]));
		Matrix[277] = detJ*(D[1]*kerB[22]*kerB[23]+D[2]*kerB[22]*kerB[23]);
		Matrix[278] = detJ*(D[1]*kerB[21]*kerB[23]+D[2]*kerB[21]*kerB[23]);
		Matrix[279] = detJ*(D[2]*kerB[18]*kerB[21]+D[0]*kerB[20]*kerB[23]+D[2]*kerB[19]*kerB[22]);
		Matrix[280] = detJ*(D[1]*kerB[19]*kerB[23]+D[2]*kerB[20]*kerB[22]);
		Matrix[281] = detJ*(D[1]*kerB[18]*kerB[23]+D[2]*kerB[20]*kerB[21]);
		Matrix[282] = detJ*(D[2]*kerB[15]*kerB[21]+D[0]*kerB[17]*kerB[23]+D[2]*kerB[16]*kerB[22]);
		Matrix[283] = detJ*(D[1]*kerB[16]*kerB[23]+D[2]*kerB[17]*kerB[22]);
		Matrix[284] = detJ*(D[1]*kerB[15]*kerB[23]+D[2]*kerB[17]*kerB[21]);
		Matrix[285] = detJ*(D[2]*kerB[12]*kerB[21]+D[0]*kerB[14]*kerB[23]+D[2]*kerB[13]*kerB[22]);
		Matrix[286] = detJ*(D[1]*kerB[13]*kerB[23]+D[2]*kerB[14]*kerB[22]);
		Matrix[287] = detJ*(D[1]*kerB[12]*kerB[23]+D[2]*kerB[14]*kerB[21]);
		Matrix[288] = detJ*(D[2]*kerB[9]*kerB[21]+D[0]*kerB[11]*kerB[23]+D[2]*kerB[10]*kerB[22]);
		Matrix[289] = detJ*(D[1]*kerB[10]*kerB[23]+D[2]*kerB[11]*kerB[22]);
		Matrix[290] = detJ*(D[1]*kerB[9]*kerB[23]+D[2]*kerB[11]*kerB[21]);
		Matrix[291] = detJ*(D[2]*kerB[6]*kerB[21]+D[0]*kerB[8]*kerB[23]+D[2]*kerB[7]*kerB[22]);
		Matrix[292] = detJ*(D[1]*kerB[7]*kerB[23]+D[2]*kerB[8]*kerB[22]);
		Matrix[293] = detJ*(D[1]*kerB[6]*kerB[23]+D[2]*kerB[8]*kerB[21]);
		Matrix[294] = detJ*(D[2]*kerB[3]*kerB[21]+D[0]*kerB[5]*kerB[23]+D[2]*kerB[4]*kerB[22]);
		Matrix[295] = detJ*(D[1]*kerB[4]*kerB[23]+D[2]*kerB[5]*kerB[22]);
		Matrix[296] = detJ*(D[1]*kerB[3]*kerB[23]+D[2]*kerB[5]*kerB[21]);
		Matrix[297] = detJ*(D[2]*kerB[0]*kerB[21]+D[0]*kerB[2]*kerB[23]+D[2]*kerB[1]*kerB[22]);
		Matrix[298] = detJ*(D[1]*kerB[1]*kerB[23]+D[2]*kerB[2]*kerB[22]);
		Matrix[299] = detJ*(D[1]*kerB[0]*kerB[23]+D[2]*kerB[2]*kerB[21]);
	}
}
	

//	Calculate element stress 
void CHex::ElementStress(double* stressHex, double* Displacement)
{
	// Get nodal displacements [LM can be used here]
	double Disp[24];
	for (unsigned int i = 0; i < 24; i++)
	{
		if (LocationMatrix[i])
			Disp[i] = Displacement[LocationMatrix[i]-1];
		else
			Disp[i] = 0.0;
	}

	// Construct constitutive matrix
	CHexMaterial* material = dynamic_cast<CHexMaterial*>(ElementMaterial);	// Pointer to material of the element
	double v = material->nu;
	double k = material->E * (1-v)/(1+v)/(1-2*v);
	double D[3];
	D[0] = k;
	D[1] = k * v / (1 - v);
	D[2] = k * (1 - 2 * v) / 2.0 / (1 - v);

	// Construct coordinate matrix
	double COORXYZ[24];
	for (unsigned int i=0;i<8;i++)
	{
		for (unsigned int j=0;j<3;j++)
		{
			COORXYZ[3 * i + j]=nodes[i]->XYZ[j];
		}
	}

	// Construct Jacobi matrix
	const double xi8[8] = { 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 };
	const double eta8[8] = { -0.577350269189626 , 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,0.577350269189626 , 0.577350269189626 ,-0.577350269189626 };
	const double zeta8[8] = { -0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 , 0.577350269189626 , 0.577350269189626 };

	double stressXYZ[6][8];	// 8 gauss points, 6 stress components
	for (unsigned p=0;p<8;p++)
	{
		double xi   = xi8[p];
		double eta  = eta8[p];
		double zeta = zeta8[p];

		double GN[12];
		GN[0] = (1-eta)*(1-zeta);
		GN[1] = (1+eta)*(1-zeta);
		GN[2] = (1-eta)*(1+zeta);
		GN[3] = (1+eta)*(1+zeta);
		GN[4] = (1+xi)*(1-zeta);
		GN[5] = (1-xi)*(1-zeta);
		GN[6] = (1+xi)*(1+zeta);
		GN[7] = (1-xi)*(1+zeta);
		GN[8] = (1+xi)*(1-eta);
		GN[9] = (1+xi)*(1+eta);
		GN[10] = (1-xi)*(1+eta);
		GN[11] = (1-xi)*(1-eta);

		double J[9];
		J[0] = COORXYZ[0]*GN[0]+COORXYZ[3]*GN[1]-COORXYZ[6]*GN[1]-COORXYZ[9]*GN[0]+COORXYZ[12]*GN[2]+COORXYZ[15]*GN[3]-COORXYZ[18]*GN[3]-COORXYZ[21]*GN[2];
		J[1] = COORXYZ[0]*GN[4]-COORXYZ[3]*GN[4]+COORXYZ[6]*GN[5]-COORXYZ[9]*GN[5]-COORXYZ[12]*GN[6]+COORXYZ[15]*GN[6]+COORXYZ[18]*GN[7]-COORXYZ[21]*GN[7];
		J[2] = -COORXYZ[0]*GN[8]-COORXYZ[3]*GN[9]-COORXYZ[6]*GN[10]-COORXYZ[9]*GN[11]+COORXYZ[12]*GN[8]+COORXYZ[15]*GN[9]+COORXYZ[18]*GN[10]+COORXYZ[21]*GN[11];
		J[3] = COORXYZ[1]*GN[0]+COORXYZ[4]*GN[1]-COORXYZ[7]*GN[1]-COORXYZ[10]*GN[0]+COORXYZ[13]*GN[2]+COORXYZ[16]*GN[3]-COORXYZ[19]*GN[3]-COORXYZ[22]*GN[2];
		J[4] = COORXYZ[1]*GN[4]-COORXYZ[4]*GN[4]+COORXYZ[7]*GN[5]-COORXYZ[10]*GN[5]-COORXYZ[13]*GN[6]+COORXYZ[16]*GN[6]+COORXYZ[19]*GN[7]-COORXYZ[22]*GN[7];
		J[5] = -COORXYZ[1]*GN[8]-COORXYZ[4]*GN[9]-COORXYZ[7]*GN[10]-COORXYZ[10]*GN[11]+COORXYZ[13]*GN[8]+COORXYZ[16]*GN[9]+COORXYZ[19]*GN[10]+COORXYZ[22]*GN[11];
		J[6] = COORXYZ[2]*GN[0]+COORXYZ[5]*GN[1]-COORXYZ[8]*GN[1]-COORXYZ[11]*GN[0]+COORXYZ[14]*GN[2]+COORXYZ[17]*GN[3]-COORXYZ[20]*GN[3]-COORXYZ[23]*GN[2];
		J[7] = COORXYZ[2]*GN[4]-COORXYZ[5]*GN[4]+COORXYZ[8]*GN[5]-COORXYZ[11]*GN[5]-COORXYZ[14]*GN[6]+COORXYZ[17]*GN[6]+COORXYZ[20]*GN[7]-COORXYZ[23]*GN[7];
		J[8] = -COORXYZ[2]*GN[8]-COORXYZ[5]*GN[9]-COORXYZ[8]*GN[10]-COORXYZ[11]*GN[11]+COORXYZ[14]*GN[8]+COORXYZ[17]*GN[9]+COORXYZ[20]*GN[10]+COORXYZ[23]*GN[11];

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
		kerB[0] = GN[0]*invJ[0]+GN[4]*invJ[3]-GN[8]*invJ[6];
		kerB[1] = GN[0]*invJ[1]+GN[4]*invJ[4]-GN[8]*invJ[7];
		kerB[2] = GN[0]*invJ[2]+GN[4]*invJ[5]-GN[8]*invJ[8];
		kerB[3] = GN[1]*invJ[0]-GN[4]*invJ[3]-GN[9]*invJ[6];
		kerB[4] = GN[1]*invJ[1]-GN[4]*invJ[4]-GN[9]*invJ[7];
		kerB[5] = GN[1]*invJ[2]-GN[4]*invJ[5]-GN[9]*invJ[8];
		kerB[6] = -GN[1]*invJ[0]+GN[5]*invJ[3]-GN[10]*invJ[6];
		kerB[7] = -GN[1]*invJ[1]+GN[5]*invJ[4]-GN[10]*invJ[7];
		kerB[8] = -GN[1]*invJ[2]+GN[5]*invJ[5]-GN[10]*invJ[8];
		kerB[9] = -GN[0]*invJ[0]-GN[5]*invJ[3]-GN[11]*invJ[6];
		kerB[10] = -GN[0]*invJ[1]-GN[5]*invJ[4]-GN[11]*invJ[7];
		kerB[11] = -GN[0]*invJ[2]-GN[5]*invJ[5]-GN[11]*invJ[8];
		kerB[12] = GN[2]*invJ[0]-GN[6]*invJ[3]+GN[8]*invJ[6];
		kerB[13] = GN[2]*invJ[1]-GN[6]*invJ[4]+GN[8]*invJ[7];
		kerB[14] = GN[2]*invJ[2]-GN[6]*invJ[5]+GN[8]*invJ[8];
		kerB[15] = GN[3]*invJ[0]+GN[6]*invJ[3]+GN[9]*invJ[6];
		kerB[16] = GN[3]*invJ[1]+GN[6]*invJ[4]+GN[9]*invJ[7];
		kerB[17] = GN[3]*invJ[2]+GN[6]*invJ[5]+GN[9]*invJ[8];
		kerB[18] = -GN[3]*invJ[0]+GN[7]*invJ[3]+GN[10]*invJ[6];
		kerB[19] = -GN[3]*invJ[1]+GN[7]*invJ[4]+GN[10]*invJ[7];
		kerB[20] = -GN[3]*invJ[2]+GN[7]*invJ[5]+GN[10]*invJ[8];
		kerB[21] = -GN[2]*invJ[0]-GN[7]*invJ[3]+GN[11]*invJ[6];
		kerB[22] = -GN[2]*invJ[1]-GN[7]*invJ[4]+GN[11]*invJ[7];
		kerB[23] = -GN[2]*invJ[2]-GN[7]*invJ[5]+GN[11]*invJ[8];

		stressXYZ[0][p] = D[1] * Disp[1] * kerB[1] + D[1] * Disp[2] * kerB[2] + D[0] * Disp[3] * kerB[3] + D[1] * Disp[4] * kerB[4] + D[1] * Disp[5] * kerB[5] + D[0] * Disp[6] * kerB[6] + D[1] * Disp[7] * kerB[7] + D[1] * Disp[8] * kerB[8] + D[0] * Disp[9] * kerB[9] + D[1] * Disp[10] * kerB[10] + D[1] * Disp[11] * kerB[11] + D[0] * Disp[12] * kerB[12] + D[1] * Disp[13] * kerB[13] + D[1] * Disp[14] * kerB[14] + D[0] * Disp[15] * kerB[15] + D[1] * Disp[16] * kerB[16] + D[1] * Disp[17] * kerB[17] + D[0] * Disp[18] * kerB[18] + D[1] * Disp[19] * kerB[19] + D[1] * Disp[20] * kerB[20] + D[0] * Disp[21] * kerB[21] + D[1] * Disp[22] * kerB[22] + D[1] * Disp[23] * kerB[23] + D[0] * Disp[0] * kerB[0];
		stressXYZ[1][p] = D[0] * Disp[1] * kerB[1] + D[1] * Disp[2] * kerB[2] + D[1] * Disp[3] * kerB[3] + D[0] * Disp[4] * kerB[4] + D[1] * Disp[5] * kerB[5] + D[1] * Disp[6] * kerB[6] + D[0] * Disp[7] * kerB[7] + D[1] * Disp[8] * kerB[8] + D[1] * Disp[9] * kerB[9] + D[0] * Disp[10] * kerB[10] + D[1] * Disp[11] * kerB[11] + D[1] * Disp[12] * kerB[12] + D[0] * Disp[13] * kerB[13] + D[1] * Disp[14] * kerB[14] + D[1] * Disp[15] * kerB[15] + D[0] * Disp[16] * kerB[16] + D[1] * Disp[17] * kerB[17] + D[1] * Disp[18] * kerB[18] + D[0] * Disp[19] * kerB[19] + D[1] * Disp[20] * kerB[20] + D[1] * Disp[21] * kerB[21] + D[0] * Disp[22] * kerB[22] + D[1] * Disp[23] * kerB[23] + D[1] * Disp[0] * kerB[0];
		stressXYZ[2][p] = D[1] * Disp[1] * kerB[1] + D[0] * Disp[2] * kerB[2] + D[1] * Disp[3] * kerB[3] + D[1] * Disp[4] * kerB[4] + D[0] * Disp[5] * kerB[5] + D[1] * Disp[6] * kerB[6] + D[1] * Disp[7] * kerB[7] + D[0] * Disp[8] * kerB[8] + D[1] * Disp[9] * kerB[9] + D[1] * Disp[10] * kerB[10] + D[0] * Disp[11] * kerB[11] + D[1] * Disp[12] * kerB[12] + D[1] * Disp[13] * kerB[13] + D[0] * Disp[14] * kerB[14] + D[1] * Disp[15] * kerB[15] + D[1] * Disp[16] * kerB[16] + D[0] * Disp[17] * kerB[17] + D[1] * Disp[18] * kerB[18] + D[1] * Disp[19] * kerB[19] + D[0] * Disp[20] * kerB[20] + D[1] * Disp[21] * kerB[21] + D[1] * Disp[22] * kerB[22] + D[0] * Disp[23] * kerB[23] + D[1] * Disp[0] * kerB[0];
		stressXYZ[3][p] = D[2] * Disp[1] * kerB[0] + D[2] * Disp[3] * kerB[4] + D[2] * Disp[4] * kerB[3] + D[2] * Disp[6] * kerB[7] + D[2] * Disp[7] * kerB[6] + D[2] * Disp[9] * kerB[10] + D[2] * Disp[10] * kerB[9] + D[2] * Disp[12] * kerB[13] + D[2] * Disp[13] * kerB[12] + D[2] * Disp[15] * kerB[16] + D[2] * Disp[16] * kerB[15] + D[2] * Disp[18] * kerB[19] + D[2] * Disp[19] * kerB[20] + D[2] * Disp[21] * kerB[22] + D[2] * Disp[22] * kerB[21] + D[2] * Disp[0] * kerB[1];
		stressXYZ[4][p] = D[2] * Disp[1] * kerB[2] + D[2] * Disp[2] * kerB[1] + D[2] * Disp[4] * kerB[5] + D[2] * Disp[5] * kerB[4] + D[2] * Disp[7] * kerB[8] + D[2] * Disp[8] * kerB[7] + D[2] * Disp[10] * kerB[11] + D[2] * Disp[11] * kerB[10] + D[2] * Disp[13] * kerB[14] + D[2] * Disp[14] * kerB[13] + D[2] * Disp[16] * kerB[17] + D[2] * Disp[17] * kerB[16] + D[2] * Disp[19] * kerB[20] + D[2] * Disp[20] * kerB[19] + D[2] * Disp[22] * kerB[23] + D[2] * Disp[23] * kerB[22];
		stressXYZ[5][p] = D[2] * Disp[2] * kerB[0] + D[2] * Disp[3] * kerB[5] + D[2] * Disp[5] * kerB[3] + D[2] * Disp[6] * kerB[8] + D[2] * Disp[8] * kerB[6] + D[2] * Disp[9] * kerB[11] + D[2] * Disp[11] * kerB[9] + D[2] * Disp[12] * kerB[14] + D[2] * Disp[14] * kerB[12] + D[2] * Disp[15] * kerB[17] + D[2] * Disp[17] * kerB[15] + D[2] * Disp[18] * kerB[20] + D[2] * Disp[20] * kerB[18] + D[2] * Disp[21] * kerB[23] + D[2] * Disp[23] * kerB[21] + D[2] * Disp[0] * kerB[2];
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
			stressHex[6 * j + i] = recovery[i];
		}
	}
}
