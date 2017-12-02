#include "Elements/Plate.h"

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

//	Constructor
CPlate::CPlate()
{
	NEN = 4;	// Each element has 2 nodes
	nodes = new CNode*[NEN];
	
	ND = 24;
    LocationMatrix = new unsigned int[24];

	ElementMaterial = nullptr;
}

//	Desconstructor
CPlate::~CPlate()
{
}

//	Read element data from stream Input
bool CPlate::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4;	// node numbers in counter-clock sequence

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial = dynamic_cast<CPlateMaterial*>(MaterialSets) + MSet - 1;
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	nodes[3] = &NodeList[N4 - 1];
	/*thinking about how to add rotating degrees of freedom*/

	/*the codes next is a proof that this code is not well developed. Other situations will be considered later.*/
	double edges[6];
	edges[0]=nodes[1]->XYZ[0]-nodes[0]->XYZ[0];
	edges[1]=nodes[1]->XYZ[1]-nodes[0]->XYZ[1];
	edges[2]=nodes[2]->XYZ[0]-nodes[0]->XYZ[0];
	edges[3]=nodes[2]->XYZ[1]-nodes[0]->XYZ[1];
	edges[4]=nodes[3]->XYZ[0]-nodes[0]->XYZ[0];
	edges[5]=nodes[3]->XYZ[1]-nodes[0]->XYZ[1];
	
	return true;
}

//Write plate element data
void CPlate::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber 
		   << setw(9) << nodes[1]->NodeNumber << setw(9) << nodes[2]->NodeNumber 
		   << setw(9) << nodes[3]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CPlate::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 6; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
	//This plate is only a test one. The ultimate plate would be released in shell.
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CPlate::SizeOfStiffnessMatrix() { return 300; }

/*it makes me think of another question that when should the rotation dof turn to zero*/ 


//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CPlate::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	
	double xdir[3];
	double ydir[3];
	double zdir[3];
	xdir[0]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0]);
	xdir[1]=(nodes[1]->XYZ[1]-nodes[0]->XYZ[1]);
	xdir[2]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2]);
	ydir[0]=(nodes[3]->XYZ[0]-nodes[0]->XYZ[0]);
	ydir[1]=(nodes[3]->XYZ[1]-nodes[0]->XYZ[1]);
	ydir[2]=(nodes[3]->XYZ[2]-nodes[0]->XYZ[2]);
	double LX2=xdir[0]*xdir[0]+xdir[1]*xdir[1]+xdir[2]*xdir[2];
	double LX=sqrt(LX2);
	double LY2=ydir[0]*ydir[0]+ydir[1]*ydir[1]+ydir[2]*ydir[2];
	double LY=sqrt(LY2);
	xdir[0]=xdir[0]/LX;
	xdir[1]=xdir[1]/LX;
	xdir[2]=xdir[2]/LX;
	ydir[0]=ydir[0]/LY;
	ydir[1]=ydir[1]/LY;
	ydir[2]=ydir[2]/LY;
	zdir[0]=xdir[1]*ydir[2]-xdir[2]*ydir[1];
	zdir[1]=xdir[2]*ydir[0]-xdir[0]*ydir[2];
	zdir[2]=xdir[0]*ydir[1]-xdir[1]*ydir[0];


	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial);	// Pointer to material of the element

	double nu=material->nu;
	const double k= material->E*material->h*material->h*material->h*0.0833333333333333/(1-nu*nu)*LX*LY/4;

	double xpsi=LX/2;
	double yeta=LY/2;

#ifdef _DEBUG_
	cout<<"Jacobian"<<setw(20)<<LX<<setw(20)<<LY;
#endif
	double MatrixLit[78];
	MatrixLit[0] = k*(1/(xpsi*xpsi*xpsi*xpsi) - ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) - (xpsi*xpsi*xpsi*xpsi))/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[1] = -k*((8*(nu/2 - 1/2))/(15*(xpsi*xpsi)) - 4/(3*(yeta*yeta)));
	MatrixLit[2] = k*(1/(yeta*yeta*yeta) + ((2*nu)/5 + 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[3] = -k*((4*(nu - 1))/(15*(yeta*yeta)) - 4/(3*(xpsi*xpsi)));
	MatrixLit[4] = -(k*nu)/(xpsi*yeta);
	MatrixLit[5] = -k*(1/(xpsi*xpsi*xpsi) + ((2*nu)/5 + 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[6] = k*(1/(xpsi*xpsi*xpsi*xpsi) - ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) - (xpsi*xpsi*xpsi*xpsi))/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[7] = k*(1/(xpsi*xpsi*xpsi) - (nu/10 - 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[8] = k*(1/(2*(yeta*yeta*yeta)) - ((2*nu)/5 + 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[9] = -k*(1/(xpsi*xpsi*xpsi*xpsi) - ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) + (xpsi*xpsi*xpsi*xpsi)/2)/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[10] = -k*((8*(nu/2 - 1/2))/(15*(xpsi*xpsi)) - 4/(3*(yeta*yeta)));
	MatrixLit[11] = k*(1/(yeta*yeta*yeta) + ((2*nu)/5 + 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[12] = 0;
	MatrixLit[13] = k*((8*(nu/2 - 1/2))/(15*(xpsi*xpsi)) + 2/(3*(yeta*yeta)));
	MatrixLit[14] = k*(1/(2*(yeta*yeta*yeta)) - ((2*nu)/5 + 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[15] = -k*((4*(nu - 1))/(15*(yeta*yeta)) - 4/(3*(xpsi*xpsi)));
	MatrixLit[16] = (k*nu)/(xpsi*yeta);
	MatrixLit[17] = k*(1/(xpsi*xpsi*xpsi) + ((2*nu)/5 + 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[18] = k*((3*nu - 3)/(45*(yeta*yeta)) + 2/(3*(xpsi*xpsi)));
	MatrixLit[19] = 0;
	MatrixLit[20] = -k*(1/(xpsi*xpsi*xpsi) - (nu/10 - 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[21] = k*(1/(xpsi*xpsi*xpsi*xpsi) - ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) - (xpsi*xpsi*xpsi*xpsi))/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[22] = k*(1/(2*(xpsi*xpsi*xpsi)) - ((2*nu)/5 + 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[23] = -k*(1/(yeta*yeta*yeta) - (nu/10 - 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[24] = -k*(1/(yeta*yeta*yeta*yeta) - ((xpsi*xpsi)*((nu*(yeta*yeta))/5 - (7*(yeta*yeta))/10) + (yeta*yeta*yeta*yeta)/2)/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[25] = (k*((xpsi*xpsi)*(nu/10 - 1/10) + (yeta*yeta)/2))/((xpsi*xpsi*xpsi)*(yeta*yeta));
	MatrixLit[26] = -k*(1/(2*(yeta*yeta*yeta)) + (nu/10 - 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[27] = -k*(1/(2*(xpsi*xpsi*xpsi*xpsi)) + ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) + (xpsi*xpsi*xpsi*xpsi)/2)/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[28] = -k*((8*(nu/2 - 1/2))/(15*(xpsi*xpsi)) - 4/(3*(yeta*yeta)));
	MatrixLit[29] = -k*(1/(yeta*yeta*yeta) + ((2*nu)/5 + 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[30] = 0;
	MatrixLit[31] = k*((nu/15 - 1/15)/(xpsi*xpsi) + 2/(3*(yeta*yeta)));
	MatrixLit[32] = k*(1/(yeta*yeta*yeta) - (nu/10 - 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[33] = 0;
	MatrixLit[34] = -k*((nu/15 - 1/15)/(xpsi*xpsi) - 1/(3*(yeta*yeta)));
	MatrixLit[35] = (k*((yeta*yeta)*(nu/10 - 1/10) + (xpsi*xpsi)/2))/((xpsi*xpsi)*(yeta*yeta*yeta));
	MatrixLit[36] = -k*((4*(nu - 1))/(15*(yeta*yeta)) - 4/(3*(xpsi*xpsi)));
	MatrixLit[37] = -(k*nu)/(xpsi*yeta);
	MatrixLit[38] = k*(1/(xpsi*xpsi*xpsi) + ((2*nu)/5 + 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[39] = k*((4*(nu - 1))/(15*(yeta*yeta)) + 2/(3*(xpsi*xpsi)));
	MatrixLit[40] = 0;
	MatrixLit[41] = k*(1/(2*(xpsi*xpsi*xpsi)) - ((2*nu)/5 + 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[42] = -k*((3*nu - 3)/(45*(yeta*yeta)) - 1/(3*(xpsi*xpsi)));
	MatrixLit[43] = 0;
	MatrixLit[44] = -k*(1/(2*(xpsi*xpsi*xpsi)) + (nu/10 - 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[45] = k*(1/(xpsi*xpsi*xpsi*xpsi) - ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) - (xpsi*xpsi*xpsi*xpsi))/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[46] = -k*(1/(xpsi*xpsi*xpsi) - (nu/10 - 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[47] = (k*((yeta*yeta)*((2*nu)/5 + 1/10) - (xpsi*xpsi)/2))/((xpsi*xpsi)*(yeta*yeta*yeta));
	MatrixLit[48] = -k*(1/(xpsi*xpsi*xpsi*xpsi) - ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) + (xpsi*xpsi*xpsi*xpsi)/2)/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[49] = -k*(1/(2*(xpsi*xpsi*xpsi)) + (nu/10 - 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[50] = -k*(1/(2*(yeta*yeta*yeta)) + (nu/10 - 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[51] = -k*(1/(2*(xpsi*xpsi*xpsi*xpsi)) + ((yeta*yeta)*((nu*(xpsi*xpsi))/5 - (7*(xpsi*xpsi))/10) + (xpsi*xpsi*xpsi*xpsi)/2)/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[52] = (k*((xpsi*xpsi)*((2*nu)/5 + 1/10) - (yeta*yeta)/2))/((xpsi*xpsi*xpsi)*(yeta*yeta));
	MatrixLit[53] = -k*(1/(yeta*yeta*yeta) - (nu/10 - 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[54] = -k*(1/(yeta*yeta*yeta*yeta) - ((xpsi*xpsi)*((nu*(yeta*yeta))/5 - (7*(yeta*yeta))/10) + (yeta*yeta*yeta*yeta)/2)/((xpsi*xpsi*xpsi*xpsi)*(yeta*yeta*yeta*yeta)));
	MatrixLit[55] = -k*((8*(nu/2 - 1/2))/(15*(xpsi*xpsi)) - 4/(3*(yeta*yeta)));
	MatrixLit[56] = -k*(1/(yeta*yeta*yeta) + ((2*nu)/5 + 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[57] = 0;
	MatrixLit[58] = k*((8*(nu/2 - 1/2))/(15*(xpsi*xpsi)) + 2/(3*(yeta*yeta)));
	MatrixLit[59] = (k*((yeta*yeta)*((2*nu)/5 + 1/10) - (xpsi*xpsi)/2))/((xpsi*xpsi)*(yeta*yeta*yeta));
	MatrixLit[60] = 0;
	MatrixLit[61] = -k*((nu/15 - 1/15)/(xpsi*xpsi) - 1/(3*(yeta*yeta)));
	MatrixLit[62] = (k*((yeta*yeta)*(nu/10 - 1/10) + (xpsi*xpsi)/2))/((xpsi*xpsi)*(yeta*yeta*yeta));
	MatrixLit[63] = 0;
	MatrixLit[64] = k*((nu/15 - 1/15)/(xpsi*xpsi) + 2/(3*(yeta*yeta)));
	MatrixLit[65] = k*(1/(yeta*yeta*yeta) - (nu/10 - 1/10)/((xpsi*xpsi)*yeta));
	MatrixLit[66] = -k*((4*(nu - 1))/(15*(yeta*yeta)) - 4/(3*(xpsi*xpsi)));
	MatrixLit[67] = (k*nu)/(xpsi*yeta);
	MatrixLit[68] = -k*(1/(xpsi*xpsi*xpsi) + ((2*nu)/5 + 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[69] = k*((3*nu - 3)/(45*(yeta*yeta)) + 2/(3*(xpsi*xpsi)));
	MatrixLit[70] = 0;
	MatrixLit[71] = k*(1/(xpsi*xpsi*xpsi) - (nu/10 - 1/10)/(xpsi*(yeta*yeta)));
	MatrixLit[72] = -k*((3*nu - 3)/(45*(yeta*yeta)) - 1/(3*(xpsi*xpsi)));
	MatrixLit[73] = 0;
	MatrixLit[74] = (k*((xpsi*xpsi)*(nu/10 - 1/10) + (yeta*yeta)/2))/((xpsi*xpsi*xpsi)*(yeta*yeta));
	MatrixLit[75] = k*((4*(nu - 1))/(15*(yeta*yeta)) + 2/(3*(xpsi*xpsi)));
	MatrixLit[76] = 0;
	MatrixLit[77] = (k*((xpsi*xpsi)*((2*nu)/5 + 1/10) - (yeta*yeta)/2))/((xpsi*xpsi*xpsi)*(yeta*yeta));

#ifdef _DEBUG_
	for (unsigned int i=0;i<78;++i)
		{
			cout<<"Matrix element "<<i<<" is "<<MatrixLit[i]<<endl;
	};
#endif

	Matrix[0] = MatrixLit[0] * zdir[0]*zdir[0];
	Matrix[1] = MatrixLit[0] * zdir[1]*zdir[1];
	Matrix[2] = MatrixLit[0] * zdir[0]*zdir[1];
	Matrix[3] = MatrixLit[0] * zdir[2]*zdir[2];
	Matrix[4] = MatrixLit[0] * zdir[1]*zdir[2];
	Matrix[5] = MatrixLit[0] * zdir[0]*zdir[2];
	Matrix[6] = MatrixLit[1] * xdir[0]*xdir[0] + MatrixLit[4] * xdir[0]*ydir[0] + MatrixLit[4] * ydir[0]*xdir[0] + MatrixLit[3] * ydir[0]*ydir[0];
	Matrix[7] = MatrixLit[2] * zdir[2]*xdir[0] + MatrixLit[5] * zdir[2]*ydir[0];
	Matrix[8] = MatrixLit[2] * zdir[1]*xdir[0] + MatrixLit[5] * zdir[1]*ydir[0];
	Matrix[9] = MatrixLit[2] * zdir[0]*xdir[0] + MatrixLit[5] * zdir[0]*ydir[0];
	Matrix[10] = MatrixLit[1] * xdir[1]*xdir[1] + MatrixLit[4] * xdir[1]*ydir[1] + MatrixLit[4] * ydir[1]*xdir[1] + MatrixLit[3] * ydir[1]*ydir[1];
	Matrix[11] = MatrixLit[1] * xdir[0]*xdir[1] + MatrixLit[4] * xdir[0]*ydir[1] + MatrixLit[4] * ydir[0]*xdir[1] + MatrixLit[3] * ydir[0]*ydir[1];
	Matrix[12] = MatrixLit[2] * zdir[2]*xdir[1] + MatrixLit[5] * zdir[2]*ydir[1];
	Matrix[13] = MatrixLit[2] * zdir[1]*xdir[1] + MatrixLit[5] * zdir[1]*ydir[1];
	Matrix[14] = MatrixLit[2] * zdir[0]*xdir[1] + MatrixLit[5] * zdir[0]*ydir[1];
	Matrix[15] = MatrixLit[10] * xdir[2]*xdir[2] + MatrixLit[16] * xdir[2]*ydir[2] + MatrixLit[16] * ydir[2]*xdir[2] + MatrixLit[15] * ydir[2]*ydir[2];
	Matrix[16] = MatrixLit[13] * xdir[1]*xdir[2] + MatrixLit[19] * xdir[1]*ydir[2] + MatrixLit[12] * ydir[1]*xdir[2] + MatrixLit[18] * ydir[1]*ydir[2];
	Matrix[17] = MatrixLit[13] * xdir[0]*xdir[2] + MatrixLit[19] * xdir[0]*ydir[2] + MatrixLit[12] * ydir[0]*xdir[2] + MatrixLit[18] * ydir[0]*ydir[2];
	Matrix[18] = MatrixLit[14] * zdir[2]*xdir[2] + MatrixLit[5] * zdir[2]*ydir[2];
	Matrix[19] = MatrixLit[14] * zdir[1]*xdir[2] + MatrixLit[5] * zdir[1]*ydir[2];
	Matrix[20] = MatrixLit[14] * zdir[0]*xdir[2] + MatrixLit[5] * zdir[0]*ydir[2];
	Matrix[21] = MatrixLit[6] * zdir[0]*zdir[0];
	Matrix[22] = MatrixLit[11] * xdir[2]*zdir[0] + MatrixLit[17] * ydir[2]*zdir[0];
	Matrix[23] = MatrixLit[8] * xdir[1]*zdir[0] + MatrixLit[7] * ydir[1]*zdir[0];
	Matrix[24] = MatrixLit[8] * xdir[0]*zdir[0] + MatrixLit[7] * ydir[0]*zdir[0];
	Matrix[25] = MatrixLit[9] * zdir[2]*zdir[0];
	Matrix[26] = MatrixLit[9] * zdir[1]*zdir[0];
	Matrix[27] = MatrixLit[9] * zdir[0]*zdir[0];
	Matrix[28] = MatrixLit[6] * zdir[1]*zdir[1];
	Matrix[29] = MatrixLit[6] * zdir[0]*zdir[1];
	Matrix[30] = MatrixLit[11] * xdir[2]*zdir[1] + MatrixLit[17] * ydir[2]*zdir[1];
	Matrix[31] = MatrixLit[8] * xdir[1]*zdir[1] + MatrixLit[7] * ydir[1]*zdir[1];
	Matrix[32] = MatrixLit[8] * xdir[0]*zdir[1] + MatrixLit[7] * ydir[0]*zdir[1];
	Matrix[33] = MatrixLit[9] * zdir[2]*zdir[1];
	Matrix[34] = MatrixLit[9] * zdir[1]*zdir[1];
	Matrix[35] = MatrixLit[9] * zdir[0]*zdir[1];
	Matrix[36] = MatrixLit[6] * zdir[2]*zdir[2];
	Matrix[37] = MatrixLit[6] * zdir[1]*zdir[2];
	Matrix[38] = MatrixLit[6] * zdir[0]*zdir[2];
	Matrix[39] = MatrixLit[11] * xdir[2]*zdir[2] + MatrixLit[17] * ydir[2]*zdir[2];
	Matrix[40] = MatrixLit[8] * xdir[1]*zdir[2] + MatrixLit[7] * ydir[1]*zdir[2];
	Matrix[41] = MatrixLit[8] * xdir[0]*zdir[2] + MatrixLit[7] * ydir[0]*zdir[2];
	Matrix[42] = MatrixLit[9] * zdir[2]*zdir[2];
	Matrix[43] = MatrixLit[9] * zdir[1]*zdir[2];
	Matrix[44] = MatrixLit[9] * zdir[0]*zdir[2];
	Matrix[45] = MatrixLit[10] * xdir[0]*xdir[0] + MatrixLit[16] * xdir[0]*ydir[0] + MatrixLit[16] * ydir[0]*xdir[0] + MatrixLit[15] * ydir[0]*ydir[0];
	Matrix[46] = MatrixLit[11] * zdir[2]*xdir[0] + MatrixLit[17] * zdir[2]*ydir[0];
	Matrix[47] = MatrixLit[11] * zdir[1]*xdir[0] + MatrixLit[17] * zdir[1]*ydir[0];
	Matrix[48] = MatrixLit[11] * zdir[0]*xdir[0] + MatrixLit[17] * zdir[0]*ydir[0];
	Matrix[49] = MatrixLit[10] * xdir[2]*xdir[0] + MatrixLit[16] * xdir[2]*ydir[0] + MatrixLit[16] * ydir[2]*xdir[0] + MatrixLit[15] * ydir[2]*ydir[0];
	Matrix[50] = MatrixLit[13] * xdir[1]*xdir[0] + MatrixLit[19] * xdir[1]*ydir[0] + MatrixLit[12] * ydir[1]*xdir[0] + MatrixLit[18] * ydir[1]*ydir[0];
	Matrix[51] = MatrixLit[13] * xdir[0]*xdir[0] + MatrixLit[19] * xdir[0]*ydir[0] + MatrixLit[12] * ydir[0]*xdir[0] + MatrixLit[18] * ydir[0]*ydir[0];
	Matrix[52] = MatrixLit[14] * zdir[2]*xdir[0] + MatrixLit[5] * zdir[2]*ydir[0];
	Matrix[53] = MatrixLit[14] * zdir[1]*xdir[0] + MatrixLit[5] * zdir[1]*ydir[0];
	Matrix[54] = MatrixLit[14] * zdir[0]*xdir[0] + MatrixLit[5] * zdir[0]*ydir[0];
	Matrix[55] = MatrixLit[10] * xdir[1]*xdir[1] + MatrixLit[16] * xdir[1]*ydir[1] + MatrixLit[16] * ydir[1]*xdir[1] + MatrixLit[15] * ydir[1]*ydir[1];
	Matrix[56] = MatrixLit[10] * xdir[0]*xdir[1] + MatrixLit[16] * xdir[0]*ydir[1] + MatrixLit[16] * ydir[0]*xdir[1] + MatrixLit[15] * ydir[0]*ydir[1];
	Matrix[57] = MatrixLit[11] * zdir[2]*xdir[1] + MatrixLit[17] * zdir[2]*ydir[1];
	Matrix[58] = MatrixLit[11] * zdir[1]*xdir[1] + MatrixLit[17] * zdir[1]*ydir[1];
	Matrix[59] = MatrixLit[11] * zdir[0]*xdir[1] + MatrixLit[17] * zdir[0]*ydir[1];
	Matrix[60] = MatrixLit[10] * xdir[2]*xdir[1] + MatrixLit[16] * xdir[2]*ydir[1] + MatrixLit[16] * ydir[2]*xdir[1] + MatrixLit[15] * ydir[2]*ydir[1];
	Matrix[61] = MatrixLit[13] * xdir[1]*xdir[1] + MatrixLit[19] * xdir[1]*ydir[1] + MatrixLit[12] * ydir[1]*xdir[1] + MatrixLit[18] * ydir[1]*ydir[1];
	Matrix[62] = MatrixLit[13] * xdir[0]*xdir[1] + MatrixLit[19] * xdir[0]*ydir[1] + MatrixLit[12] * ydir[0]*xdir[1] + MatrixLit[18] * ydir[0]*ydir[1];
	Matrix[63] = MatrixLit[14] * zdir[2]*xdir[1] + MatrixLit[5] * zdir[2]*ydir[1];
	Matrix[64] = MatrixLit[14] * zdir[1]*xdir[1] + MatrixLit[5] * zdir[1]*ydir[1];
	Matrix[65] = MatrixLit[14] * zdir[0]*xdir[1] + MatrixLit[5] * zdir[0]*ydir[1];
	Matrix[66] = MatrixLit[28] * xdir[2]*xdir[2] + MatrixLit[37] * xdir[2]*ydir[2] + MatrixLit[37] * ydir[2]*xdir[2] + MatrixLit[36] * ydir[2]*ydir[2];
	Matrix[67] = MatrixLit[31] * xdir[1]*xdir[2] + MatrixLit[40] * xdir[1]*ydir[2] + MatrixLit[30] * ydir[1]*xdir[2] + MatrixLit[39] * ydir[1]*ydir[2];
	Matrix[68] = MatrixLit[31] * xdir[0]*xdir[2] + MatrixLit[40] * xdir[0]*ydir[2] + MatrixLit[30] * ydir[0]*xdir[2] + MatrixLit[39] * ydir[0]*ydir[2];
	Matrix[69] = MatrixLit[32] * zdir[2]*xdir[2] + MatrixLit[17] * zdir[2]*ydir[2];
	Matrix[70] = MatrixLit[32] * zdir[1]*xdir[2] + MatrixLit[17] * zdir[1]*ydir[2];
	Matrix[71] = MatrixLit[32] * zdir[0]*xdir[2] + MatrixLit[17] * zdir[0]*ydir[2];
	Matrix[72] = MatrixLit[31] * xdir[2]*xdir[2] + MatrixLit[40] * xdir[2]*ydir[2] + MatrixLit[30] * ydir[2]*xdir[2] + MatrixLit[39] * ydir[2]*ydir[2];
	Matrix[73] = MatrixLit[34] * xdir[1]*xdir[2] + MatrixLit[43] * xdir[1]*ydir[2] + MatrixLit[33] * ydir[1]*xdir[2] + MatrixLit[42] * ydir[1]*ydir[2];
	Matrix[74] = MatrixLit[34] * xdir[0]*xdir[2] + MatrixLit[43] * xdir[0]*ydir[2] + MatrixLit[33] * ydir[0]*xdir[2] + MatrixLit[42] * ydir[0]*ydir[2];
	Matrix[75] = MatrixLit[35] * zdir[2]*xdir[2] + MatrixLit[5] * zdir[2]*ydir[2];
	Matrix[76] = MatrixLit[35] * zdir[1]*xdir[2] + MatrixLit[5] * zdir[1]*ydir[2];
	Matrix[77] = MatrixLit[35] * zdir[0]*xdir[2] + MatrixLit[5] * zdir[0]*ydir[2];
	Matrix[78] = MatrixLit[21] * zdir[0]*zdir[0];
	Matrix[79] = MatrixLit[29] * xdir[2]*zdir[0] + MatrixLit[38] * ydir[2]*zdir[0];
	Matrix[80] = MatrixLit[23] * xdir[1]*zdir[0] + MatrixLit[22] * ydir[1]*zdir[0];
	Matrix[81] = MatrixLit[23] * xdir[0]*zdir[0] + MatrixLit[22] * ydir[0]*zdir[0];
	Matrix[82] = MatrixLit[24] * zdir[2]*zdir[0];
	Matrix[83] = MatrixLit[24] * zdir[1]*zdir[0];
	Matrix[84] = MatrixLit[24] * zdir[0]*zdir[0];
	Matrix[85] = MatrixLit[23] * xdir[2]*zdir[0] + MatrixLit[22] * ydir[2]*zdir[0];
	Matrix[86] = MatrixLit[26] * xdir[1]*zdir[0] + MatrixLit[25] * ydir[1]*zdir[0];
	Matrix[87] = MatrixLit[26] * xdir[0]*zdir[0] + MatrixLit[25] * ydir[0]*zdir[0];
	Matrix[88] = MatrixLit[27] * zdir[2]*zdir[0];
	Matrix[89] = MatrixLit[27] * zdir[1]*zdir[0];
	Matrix[90] = MatrixLit[27] * zdir[0]*zdir[0];
	Matrix[91] = MatrixLit[21] * zdir[1]*zdir[1];
	Matrix[92] = MatrixLit[21] * zdir[0]*zdir[1];
	Matrix[93] = MatrixLit[29] * xdir[2]*zdir[1] + MatrixLit[38] * ydir[2]*zdir[1];
	Matrix[94] = MatrixLit[23] * xdir[1]*zdir[1] + MatrixLit[22] * ydir[1]*zdir[1];
	Matrix[95] = MatrixLit[23] * xdir[0]*zdir[1] + MatrixLit[22] * ydir[0]*zdir[1];
	Matrix[96] = MatrixLit[24] * zdir[2]*zdir[1];
	Matrix[97] = MatrixLit[24] * zdir[1]*zdir[1];
	Matrix[98] = MatrixLit[24] * zdir[0]*zdir[1];
	Matrix[99] = MatrixLit[23] * xdir[2]*zdir[1] + MatrixLit[22] * ydir[2]*zdir[1];
	Matrix[100] = MatrixLit[26] * xdir[1]*zdir[1] + MatrixLit[25] * ydir[1]*zdir[1];
	Matrix[101] = MatrixLit[26] * xdir[0]*zdir[1] + MatrixLit[25] * ydir[0]*zdir[1];
	Matrix[102] = MatrixLit[27] * zdir[2]*zdir[1];
	Matrix[103] = MatrixLit[27] * zdir[1]*zdir[1];
	Matrix[104] = MatrixLit[27] * zdir[0]*zdir[1];
	Matrix[105] = MatrixLit[21] * zdir[2]*zdir[2];
	Matrix[106] = MatrixLit[21] * zdir[1]*zdir[2];
	Matrix[107] = MatrixLit[21] * zdir[0]*zdir[2];
	Matrix[108] = MatrixLit[29] * xdir[2]*zdir[2] + MatrixLit[38] * ydir[2]*zdir[2];
	Matrix[109] = MatrixLit[23] * xdir[1]*zdir[2] + MatrixLit[22] * ydir[1]*zdir[2];
	Matrix[110] = MatrixLit[23] * xdir[0]*zdir[2] + MatrixLit[22] * ydir[0]*zdir[2];
	Matrix[111] = MatrixLit[24] * zdir[2]*zdir[2];
	Matrix[112] = MatrixLit[24] * zdir[1]*zdir[2];
	Matrix[113] = MatrixLit[24] * zdir[0]*zdir[2];
	Matrix[114] = MatrixLit[23] * xdir[2]*zdir[2] + MatrixLit[22] * ydir[2]*zdir[2];
	Matrix[115] = MatrixLit[26] * xdir[1]*zdir[2] + MatrixLit[25] * ydir[1]*zdir[2];
	Matrix[116] = MatrixLit[26] * xdir[0]*zdir[2] + MatrixLit[25] * ydir[0]*zdir[2];
	Matrix[117] = MatrixLit[27] * zdir[2]*zdir[2];
	Matrix[118] = MatrixLit[27] * zdir[1]*zdir[2];
	Matrix[119] = MatrixLit[27] * zdir[0]*zdir[2];
	Matrix[120] = MatrixLit[28] * xdir[0]*xdir[0] + MatrixLit[37] * xdir[0]*ydir[0] + MatrixLit[37] * ydir[0]*xdir[0] + MatrixLit[36] * ydir[0]*ydir[0];
	Matrix[121] = MatrixLit[29] * zdir[2]*xdir[0] + MatrixLit[38] * zdir[2]*ydir[0];
	Matrix[122] = MatrixLit[29] * zdir[1]*xdir[0] + MatrixLit[38] * zdir[1]*ydir[0];
	Matrix[123] = MatrixLit[29] * zdir[0]*xdir[0] + MatrixLit[38] * zdir[0]*ydir[0];
	Matrix[124] = MatrixLit[28] * xdir[2]*xdir[0] + MatrixLit[37] * xdir[2]*ydir[0] + MatrixLit[37] * ydir[2]*xdir[0] + MatrixLit[36] * ydir[2]*ydir[0];
	Matrix[125] = MatrixLit[31] * xdir[1]*xdir[0] + MatrixLit[40] * xdir[1]*ydir[0] + MatrixLit[30] * ydir[1]*xdir[0] + MatrixLit[39] * ydir[1]*ydir[0];
	Matrix[126] = MatrixLit[31] * xdir[0]*xdir[0] + MatrixLit[40] * xdir[0]*ydir[0] + MatrixLit[30] * ydir[0]*xdir[0] + MatrixLit[39] * ydir[0]*ydir[0];
	Matrix[127] = MatrixLit[32] * zdir[2]*xdir[0] + MatrixLit[17] * zdir[2]*ydir[0];
	Matrix[128] = MatrixLit[32] * zdir[1]*xdir[0] + MatrixLit[17] * zdir[1]*ydir[0];
	Matrix[129] = MatrixLit[32] * zdir[0]*xdir[0] + MatrixLit[17] * zdir[0]*ydir[0];
	Matrix[130] = MatrixLit[31] * xdir[2]*xdir[0] + MatrixLit[40] * xdir[2]*ydir[0] + MatrixLit[30] * ydir[2]*xdir[0] + MatrixLit[39] * ydir[2]*ydir[0];
	Matrix[131] = MatrixLit[34] * xdir[1]*xdir[0] + MatrixLit[43] * xdir[1]*ydir[0] + MatrixLit[33] * ydir[1]*xdir[0] + MatrixLit[42] * ydir[1]*ydir[0];
	Matrix[132] = MatrixLit[34] * xdir[0]*xdir[0] + MatrixLit[43] * xdir[0]*ydir[0] + MatrixLit[33] * ydir[0]*xdir[0] + MatrixLit[42] * ydir[0]*ydir[0];
	Matrix[133] = MatrixLit[35] * zdir[2]*xdir[0] + MatrixLit[5] * zdir[2]*ydir[0];
	Matrix[134] = MatrixLit[35] * zdir[1]*xdir[0] + MatrixLit[5] * zdir[1]*ydir[0];
	Matrix[135] = MatrixLit[35] * zdir[0]*xdir[0] + MatrixLit[5] * zdir[0]*ydir[0];
	Matrix[136] = MatrixLit[28] * xdir[1]*xdir[1] + MatrixLit[37] * xdir[1]*ydir[1] + MatrixLit[37] * ydir[1]*xdir[1] + MatrixLit[36] * ydir[1]*ydir[1];
	Matrix[137] = MatrixLit[28] * xdir[0]*xdir[1] + MatrixLit[37] * xdir[0]*ydir[1] + MatrixLit[37] * ydir[0]*xdir[1] + MatrixLit[36] * ydir[0]*ydir[1];
	Matrix[138] = MatrixLit[29] * zdir[2]*xdir[1] + MatrixLit[38] * zdir[2]*ydir[1];
	Matrix[139] = MatrixLit[29] * zdir[1]*xdir[1] + MatrixLit[38] * zdir[1]*ydir[1];
	Matrix[140] = MatrixLit[29] * zdir[0]*xdir[1] + MatrixLit[38] * zdir[0]*ydir[1];
	Matrix[141] = MatrixLit[28] * xdir[2]*xdir[1] + MatrixLit[37] * xdir[2]*ydir[1] + MatrixLit[37] * ydir[2]*xdir[1] + MatrixLit[36] * ydir[2]*ydir[1];
	Matrix[142] = MatrixLit[31] * xdir[1]*xdir[1] + MatrixLit[40] * xdir[1]*ydir[1] + MatrixLit[30] * ydir[1]*xdir[1] + MatrixLit[39] * ydir[1]*ydir[1];
	Matrix[143] = MatrixLit[31] * xdir[0]*xdir[1] + MatrixLit[40] * xdir[0]*ydir[1] + MatrixLit[30] * ydir[0]*xdir[1] + MatrixLit[39] * ydir[0]*ydir[1];
	Matrix[144] = MatrixLit[32] * zdir[2]*xdir[1] + MatrixLit[17] * zdir[2]*ydir[1];
	Matrix[145] = MatrixLit[32] * zdir[1]*xdir[1] + MatrixLit[17] * zdir[1]*ydir[1];
	Matrix[146] = MatrixLit[32] * zdir[0]*xdir[1] + MatrixLit[17] * zdir[0]*ydir[1];
	Matrix[147] = MatrixLit[31] * xdir[2]*xdir[1] + MatrixLit[40] * xdir[2]*ydir[1] + MatrixLit[30] * ydir[2]*xdir[1] + MatrixLit[39] * ydir[2]*ydir[1];
	Matrix[148] = MatrixLit[34] * xdir[1]*xdir[1] + MatrixLit[43] * xdir[1]*ydir[1] + MatrixLit[33] * ydir[1]*xdir[1] + MatrixLit[42] * ydir[1]*ydir[1];
	Matrix[149] = MatrixLit[34] * xdir[0]*xdir[1] + MatrixLit[43] * xdir[0]*ydir[1] + MatrixLit[33] * ydir[0]*xdir[1] + MatrixLit[42] * ydir[0]*ydir[1];
	Matrix[150] = MatrixLit[35] * zdir[2]*xdir[1] + MatrixLit[5] * zdir[2]*ydir[1];
	Matrix[151] = MatrixLit[35] * zdir[1]*xdir[1] + MatrixLit[5] * zdir[1]*ydir[1];
	Matrix[152] = MatrixLit[35] * zdir[0]*xdir[1] + MatrixLit[5] * zdir[0]*ydir[1];
	Matrix[153] = MatrixLit[55] * xdir[2]*xdir[2] + MatrixLit[67] * xdir[2]*ydir[2] + MatrixLit[67] * ydir[2]*xdir[2] + MatrixLit[66] * ydir[2]*ydir[2];
	Matrix[154] = MatrixLit[58] * xdir[1]*xdir[2] + MatrixLit[70] * xdir[1]*ydir[2] + MatrixLit[57] * ydir[1]*xdir[2] + MatrixLit[69] * ydir[1]*ydir[2];
	Matrix[155] = MatrixLit[58] * xdir[0]*xdir[2] + MatrixLit[70] * xdir[0]*ydir[2] + MatrixLit[57] * ydir[0]*xdir[2] + MatrixLit[69] * ydir[0]*ydir[2];
	Matrix[156] = MatrixLit[59] * zdir[2]*xdir[2] + MatrixLit[38] * zdir[2]*ydir[2];
	Matrix[157] = MatrixLit[59] * zdir[1]*xdir[2] + MatrixLit[38] * zdir[1]*ydir[2];
	Matrix[158] = MatrixLit[59] * zdir[0]*xdir[2] + MatrixLit[38] * zdir[0]*ydir[2];
	Matrix[159] = MatrixLit[58] * xdir[2]*xdir[2] + MatrixLit[70] * xdir[2]*ydir[2] + MatrixLit[57] * ydir[2]*xdir[2] + MatrixLit[69] * ydir[2]*ydir[2];
	Matrix[160] = MatrixLit[61] * xdir[1]*xdir[2] + MatrixLit[73] * xdir[1]*ydir[2] + MatrixLit[60] * ydir[1]*xdir[2] + MatrixLit[72] * ydir[1]*ydir[2];
	Matrix[161] = MatrixLit[61] * xdir[0]*xdir[2] + MatrixLit[73] * xdir[0]*ydir[2] + MatrixLit[60] * ydir[0]*xdir[2] + MatrixLit[72] * ydir[0]*ydir[2];
	Matrix[162] = MatrixLit[62] * zdir[2]*xdir[2] + MatrixLit[17] * zdir[2]*ydir[2];
	Matrix[163] = MatrixLit[62] * zdir[1]*xdir[2] + MatrixLit[17] * zdir[1]*ydir[2];
	Matrix[164] = MatrixLit[62] * zdir[0]*xdir[2] + MatrixLit[17] * zdir[0]*ydir[2];
	Matrix[165] = MatrixLit[61] * xdir[2]*xdir[2] + MatrixLit[73] * xdir[2]*ydir[2] + MatrixLit[60] * ydir[2]*xdir[2] + MatrixLit[72] * ydir[2]*ydir[2];
	Matrix[166] = MatrixLit[64] * xdir[1]*xdir[2] + MatrixLit[76] * xdir[1]*ydir[2] + MatrixLit[63] * ydir[1]*xdir[2] + MatrixLit[75] * ydir[1]*ydir[2];
	Matrix[167] = MatrixLit[64] * xdir[0]*xdir[2] + MatrixLit[76] * xdir[0]*ydir[2] + MatrixLit[63] * ydir[0]*xdir[2] + MatrixLit[75] * ydir[0]*ydir[2];
	Matrix[168] = MatrixLit[65] * zdir[2]*xdir[2] + MatrixLit[5] * zdir[2]*ydir[2];
	Matrix[169] = MatrixLit[65] * zdir[1]*xdir[2] + MatrixLit[5] * zdir[1]*ydir[2];
	Matrix[170] = MatrixLit[65] * zdir[0]*xdir[2] + MatrixLit[5] * zdir[0]*ydir[2];
	Matrix[171] = MatrixLit[45] * zdir[0]*zdir[0];
	Matrix[172] = MatrixLit[56] * xdir[2]*zdir[0] + MatrixLit[68] * ydir[2]*zdir[0];
	Matrix[173] = MatrixLit[47] * xdir[1]*zdir[0] + MatrixLit[46] * ydir[1]*zdir[0];
	Matrix[174] = MatrixLit[47] * xdir[0]*zdir[0] + MatrixLit[46] * ydir[0]*zdir[0];
	Matrix[175] = MatrixLit[48] * zdir[2]*zdir[0];
	Matrix[176] = MatrixLit[48] * zdir[1]*zdir[0];
	Matrix[177] = MatrixLit[48] * zdir[0]*zdir[0];
	Matrix[178] = MatrixLit[47] * xdir[2]*zdir[0] + MatrixLit[46] * ydir[2]*zdir[0];
	Matrix[179] = MatrixLit[50] * xdir[1]*zdir[0] + MatrixLit[49] * ydir[1]*zdir[0];
	Matrix[180] = MatrixLit[50] * xdir[0]*zdir[0] + MatrixLit[49] * ydir[0]*zdir[0];
	Matrix[181] = MatrixLit[51] * zdir[2]*zdir[0];
	Matrix[182] = MatrixLit[51] * zdir[1]*zdir[0];
	Matrix[183] = MatrixLit[51] * zdir[0]*zdir[0];
	Matrix[184] = MatrixLit[50] * xdir[2]*zdir[0] + MatrixLit[49] * ydir[2]*zdir[0];
	Matrix[185] = MatrixLit[53] * xdir[1]*zdir[0] + MatrixLit[52] * ydir[1]*zdir[0];
	Matrix[186] = MatrixLit[53] * xdir[0]*zdir[0] + MatrixLit[52] * ydir[0]*zdir[0];
	Matrix[187] = MatrixLit[54] * zdir[2]*zdir[0];
	Matrix[188] = MatrixLit[54] * zdir[1]*zdir[0];
	Matrix[189] = MatrixLit[54] * zdir[0]*zdir[0];
	Matrix[190] = MatrixLit[45] * zdir[1]*zdir[1];
	Matrix[191] = MatrixLit[45] * zdir[0]*zdir[1];
	Matrix[192] = MatrixLit[56] * xdir[2]*zdir[1] + MatrixLit[68] * ydir[2]*zdir[1];
	Matrix[193] = MatrixLit[47] * xdir[1]*zdir[1] + MatrixLit[46] * ydir[1]*zdir[1];
	Matrix[194] = MatrixLit[47] * xdir[0]*zdir[1] + MatrixLit[46] * ydir[0]*zdir[1];
	Matrix[195] = MatrixLit[48] * zdir[2]*zdir[1];
	Matrix[196] = MatrixLit[48] * zdir[1]*zdir[1];
	Matrix[197] = MatrixLit[48] * zdir[0]*zdir[1];
	Matrix[198] = MatrixLit[47] * xdir[2]*zdir[1] + MatrixLit[46] * ydir[2]*zdir[1];
	Matrix[199] = MatrixLit[50] * xdir[1]*zdir[1] + MatrixLit[49] * ydir[1]*zdir[1];
	Matrix[200] = MatrixLit[50] * xdir[0]*zdir[1] + MatrixLit[49] * ydir[0]*zdir[1];
	Matrix[201] = MatrixLit[51] * zdir[2]*zdir[1];
	Matrix[202] = MatrixLit[51] * zdir[1]*zdir[1];
	Matrix[203] = MatrixLit[51] * zdir[0]*zdir[1];
	Matrix[204] = MatrixLit[50] * xdir[2]*zdir[1] + MatrixLit[49] * ydir[2]*zdir[1];
	Matrix[205] = MatrixLit[53] * xdir[1]*zdir[1] + MatrixLit[52] * ydir[1]*zdir[1];
	Matrix[206] = MatrixLit[53] * xdir[0]*zdir[1] + MatrixLit[52] * ydir[0]*zdir[1];
	Matrix[207] = MatrixLit[54] * zdir[2]*zdir[1];
	Matrix[208] = MatrixLit[54] * zdir[1]*zdir[1];
	Matrix[209] = MatrixLit[54] * zdir[0]*zdir[1];
	Matrix[210] = MatrixLit[45] * zdir[2]*zdir[2];
	Matrix[211] = MatrixLit[45] * zdir[1]*zdir[2];
	Matrix[212] = MatrixLit[45] * zdir[0]*zdir[2];
	Matrix[213] = MatrixLit[56] * xdir[2]*zdir[2] + MatrixLit[68] * ydir[2]*zdir[2];
	Matrix[214] = MatrixLit[47] * xdir[1]*zdir[2] + MatrixLit[46] * ydir[1]*zdir[2];
	Matrix[215] = MatrixLit[47] * xdir[0]*zdir[2] + MatrixLit[46] * ydir[0]*zdir[2];
	Matrix[216] = MatrixLit[48] * zdir[2]*zdir[2];
	Matrix[217] = MatrixLit[48] * zdir[1]*zdir[2];
	Matrix[218] = MatrixLit[48] * zdir[0]*zdir[2];
	Matrix[219] = MatrixLit[47] * xdir[2]*zdir[2] + MatrixLit[46] * ydir[2]*zdir[2];
	Matrix[220] = MatrixLit[50] * xdir[1]*zdir[2] + MatrixLit[49] * ydir[1]*zdir[2];
	Matrix[221] = MatrixLit[50] * xdir[0]*zdir[2] + MatrixLit[49] * ydir[0]*zdir[2];
	Matrix[222] = MatrixLit[51] * zdir[2]*zdir[2];
	Matrix[223] = MatrixLit[51] * zdir[1]*zdir[2];
	Matrix[224] = MatrixLit[51] * zdir[0]*zdir[2];
	Matrix[225] = MatrixLit[50] * xdir[2]*zdir[2] + MatrixLit[49] * ydir[2]*zdir[2];
	Matrix[226] = MatrixLit[53] * xdir[1]*zdir[2] + MatrixLit[52] * ydir[1]*zdir[2];
	Matrix[227] = MatrixLit[53] * xdir[0]*zdir[2] + MatrixLit[52] * ydir[0]*zdir[2];
	Matrix[228] = MatrixLit[54] * zdir[2]*zdir[2];
	Matrix[229] = MatrixLit[54] * zdir[1]*zdir[2];
	Matrix[230] = MatrixLit[54] * zdir[0]*zdir[2];
	Matrix[231] = MatrixLit[55] * xdir[0]*xdir[0] + MatrixLit[67] * xdir[0]*ydir[0] + MatrixLit[67] * ydir[0]*xdir[0] + MatrixLit[66] * ydir[0]*ydir[0];
	Matrix[232] = MatrixLit[56] * zdir[2]*xdir[0] + MatrixLit[68] * zdir[2]*ydir[0];
	Matrix[233] = MatrixLit[56] * zdir[1]*xdir[0] + MatrixLit[68] * zdir[1]*ydir[0];
	Matrix[234] = MatrixLit[56] * zdir[0]*xdir[0] + MatrixLit[68] * zdir[0]*ydir[0];
	Matrix[235] = MatrixLit[55] * xdir[2]*xdir[0] + MatrixLit[67] * xdir[2]*ydir[0] + MatrixLit[67] * ydir[2]*xdir[0] + MatrixLit[66] * ydir[2]*ydir[0];
	Matrix[236] = MatrixLit[58] * xdir[1]*xdir[0] + MatrixLit[70] * xdir[1]*ydir[0] + MatrixLit[57] * ydir[1]*xdir[0] + MatrixLit[69] * ydir[1]*ydir[0];
	Matrix[237] = MatrixLit[58] * xdir[0]*xdir[0] + MatrixLit[70] * xdir[0]*ydir[0] + MatrixLit[57] * ydir[0]*xdir[0] + MatrixLit[69] * ydir[0]*ydir[0];
	Matrix[238] = MatrixLit[59] * zdir[2]*xdir[0] + MatrixLit[38] * zdir[2]*ydir[0];
	Matrix[239] = MatrixLit[59] * zdir[1]*xdir[0] + MatrixLit[38] * zdir[1]*ydir[0];
	Matrix[240] = MatrixLit[59] * zdir[0]*xdir[0] + MatrixLit[38] * zdir[0]*ydir[0];
	Matrix[241] = MatrixLit[58] * xdir[2]*xdir[0] + MatrixLit[70] * xdir[2]*ydir[0] + MatrixLit[57] * ydir[2]*xdir[0] + MatrixLit[69] * ydir[2]*ydir[0];
	Matrix[242] = MatrixLit[61] * xdir[1]*xdir[0] + MatrixLit[73] * xdir[1]*ydir[0] + MatrixLit[60] * ydir[1]*xdir[0] + MatrixLit[72] * ydir[1]*ydir[0];
	Matrix[243] = MatrixLit[61] * xdir[0]*xdir[0] + MatrixLit[73] * xdir[0]*ydir[0] + MatrixLit[60] * ydir[0]*xdir[0] + MatrixLit[72] * ydir[0]*ydir[0];
	Matrix[244] = MatrixLit[62] * zdir[2]*xdir[0] + MatrixLit[17] * zdir[2]*ydir[0];
	Matrix[245] = MatrixLit[62] * zdir[1]*xdir[0] + MatrixLit[17] * zdir[1]*ydir[0];
	Matrix[246] = MatrixLit[62] * zdir[0]*xdir[0] + MatrixLit[17] * zdir[0]*ydir[0];
	Matrix[247] = MatrixLit[61] * xdir[2]*xdir[0] + MatrixLit[73] * xdir[2]*ydir[0] + MatrixLit[60] * ydir[2]*xdir[0] + MatrixLit[72] * ydir[2]*ydir[0];
	Matrix[248] = MatrixLit[64] * xdir[1]*xdir[0] + MatrixLit[76] * xdir[1]*ydir[0] + MatrixLit[63] * ydir[1]*xdir[0] + MatrixLit[75] * ydir[1]*ydir[0];
	Matrix[249] = MatrixLit[64] * xdir[0]*xdir[0] + MatrixLit[76] * xdir[0]*ydir[0] + MatrixLit[63] * ydir[0]*xdir[0] + MatrixLit[75] * ydir[0]*ydir[0];
	Matrix[250] = MatrixLit[65] * zdir[2]*xdir[0] + MatrixLit[5] * zdir[2]*ydir[0];
	Matrix[251] = MatrixLit[65] * zdir[1]*xdir[0] + MatrixLit[5] * zdir[1]*ydir[0];
	Matrix[252] = MatrixLit[65] * zdir[0]*xdir[0] + MatrixLit[5] * zdir[0]*ydir[0];
	Matrix[253] = MatrixLit[55] * xdir[1]*xdir[1] + MatrixLit[67] * xdir[1]*ydir[1] + MatrixLit[67] * ydir[1]*xdir[1] + MatrixLit[66] * ydir[1]*ydir[1];
	Matrix[254] = MatrixLit[55] * xdir[0]*xdir[1] + MatrixLit[67] * xdir[0]*ydir[1] + MatrixLit[67] * ydir[0]*xdir[1] + MatrixLit[66] * ydir[0]*ydir[1];
	Matrix[255] = MatrixLit[56] * zdir[2]*xdir[1] + MatrixLit[68] * zdir[2]*ydir[1];
	Matrix[256] = MatrixLit[56] * zdir[1]*xdir[1] + MatrixLit[68] * zdir[1]*ydir[1];
	Matrix[257] = MatrixLit[56] * zdir[0]*xdir[1] + MatrixLit[68] * zdir[0]*ydir[1];
	Matrix[258] = MatrixLit[55] * xdir[2]*xdir[1] + MatrixLit[67] * xdir[2]*ydir[1] + MatrixLit[67] * ydir[2]*xdir[1] + MatrixLit[66] * ydir[2]*ydir[1];
	Matrix[259] = MatrixLit[58] * xdir[1]*xdir[1] + MatrixLit[70] * xdir[1]*ydir[1] + MatrixLit[57] * ydir[1]*xdir[1] + MatrixLit[69] * ydir[1]*ydir[1];
	Matrix[260] = MatrixLit[58] * xdir[0]*xdir[1] + MatrixLit[70] * xdir[0]*ydir[1] + MatrixLit[57] * ydir[0]*xdir[1] + MatrixLit[69] * ydir[0]*ydir[1];
	Matrix[261] = MatrixLit[59] * zdir[2]*xdir[1] + MatrixLit[38] * zdir[2]*ydir[1];
	Matrix[262] = MatrixLit[59] * zdir[1]*xdir[1] + MatrixLit[38] * zdir[1]*ydir[1];
	Matrix[263] = MatrixLit[59] * zdir[0]*xdir[1] + MatrixLit[38] * zdir[0]*ydir[1];
	Matrix[264] = MatrixLit[58] * xdir[2]*xdir[1] + MatrixLit[70] * xdir[2]*ydir[1] + MatrixLit[57] * ydir[2]*xdir[1] + MatrixLit[69] * ydir[2]*ydir[1];
	Matrix[265] = MatrixLit[61] * xdir[1]*xdir[1] + MatrixLit[73] * xdir[1]*ydir[1] + MatrixLit[60] * ydir[1]*xdir[1] + MatrixLit[72] * ydir[1]*ydir[1];
	Matrix[266] = MatrixLit[61] * xdir[0]*xdir[1] + MatrixLit[73] * xdir[0]*ydir[1] + MatrixLit[60] * ydir[0]*xdir[1] + MatrixLit[72] * ydir[0]*ydir[1];
	Matrix[267] = MatrixLit[62] * zdir[2]*xdir[1] + MatrixLit[17] * zdir[2]*ydir[1];
	Matrix[268] = MatrixLit[62] * zdir[1]*xdir[1] + MatrixLit[17] * zdir[1]*ydir[1];
	Matrix[269] = MatrixLit[62] * zdir[0]*xdir[1] + MatrixLit[17] * zdir[0]*ydir[1];
	Matrix[270] = MatrixLit[61] * xdir[2]*xdir[1] + MatrixLit[73] * xdir[2]*ydir[1] + MatrixLit[60] * ydir[2]*xdir[1] + MatrixLit[72] * ydir[2]*ydir[1];
	Matrix[271] = MatrixLit[64] * xdir[1]*xdir[1] + MatrixLit[76] * xdir[1]*ydir[1] + MatrixLit[63] * ydir[1]*xdir[1] + MatrixLit[75] * ydir[1]*ydir[1];
	Matrix[272] = MatrixLit[64] * xdir[0]*xdir[1] + MatrixLit[76] * xdir[0]*ydir[1] + MatrixLit[63] * ydir[0]*xdir[1] + MatrixLit[75] * ydir[0]*ydir[1];
	Matrix[273] = MatrixLit[65] * zdir[2]*xdir[1] + MatrixLit[5] * zdir[2]*ydir[1];
	Matrix[274] = MatrixLit[65] * zdir[1]*xdir[1] + MatrixLit[5] * zdir[1]*ydir[1];
	Matrix[275] = MatrixLit[65] * zdir[0]*xdir[1] + MatrixLit[5] * zdir[0]*ydir[1];
	Matrix[276] = MatrixLit[91] * xdir[2]*xdir[2] + MatrixLit[106] * xdir[2]*ydir[2] + MatrixLit[106] * ydir[2]*xdir[2] + MatrixLit[105] * ydir[2]*ydir[2];
	Matrix[277] = MatrixLit[94] * xdir[1]*xdir[2] + MatrixLit[109] * xdir[1]*ydir[2] + MatrixLit[93] * ydir[1]*xdir[2] + MatrixLit[108] * ydir[1]*ydir[2];
	Matrix[278] = MatrixLit[94] * xdir[0]*xdir[2] + MatrixLit[109] * xdir[0]*ydir[2] + MatrixLit[93] * ydir[0]*xdir[2] + MatrixLit[108] * ydir[0]*ydir[2];
	Matrix[279] = MatrixLit[95] * zdir[2]*xdir[2] + MatrixLit[68] * zdir[2]*ydir[2];
	Matrix[280] = MatrixLit[95] * zdir[1]*xdir[2] + MatrixLit[68] * zdir[1]*ydir[2];
	Matrix[281] = MatrixLit[95] * zdir[0]*xdir[2] + MatrixLit[68] * zdir[0]*ydir[2];
	Matrix[282] = MatrixLit[94] * xdir[2]*xdir[2] + MatrixLit[109] * xdir[2]*ydir[2] + MatrixLit[93] * ydir[2]*xdir[2] + MatrixLit[108] * ydir[2]*ydir[2];
	Matrix[283] = MatrixLit[97] * xdir[1]*xdir[2] + MatrixLit[112] * xdir[1]*ydir[2] + MatrixLit[96] * ydir[1]*xdir[2] + MatrixLit[111] * ydir[1]*ydir[2];
	Matrix[284] = MatrixLit[97] * xdir[0]*xdir[2] + MatrixLit[112] * xdir[0]*ydir[2] + MatrixLit[96] * ydir[0]*xdir[2] + MatrixLit[111] * ydir[0]*ydir[2];
	Matrix[285] = MatrixLit[98] * zdir[2]*xdir[2] + MatrixLit[38] * zdir[2]*ydir[2];
	Matrix[286] = MatrixLit[98] * zdir[1]*xdir[2] + MatrixLit[38] * zdir[1]*ydir[2];
	Matrix[287] = MatrixLit[98] * zdir[0]*xdir[2] + MatrixLit[38] * zdir[0]*ydir[2];
	Matrix[288] = MatrixLit[97] * xdir[2]*xdir[2] + MatrixLit[112] * xdir[2]*ydir[2] + MatrixLit[96] * ydir[2]*xdir[2] + MatrixLit[111] * ydir[2]*ydir[2];
	Matrix[289] = MatrixLit[100] * xdir[1]*xdir[2] + MatrixLit[115] * xdir[1]*ydir[2] + MatrixLit[99] * ydir[1]*xdir[2] + MatrixLit[114] * ydir[1]*ydir[2];
	Matrix[290] = MatrixLit[100] * xdir[0]*xdir[2] + MatrixLit[115] * xdir[0]*ydir[2] + MatrixLit[99] * ydir[0]*xdir[2] + MatrixLit[114] * ydir[0]*ydir[2];
	Matrix[291] = MatrixLit[101] * zdir[2]*xdir[2] + MatrixLit[17] * zdir[2]*ydir[2];
	Matrix[292] = MatrixLit[101] * zdir[1]*xdir[2] + MatrixLit[17] * zdir[1]*ydir[2];
	Matrix[293] = MatrixLit[101] * zdir[0]*xdir[2] + MatrixLit[17] * zdir[0]*ydir[2];
	Matrix[294] = MatrixLit[100] * xdir[2]*xdir[2] + MatrixLit[115] * xdir[2]*ydir[2] + MatrixLit[99] * ydir[2]*xdir[2] + MatrixLit[114] * ydir[2]*ydir[2];
	Matrix[295] = MatrixLit[103] * xdir[1]*xdir[2] + MatrixLit[118] * xdir[1]*ydir[2] + MatrixLit[102] * ydir[1]*xdir[2] + MatrixLit[117] * ydir[1]*ydir[2];
	Matrix[296] = MatrixLit[103] * xdir[0]*xdir[2] + MatrixLit[118] * xdir[0]*ydir[2] + MatrixLit[102] * ydir[0]*xdir[2] + MatrixLit[117] * ydir[0]*ydir[2];
	Matrix[297] = MatrixLit[104] * zdir[2]*xdir[2] + MatrixLit[5] * zdir[2]*ydir[2];
	Matrix[298] = MatrixLit[104] * zdir[1]*xdir[2] + MatrixLit[5] * zdir[1]*ydir[2];
	Matrix[299] = MatrixLit[104] * zdir[0]*xdir[2] + MatrixLit[5] * zdir[0]*ydir[2];
/*#ifdef _DEBUG_
	for (unsigned int i=0;i<78;++i)
		{
			cout<<"Matrix element "<<i<<" is "<<Matrix[i]<<endl;
	};
#endif*/


}

//	Calculate element stress 
//Stress regulation is put aside for now
void CPlate::ElementStress(double* stress, double* Displacement,double* position)
{
	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial);	// Pointer to material of the element

	double xdir[3];
	double ydir[3];
	double zdir[3];
	xdir[0]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0]);
	xdir[1]=(nodes[1]->XYZ[1]-nodes[0]->XYZ[1]);
	xdir[2]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2]);
	ydir[0]=(nodes[3]->XYZ[0]-nodes[0]->XYZ[0]);
	ydir[1]=(nodes[3]->XYZ[1]-nodes[0]->XYZ[1]);
	ydir[2]=(nodes[3]->XYZ[2]-nodes[0]->XYZ[2]);
	double LX2=xdir[0]*xdir[0]+xdir[1]*xdir[1]+xdir[2]*xdir[2];
	double LX=sqrt(LX2);
	double LY2=ydir[0]*ydir[0]+ydir[1]*ydir[1]+ydir[2]*ydir[2];
	double LY=sqrt(LY2);
	xdir[0]=xdir[0]/LX;
	xdir[1]=xdir[1]/LX;
	xdir[2]=xdir[2]/LX;
	ydir[0]=ydir[0]/LY;
	ydir[1]=ydir[1]/LY;
	ydir[2]=ydir[2]/LY;
	zdir[0]=xdir[1]*ydir[2]-xdir[2]*ydir[1];
	zdir[1]=xdir[2]*ydir[0]-xdir[0]*ydir[2];
	zdir[2]=xdir[0]*ydir[1]-xdir[1]*ydir[0];

	double xpsi(LX);
	double yeta(LY);

	double dis[24];//displacement of nodes
	for (unsigned i=0;i<=11;i++){
		if (LocationMatrix[i])
		dis[i+1]=Displacement[LocationMatrix[i]-1];
		else
			dis[i]=0.0;

	};
	

	double nu=material->nu;
	double Jacobian=xpsi*yeta;
	const double k= material->E*material->h*material->h*material->h*material->h*0.5*0.0833333333333333/(1-nu*nu);

	double psix=yeta/Jacobian;
	double etay=xpsi/Jacobian;

	/*well it's just too complex a set of formulas...*/
stress[0] = dis[0]*(psix*psix*  0.6830127018922193+nu*etay*etay*  0.6830127018922193)+dis[1]*(psix*psix*  0.0000000000000000+nu*etay*etay*  1.0773502691896257)+dis[2]*(psix*psix* -1.0773502691896257+nu*etay*etay*  0.0000000000000000)+dis[3]*(psix*psix* -0.6830127018922193+nu*etay*etay*  0.1830127018922193)+dis[4]*(psix*psix*  0.0000000000000000+nu*etay*etay*  0.2886751345948129)+dis[5]*(psix*psix* -0.2886751345948129+nu*etay*etay*  0.0000000000000000)+dis[6]*(psix*psix* -0.1830127018922193+nu*etay*etay* -0.1830127018922193)+dis[7]*(psix*psix*  0.0000000000000000+nu*etay*etay*  0.0773502691896258)+dis[8]*(psix*psix* -0.0773502691896258+nu*etay*etay*  0.0000000000000000)+dis[9]*(psix*psix*  0.1830127018922193+nu*etay*etay* -0.6830127018922193)+dis[10]*(psix*psix*  0.0000000000000000+nu*etay*etay*  0.2886751345948129)+dis[11]*(psix*psix* -0.2886751345948129+nu*etay*etay*  0.0000000000000000);
stress[1] = dis[0]*(nu*psix*psix*  0.6830127018922193+etay*etay*  0.6830127018922193)+dis[1]*(nu*psix*psix*  0.0000000000000000+etay*etay*  1.0773502691896257)+dis[2]*(nu*psix*psix* -1.0773502691896257+etay*etay*  0.0000000000000000)+dis[3]*(nu*psix*psix* -0.6830127018922193+etay*etay*  0.1830127018922193)+dis[4]*(nu*psix*psix*  0.0000000000000000+etay*etay*  0.2886751345948129)+dis[5]*(nu*psix*psix* -0.2886751345948129+etay*etay*  0.0000000000000000)+dis[6]*(nu*psix*psix* -0.1830127018922193+etay*etay* -0.1830127018922193)+dis[7]*(nu*psix*psix*  0.0000000000000000+etay*etay*  0.0773502691896258)+dis[8]*(nu*psix*psix* -0.0773502691896258+etay*etay*  0.0000000000000000)+dis[9]*(nu*psix*psix*  0.1830127018922193+etay*etay* -0.6830127018922193)+dis[10]*(nu*psix*psix*  0.0000000000000000+etay*etay*  0.2886751345948129)+dis[11]*(nu*psix*psix* -0.2886751345948129+etay*etay*  0.0000000000000000);
stress[2] = dis[0]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[2]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[3]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[5]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[6]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[8]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[9]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[11]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129);
stress[3] = dis[0]*(psix*psix* -0.6830127018922193+nu*etay*etay*  0.1830127018922193)+dis[1]*(psix*psix*  0.0000000000000000+nu*etay*etay*  0.2886751345948129)+dis[2]*(psix*psix*  0.2886751345948129+nu*etay*etay*  0.0000000000000000)+dis[3]*(psix*psix*  0.6830127018922193+nu*etay*etay*  0.6830127018922193)+dis[4]*(psix*psix*  0.0000000000000000+nu*etay*etay*  1.0773502691896257)+dis[5]*(psix*psix*  1.0773502691896257+nu*etay*etay*  0.0000000000000000)+dis[6]*(psix*psix*  0.1830127018922193+nu*etay*etay* -0.6830127018922193)+dis[7]*(psix*psix*  0.0000000000000000+nu*etay*etay*  0.2886751345948129)+dis[8]*(psix*psix*  0.2886751345948129+nu*etay*etay*  0.0000000000000000)+dis[9]*(psix*psix* -0.1830127018922193+nu*etay*etay* -0.1830127018922193)+dis[10]*(psix*psix*  0.0000000000000000+nu*etay*etay*  0.0773502691896258)+dis[11]*(psix*psix*  0.0773502691896258+nu*etay*etay*  0.0000000000000000);
stress[4] = dis[0]*(nu*psix*psix* -0.6830127018922193+etay*etay*  0.1830127018922193)+dis[1]*(nu*psix*psix*  0.0000000000000000+etay*etay*  0.2886751345948129)+dis[2]*(nu*psix*psix*  0.2886751345948129+etay*etay*  0.0000000000000000)+dis[3]*(nu*psix*psix*  0.6830127018922193+etay*etay*  0.6830127018922193)+dis[4]*(nu*psix*psix*  0.0000000000000000+etay*etay*  1.0773502691896257)+dis[5]*(nu*psix*psix*  1.0773502691896257+etay*etay*  0.0000000000000000)+dis[6]*(nu*psix*psix*  0.1830127018922193+etay*etay* -0.6830127018922193)+dis[7]*(nu*psix*psix*  0.0000000000000000+etay*etay*  0.2886751345948129)+dis[8]*(nu*psix*psix*  0.2886751345948129+etay*etay*  0.0000000000000000)+dis[9]*(nu*psix*psix* -0.1830127018922193+etay*etay* -0.1830127018922193)+dis[10]*(nu*psix*psix*  0.0000000000000000+etay*etay*  0.0773502691896258)+dis[11]*(nu*psix*psix*  0.0773502691896258+etay*etay*  0.0000000000000000);
stress[5] = dis[0]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[2]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[3]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[5]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[6]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[8]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[9]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[11]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129);
stress[6] = dis[0]*(psix*psix* -0.1830127018922193+nu*etay*etay* -0.1830127018922193)+dis[1]*(psix*psix*  0.0000000000000000+nu*etay*etay* -0.0773502691896258)+dis[2]*(psix*psix*  0.0773502691896258+nu*etay*etay*  0.0000000000000000)+dis[3]*(psix*psix*  0.1830127018922193+nu*etay*etay* -0.6830127018922193)+dis[4]*(psix*psix*  0.0000000000000000+nu*etay*etay* -0.2886751345948129)+dis[5]*(psix*psix*  0.2886751345948129+nu*etay*etay*  0.0000000000000000)+dis[6]*(psix*psix*  0.6830127018922193+nu*etay*etay*  0.6830127018922193)+dis[7]*(psix*psix*  0.0000000000000000+nu*etay*etay* -1.0773502691896257)+dis[8]*(psix*psix*  1.0773502691896257+nu*etay*etay*  0.0000000000000000)+dis[9]*(psix*psix* -0.6830127018922193+nu*etay*etay*  0.1830127018922193)+dis[10]*(psix*psix*  0.0000000000000000+nu*etay*etay* -0.2886751345948129)+dis[11]*(psix*psix*  0.2886751345948129+nu*etay*etay*  0.0000000000000000);
stress[7] = dis[0]*(nu*psix*psix* -0.1830127018922193+etay*etay* -0.1830127018922193)+dis[1]*(nu*psix*psix*  0.0000000000000000+etay*etay* -0.0773502691896258)+dis[2]*(nu*psix*psix*  0.0773502691896258+etay*etay*  0.0000000000000000)+dis[3]*(nu*psix*psix*  0.1830127018922193+etay*etay* -0.6830127018922193)+dis[4]*(nu*psix*psix*  0.0000000000000000+etay*etay* -0.2886751345948129)+dis[5]*(nu*psix*psix*  0.2886751345948129+etay*etay*  0.0000000000000000)+dis[6]*(nu*psix*psix*  0.6830127018922193+etay*etay*  0.6830127018922193)+dis[7]*(nu*psix*psix*  0.0000000000000000+etay*etay* -1.0773502691896257)+dis[8]*(nu*psix*psix*  1.0773502691896257+etay*etay*  0.0000000000000000)+dis[9]*(nu*psix*psix* -0.6830127018922193+etay*etay*  0.1830127018922193)+dis[10]*(nu*psix*psix*  0.0000000000000000+etay*etay* -0.2886751345948129)+dis[11]*(nu*psix*psix*  0.2886751345948129+etay*etay*  0.0000000000000000);
stress[8] = dis[0]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[2]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[3]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[5]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[6]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[8]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[9]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[11]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129);
stress[9] = dis[0]*(psix*psix*  0.1830127018922193+nu*etay*etay* -0.6830127018922193)+dis[1]*(psix*psix*  0.0000000000000000+nu*etay*etay* -0.2886751345948129)+dis[2]*(psix*psix* -0.2886751345948129+nu*etay*etay*  0.0000000000000000)+dis[3]*(psix*psix* -0.1830127018922193+nu*etay*etay* -0.1830127018922193)+dis[4]*(psix*psix*  0.0000000000000000+nu*etay*etay* -0.0773502691896258)+dis[5]*(psix*psix* -0.0773502691896258+nu*etay*etay*  0.0000000000000000)+dis[6]*(psix*psix* -0.6830127018922193+nu*etay*etay*  0.1830127018922193)+dis[7]*(psix*psix*  0.0000000000000000+nu*etay*etay* -0.2886751345948129)+dis[8]*(psix*psix* -0.2886751345948129+nu*etay*etay*  0.0000000000000000)+dis[9]*(psix*psix*  0.6830127018922193+nu*etay*etay*  0.6830127018922193)+dis[10]*(psix*psix*  0.0000000000000000+nu*etay*etay* -1.0773502691896257)+dis[11]*(psix*psix* -1.0773502691896257+nu*etay*etay*  0.0000000000000000);
stress[10] = dis[0]*(nu*psix*psix*  0.1830127018922193+etay*etay* -0.6830127018922193)+dis[1]*(nu*psix*psix*  0.0000000000000000+etay*etay* -0.2886751345948129)+dis[2]*(nu*psix*psix* -0.2886751345948129+etay*etay*  0.0000000000000000)+dis[3]*(nu*psix*psix* -0.1830127018922193+etay*etay* -0.1830127018922193)+dis[4]*(nu*psix*psix*  0.0000000000000000+etay*etay* -0.0773502691896258)+dis[5]*(nu*psix*psix* -0.0773502691896258+etay*etay*  0.0000000000000000)+dis[6]*(nu*psix*psix* -0.6830127018922193+etay*etay*  0.1830127018922193)+dis[7]*(nu*psix*psix*  0.0000000000000000+etay*etay* -0.2886751345948129)+dis[8]*(nu*psix*psix* -0.2886751345948129+etay*etay*  0.0000000000000000)+dis[9]*(nu*psix*psix*  0.6830127018922193+etay*etay*  0.6830127018922193)+dis[10]*(nu*psix*psix*  0.0000000000000000+etay*etay* -1.0773502691896257)+dis[11]*(nu*psix*psix* -1.0773502691896257+etay*etay*  0.0000000000000000);
stress[11] = dis[0]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[2]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[3]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[5]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[6]*(+psix*etay*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[8]*(+psix*etay*(1-nu)*0.5* -0.2886751345948129)+dis[9]*(+psix*etay*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129)+dis[11]*(+psix*etay*(1-nu)*0.5*  0.2886751345948129);
for (unsigned int i=0;i<11;++i){
	stress[i]=stress[i]*k;
}

//Here a transit to 3d is needed.
	double xmid=(nodes[0]->XYZ[0]+nodes[1]->XYZ[0]+nodes[2]->XYZ[0]+nodes[3]->XYZ[0])*0.25;
	double ymid=(nodes[0]->XYZ[1]+nodes[1]->XYZ[1]+nodes[2]->XYZ[1]+nodes[3]->XYZ[1])*0.25;
	double zmid=(nodes[0]->XYZ[2]+nodes[1]->XYZ[2]+nodes[2]->XYZ[2]+nodes[3]->XYZ[2])*0.25;
	position[0]=xmid+sqrt(1.0/3)*(-xpsi);
	position[1]=ymid+sqrt(1.0/3)*(-yeta);
	position[2]=zmid;
	position[3]=xmid+sqrt(1.0/3)*(+xpsi);
	position[4]=ymid+sqrt(1.0/3)*(-yeta);
	position[5]=zmid;
	position[6]=xmid+sqrt(1.0/3)*(+xpsi);
	position[7]=ymid+sqrt(1.0/3)*(yeta);
	position[8]=zmid;
	position[9]=xmid+sqrt(1.0/3)*(-xpsi);
	position[10]=ymid+sqrt(1.0/3)*(yeta);
	position[11]=zmid;

}
