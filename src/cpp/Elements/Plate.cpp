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
	
	ND = 12;
    LocationMatrix = new unsigned int[12];

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
	if (abs(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])>FLT_EPSILON||abs(nodes[2]->XYZ[2]-nodes[0]->XYZ[2])>FLT_EPSILON||abs(nodes[3]->XYZ[2]-nodes[0]->XYZ[2])>FLT_EPSILON){
		cout<<"A plate out of XY plane is not supported currently"<<endl;
#ifdef _DEBUG_
		cout<<setw(20)<<nodes[0]->XYZ[2]<<setw(20)<<nodes[1]->XYZ[2]<<setw(20)<<nodes[2]->XYZ[2]<<setw(20)<<nodes[3]->XYZ[2]<<endl;
#endif
		return false;
	}
	if (abs(nodes[3]->XYZ[0]-nodes[0]->XYZ[0]+nodes[1]->XYZ[0]-nodes[2]->XYZ[0])>FLT_EPSILON||abs(nodes[3]->XYZ[1]-nodes[0]->XYZ[1]+nodes[1]->XYZ[1]-nodes[2]->XYZ[1])>FLT_EPSILON){
		cout<<"A shape of parallelogram is needed to guarantee convergence"<<endl;
		return false;
	}
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
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D+2];
	//This plate is only a test one. The ultimate plate would be released in shell.
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CPlate::SizeOfStiffnessMatrix() { return 78; }

/*it makes me think of another question that when should the rotation dof turn to zero*/ 


//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CPlate::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	
	double xpsi=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])*0.5;
	double xeta=(nodes[1]->XYZ[1]-nodes[0]->XYZ[1])*0.5;
	double ypsi=(nodes[3]->XYZ[0]-nodes[0]->XYZ[0])*0.5;
	double yeta=(nodes[3]->XYZ[1]-nodes[0]->XYZ[1])*0.5;

	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial);	// Pointer to material of the element

	double nu=material->nu;
	double Jacobian=xpsi*yeta-ypsi*xeta;
	const double k= material->E*material->h*material->h*material->h*0.0833333333333333/(1-nu*nu)*Jacobian;

	double psix=yeta/Jacobian;
	double psiy=-xeta/Jacobian;
	double etax=-ypsi/Jacobian;
	double etay=xpsi/Jacobian;

#ifdef _DEBUG_
	cout<<"Jacobian"<<setw(20)<<psix<<setw(20)<<psiy<<setw(20)<<etax<<setw(5)<<etay;
#endif
	Matrix[0] = k*(etax*etax*etax*etax + etay*etay*etay*etay + psix*psix*psix*psix + psiy*psiy*psiy*psiy + 2*etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[1] = k*(etax*etax*etax*psix + etay*etay*etay*psiy + (4*etax*etax*etax*etax)/3 + (4*etay*etay*etay*etay)/3 + (8*etax*etax*etay*etay)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + etax*etay*etay*psix + etax*etax*etay*psiy + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[2] = k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[3] = k*(etax*psix*psix*psix + etay*psiy*psiy*psiy + (4*psix*psix*psix*psix)/3 + (4*psiy*psiy*psiy*psiy)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 + (8*psix*psix*psiy*psiy)/3 + etax*psix*psiy*psiy + etay*psix*psix*psiy - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[4] = -k*((etax*psix*psix*psix)/6 + (etax*etax*etax*psix)/6 + (etay*psiy*psiy*psiy)/6 + (etay*etay*etay*psiy)/6 + etax*etax*psix*psix + etay*etay*psiy*psiy + (etax*psix*psiy*psiy)/6 + (etay*psix*psix*psiy)/6 + etax*etax*nu*psiy*psiy + etay*etay*nu*psix*psix + (etax*etay*etay*psix)/6 + (etax*etax*etay*psiy)/6 + 2*etax*etay*psix*psiy - 2*etax*etay*nu*psix*psiy);
	Matrix[5] = -k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[6] = k*(etax*etax*etax*etax + etay*etay*etay*etay + psix*psix*psix*psix + psiy*psiy*psiy*psiy + 2*etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[7] = k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[8] = k*(etax*etax*etax*psix + etay*etay*etay*psiy + etax*etax*etax*etax/2 + etay*etay*etay*etay/2 + etax*etax*etay*etay - (7*etax*etax*psix*psix)/10 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (7*etay*etay*psiy*psiy)/10 - (2*etax*etax*nu*psiy*psiy)/5 - (2*etay*etay*nu*psix*psix)/5 + etax*etay*etay*psix + etax*etax*etay*psiy - (6*etax*etay*psix*psiy)/5 + (4*etax*etay*nu*psix*psiy)/5);
	Matrix[9] = -k*(psix*psix*psix*psix - etay*etay*etay*etay/2 - etax*etax*etax*etax/2 + psiy*psiy*psiy*psiy - etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[10] = k*((4*etax*etax*etax*etax)/3 - etay*etay*etay*psiy - etax*etax*etax*psix + (4*etay*etay*etay*etay)/3 + (8*etax*etax*etay*etay)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 - etax*etay*etay*psix - etax*etax*etay*psiy + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[11] = k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[12] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[13] = -k*((8*etax*etax*psix*psix)/15 - (2*etay*etay*etay*etay)/3 - (4*etax*etax*etay*etay)/3 - (2*etax*etax*etax*etax)/3 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[14] = -k*(etax*etax*etax*psix + etay*etay*etay*psiy - etax*etax*etax*etax/2 - etay*etay*etay*etay/2 - etax*etax*etay*etay + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + etax*etay*etay*psix + etax*etax*etay*psiy + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[15] = k*((4*psix*psix*psix*psix)/3 - etay*psiy*psiy*psiy - etax*psix*psix*psix + (4*psiy*psiy*psiy*psiy)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 + (8*psix*psix*psiy*psiy)/3 - etax*psix*psiy*psiy - etay*psix*psix*psiy - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[16] = -k*((etax*psix*psix*psix)/6 + (etax*etax*etax*psix)/6 + (etay*psiy*psiy*psiy)/6 + (etay*etay*etay*psiy)/6 - etax*etax*psix*psix - etay*etay*psiy*psiy + (etax*psix*psiy*psiy)/6 + (etay*psix*psix*psiy)/6 - etax*etax*nu*psiy*psiy - etay*etay*nu*psix*psix + (etax*etay*etay*psix)/6 + (etax*etax*etay*psiy)/6 - 2*etax*etay*psix*psiy + 2*etax*etay*nu*psix*psiy);
	Matrix[17] = k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[18] = -k*((2*etax*etax*psix*psix)/15 - (2*psiy*psiy*psiy*psiy)/3 - (2*psix*psix*psix*psix)/3 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 - (4*psix*psix*psiy*psiy)/3 - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[19] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[20] = -k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[21] = k*(etax*etax*etax*etax + etay*etay*etay*etay + psix*psix*psix*psix + psiy*psiy*psiy*psiy + 2*etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[22] = -k*(etax*psix*psix*psix + etay*psiy*psiy*psiy - psix*psix*psix*psix/2 - psiy*psiy*psiy*psiy/2 + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 - psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[23] = -k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[24] = -k*(etax*etax*etax*etax + etay*etay*etay*etay - psix*psix*psix*psix/2 - psiy*psiy*psiy*psiy/2 + 2*etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 - psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[25] = k*(etax*psix*psix*psix + etay*psiy*psiy*psiy + psix*psix*psix*psix/2 + psiy*psiy*psiy*psiy/2 - (etax*etax*psix*psix)/5 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (etay*etay*psiy*psiy)/5 + psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy + (etax*etax*nu*psiy*psiy)/10 + (etay*etay*nu*psix*psix)/10 - (etax*etay*psix*psiy)/5 - (etax*etay*nu*psix*psiy)/5);
	Matrix[26] = -k*(etax*etax*etax*psix + etay*etay*etay*psiy + etax*etax*etax*etax/2 + etay*etay*etay*etay/2 + etax*etax*etay*etay - (etax*etax*psix*psix)/5 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (etay*etay*psiy*psiy)/5 + (etax*etax*nu*psiy*psiy)/10 + (etay*etay*nu*psix*psix)/10 + etax*etay*etay*psix + etax*etax*etay*psiy - (etax*etay*psix*psiy)/5 - (etax*etay*nu*psix*psiy)/5);
	Matrix[27] = -k*(etax*etax*etax*etax/2 + etay*etay*etay*etay/2 + psix*psix*psix*psix/2 + psiy*psiy*psiy*psiy/2 + etax*etax*etay*etay - (19*etax*etax*psix*psix)/10 - (7*etax*etax*psiy*psiy)/10 - (7*etay*etay*psix*psix)/10 - (19*etay*etay*psiy*psiy)/10 + psix*psix*psiy*psiy + (etax*etax*nu*psiy*psiy)/5 + (etay*etay*nu*psix*psix)/5 - (12*etax*etay*psix*psiy)/5 - (2*etax*etay*nu*psix*psiy)/5);
	Matrix[28] = k*(etax*etax*etax*psix + etay*etay*etay*psiy + (4*etax*etax*etax*etax)/3 + (4*etay*etay*etay*etay)/3 + (8*etax*etax*etay*etay)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + etax*etay*etay*psix + etax*etax*etay*psiy + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[29] = -k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[30] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[31] = -k*((2*etax*etax*psix*psix)/15 - (2*etay*etay*etay*etay)/3 - (4*etax*etax*etay*etay)/3 - (2*etax*etax*etax*etax)/3 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[32] = k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[33] = -(k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[34] = k*(etax*etax*etax*psix + etay*etay*etay*psiy + etax*etax*etax*etax/3 + etay*etay*etay*etay/3 + (2*etax*etax*etay*etay)/3 + (2*etax*etax*psix*psix)/15 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 + etax*etay*etay*psix + etax*etax*etay*psiy + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[35] = k*(etax*etax*etax*psix + etay*etay*etay*psiy + etax*etax*etax*etax/2 + etay*etay*etay*etay/2 + etax*etax*etay*etay - (etax*etax*psix*psix)/5 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (etay*etay*psiy*psiy)/5 + (etax*etax*nu*psiy*psiy)/10 + (etay*etay*nu*psix*psix)/10 + etax*etay*etay*psix + etax*etax*etay*psiy - (etax*etay*psix*psiy)/5 - (etax*etay*nu*psix*psiy)/5);
	Matrix[36] = k*(etax*psix*psix*psix + etay*psiy*psiy*psiy + (4*psix*psix*psix*psix)/3 + (4*psiy*psiy*psiy*psiy)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 + (8*psix*psix*psiy*psiy)/3 + etax*psix*psiy*psiy + etay*psix*psix*psiy - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[37] = -k*((etax*psix*psix*psix)/6 + (etax*etax*etax*psix)/6 + (etay*psiy*psiy*psiy)/6 + (etay*etay*etay*psiy)/6 + etax*etax*psix*psix + etay*etay*psiy*psiy + (etax*psix*psiy*psiy)/6 + (etay*psix*psix*psiy)/6 + etax*etax*nu*psiy*psiy + etay*etay*nu*psix*psix + (etax*etay*etay*psix)/6 + (etax*etax*etay*psiy)/6 + 2*etax*etay*psix*psiy - 2*etax*etay*nu*psix*psiy);
	Matrix[38] = k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[39] = -k*((8*etax*etax*psix*psix)/15 - (2*psiy*psiy*psiy*psiy)/3 - (2*psix*psix*psix*psix)/3 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*psix*psix*psiy*psiy)/3 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[40] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[41] = k*(etax*psix*psix*psix + etay*psiy*psiy*psiy + psix*psix*psix*psix/2 + psiy*psiy*psiy*psiy/2 - (7*etax*etax*psix*psix)/10 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (7*etay*etay*psiy*psiy)/10 + psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy - (2*etax*etax*nu*psiy*psiy)/5 - (2*etay*etay*nu*psix*psix)/5 - (6*etax*etay*psix*psiy)/5 + (4*etax*etay*nu*psix*psiy)/5);
	Matrix[42] = k*(etax*psix*psix*psix + etay*psiy*psiy*psiy + psix*psix*psix*psix/3 + psiy*psiy*psiy*psiy/3 + (2*etax*etax*psix*psix)/15 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 + (2*psix*psix*psiy*psiy)/3 + etax*psix*psiy*psiy + etay*psix*psix*psiy - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[43] = -(k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[44] = -k*(etax*psix*psix*psix + etay*psiy*psiy*psiy + psix*psix*psix*psix/2 + psiy*psiy*psiy*psiy/2 - (etax*etax*psix*psix)/5 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (etay*etay*psiy*psiy)/5 + psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy + (etax*etax*nu*psiy*psiy)/10 + (etay*etay*nu*psix*psix)/10 - (etax*etay*psix*psiy)/5 - (etax*etay*nu*psix*psiy)/5);
	Matrix[45] = k*(etax*etax*etax*etax + etay*etay*etay*etay + psix*psix*psix*psix + psiy*psiy*psiy*psiy + 2*etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[46] = -k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[47] = -k*(etax*etax*etax*psix + etay*etay*etay*psiy + etax*etax*etax*etax/2 + etay*etay*etay*etay/2 + etax*etax*etay*etay - (7*etax*etax*psix*psix)/10 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (7*etay*etay*psiy*psiy)/10 - (2*etax*etax*nu*psiy*psiy)/5 - (2*etay*etay*nu*psix*psix)/5 + etax*etay*etay*psix + etax*etax*etay*psiy - (6*etax*etay*psix*psiy)/5 + (4*etax*etay*nu*psix*psiy)/5);
	Matrix[48] = -k*(psix*psix*psix*psix - etay*etay*etay*etay/2 - etax*etax*etax*etax/2 + psiy*psiy*psiy*psiy - etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[49] = k*(etax*psix*psix*psix + etay*psiy*psiy*psiy - psix*psix*psix*psix/2 - psiy*psiy*psiy*psiy/2 + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[50] = k*(etax*etax*etax*psix + etay*etay*etay*psiy - etax*etax*etax*etax/2 - etay*etay*etay*etay/2 - etax*etax*etay*etay + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + etax*etay*etay*psix + etax*etax*etay*psiy + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[51] = -k*(etax*etax*etax*etax/2 + etay*etay*etay*etay/2 + psix*psix*psix*psix/2 + psiy*psiy*psiy*psiy/2 + etax*etax*etay*etay - (19*etax*etax*psix*psix)/10 - (7*etax*etax*psiy*psiy)/10 - (7*etay*etay*psix*psix)/10 - (19*etay*etay*psiy*psiy)/10 + psix*psix*psiy*psiy + (etax*etax*nu*psiy*psiy)/5 + (etay*etay*nu*psix*psix)/5 - (12*etax*etay*psix*psiy)/5 - (2*etax*etay*nu*psix*psiy)/5);
	Matrix[52] = -k*(etax*psix*psix*psix + etay*psiy*psiy*psiy + psix*psix*psix*psix/2 + psiy*psiy*psiy*psiy/2 - (7*etax*etax*psix*psix)/10 - (etax*etax*psiy*psiy)/10 - (etay*etay*psix*psix)/10 - (7*etay*etay*psiy*psiy)/10 + psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy - (2*etax*etax*nu*psiy*psiy)/5 - (2*etay*etay*nu*psix*psix)/5 - (6*etax*etay*psix*psiy)/5 + (4*etax*etay*nu*psix*psiy)/5);
	Matrix[53] = -k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[54] = -k*(etax*etax*etax*etax + etay*etay*etay*etay - psix*psix*psix*psix/2 - psiy*psiy*psiy*psiy/2 + 2*etax*etax*etay*etay + (19*etax*etax*psix*psix)/10 + (7*etax*etax*psiy*psiy)/10 + (7*etay*etay*psix*psix)/10 + (19*etay*etay*psiy*psiy)/10 - psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/5 - (etay*etay*nu*psix*psix)/5 + (12*etax*etay*psix*psiy)/5 + (2*etax*etay*nu*psix*psiy)/5);
	Matrix[55] = k*((4*etax*etax*etax*etax)/3 - etay*etay*etay*psiy - etax*etax*etax*psix + (4*etay*etay*etay*etay)/3 + (8*etax*etax*etay*etay)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 - etax*etay*etay*psix - etax*etax*etay*psiy + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[56] = -k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[57] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[58] = -k*((8*etax*etax*psix*psix)/15 - (2*etay*etay*etay*etay)/3 - (4*etax*etax*etay*etay)/3 - (2*etax*etax*etax*etax)/3 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[59] = k*(etax*etax*etax*psix + etay*etay*etay*psiy - etax*etax*etax*etax/2 - etay*etay*etay*etay/2 - etax*etax*etay*etay + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + etax*etay*etay*psix + etax*etax*etay*psiy + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[60] = -(k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[61] = k*(etax*etax*etax*etax/3 - etay*etay*etay*psiy - etax*etax*etax*psix + etay*etay*etay*etay/3 + (2*etax*etax*etay*etay)/3 + (2*etax*etax*psix*psix)/15 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 - etax*etay*etay*psix - etax*etax*etay*psiy + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[62] = -k*(etax*etax*etax*psix + etay*etay*etay*psiy - etax*etax*etax*etax/2 - etay*etay*etay*etay/2 - etax*etax*etay*etay + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + etax*etay*etay*psix + etax*etax*etay*psiy + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[63] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[64] = -k*((2*etax*etax*psix*psix)/15 - (2*etay*etay*etay*etay)/3 - (4*etax*etax*etay*etay)/3 - (2*etax*etax*etax*etax)/3 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[65] = k*(etax*etax*etax*etax + etay*etay*etay*etay + 2*etax*etax*etay*etay + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[66] = k*((4*psix*psix*psix*psix)/3 - etay*psiy*psiy*psiy - etax*psix*psix*psix + (4*psiy*psiy*psiy*psiy)/3 + (8*etax*etax*psix*psix)/15 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 + (8*psix*psix*psiy*psiy)/3 - etax*psix*psiy*psiy - etay*psix*psix*psiy - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[67] = -k*((etax*psix*psix*psix)/6 + (etax*etax*etax*psix)/6 + (etay*psiy*psiy*psiy)/6 + (etay*etay*etay*psiy)/6 - etax*etax*psix*psix - etay*etay*psiy*psiy + (etax*psix*psiy*psiy)/6 + (etay*psix*psix*psiy)/6 - etax*etax*nu*psiy*psiy - etay*etay*nu*psix*psix + (etax*etay*etay*psix)/6 + (etax*etax*etay*psiy)/6 - 2*etax*etay*psix*psiy + 2*etax*etay*nu*psix*psiy);
	Matrix[68] = -k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 + 2*psix*psix*psiy*psiy + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);
	Matrix[69] = -k*((2*etax*etax*psix*psix)/15 - (2*psiy*psiy*psiy*psiy)/3 - (2*psix*psix*psix*psix)/3 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 - (4*psix*psix*psiy*psiy)/3 - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[70] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[71] = k*(psix*psix*psix*psix + psiy*psiy*psiy*psiy + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 + 2*psix*psix*psiy*psiy - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[72] = k*(psix*psix*psix*psix/3 - etay*psiy*psiy*psiy - etax*psix*psix*psix + psiy*psiy*psiy*psiy/3 + (2*etax*etax*psix*psix)/15 + (etax*etax*psiy*psiy)/15 + (etay*etay*psix*psix)/15 + (2*etay*etay*psiy*psiy)/15 + (2*psix*psix*psiy*psiy)/3 - etax*psix*psiy*psiy - etay*psix*psix*psiy - (etax*etax*nu*psiy*psiy)/15 - (etay*etay*nu*psix*psix)/15 + (2*etax*etay*psix*psiy)/15 + (2*etax*etay*nu*psix*psiy)/15);
	Matrix[73] = -(k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[74] = -k*(etax*psix*psix*psix + etay*psiy*psiy*psiy - psix*psix*psix*psix/2 - psiy*psiy*psiy*psiy/2 + (etax*etax*psix*psix)/5 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (etay*etay*psiy*psiy)/5 - psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy - (etax*etax*nu*psiy*psiy)/10 - (etay*etay*nu*psix*psix)/10 + (etax*etay*psix*psiy)/5 + (etax*etay*nu*psix*psiy)/5);
	Matrix[75] = -k*((8*etax*etax*psix*psix)/15 - (2*psiy*psiy*psiy*psiy)/3 - (2*psix*psix*psix*psix)/3 + (4*etax*etax*psiy*psiy)/15 + (4*etay*etay*psix*psix)/15 + (8*etay*etay*psiy*psiy)/15 - (4*psix*psix*psiy*psiy)/3 - (4*etax*etax*nu*psiy*psiy)/15 - (4*etay*etay*nu*psix*psix)/15 + (8*etax*etay*psix*psiy)/15 + (8*etax*etay*nu*psix*psiy)/15);
	Matrix[76] = (k*(etax*psix + etay*psiy)*(etax*etax + etay*etay + psix*psix + psiy*psiy))/6;
	Matrix[77] = k*(etax*psix*psix*psix + etay*psiy*psiy*psiy - psix*psix*psix*psix/2 - psiy*psiy*psiy*psiy/2 + (7*etax*etax*psix*psix)/10 + (etax*etax*psiy*psiy)/10 + (etay*etay*psix*psix)/10 + (7*etay*etay*psiy*psiy)/10 - psix*psix*psiy*psiy + etax*psix*psiy*psiy + etay*psix*psix*psiy + (2*etax*etax*nu*psiy*psiy)/5 + (2*etay*etay*nu*psix*psix)/5 + (6*etax*etay*psix*psiy)/5 - (4*etax*etay*nu*psix*psiy)/5);


#ifdef _DEBUG_
	for (unsigned int i=0;i<78;++i)
		{
			cout<<"Matrix element "<<i<<" is "<<Matrix[i]<<endl;
	};
#endif


}

//	Calculate element stress 
void CPlate::ElementStress(double* stress, double* Displacement,double* position)
{
	CPlateMaterial* material = dynamic_cast<CPlateMaterial*>(ElementMaterial);	// Pointer to material of the element

	//remember to change the size of stress
	double xpsi=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])*0.5;
	double xeta=(nodes[1]->XYZ[1]-nodes[0]->XYZ[1])*0.5;
	double ypsi=(nodes[3]->XYZ[0]-nodes[0]->XYZ[0])*0.5;
	double yeta=(nodes[3]->XYZ[1]-nodes[0]->XYZ[1])*0.5;
	double dis[12];//displacement of nodes
	for (unsigned i=0;i<=11;i++){
		if (LocationMatrix[i])
		dis[i+1]=Displacement[LocationMatrix[i]-1];
		else
			dis[i]=0.0;

	};
	

	double nu=material->nu;
	double Jacobian=xpsi*yeta-ypsi*xeta;
	const double k= material->E*material->h*material->h*material->h*material->h*0.5*0.0833333333333333/(1-nu*nu);

	double psix=yeta/Jacobian;
	double psiy=-xeta/Jacobian;
	double etax=-ypsi/Jacobian;
	double etay=xpsi/Jacobian;

	/*well it's just too complex a set of formulas...*/
stress[0] = dis[0]*((psix*psix+nu*psiy*psiy)*  0.6830127018922193+(etax*etax+nu*etay*etay)*  0.6830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[1]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  1.0773502691896257+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[2]*((psix*psix+nu*psiy*psiy)* -1.0773502691896257+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[3]*((psix*psix+nu*psiy*psiy)* -0.6830127018922193+(etax*etax+nu*etay*etay)*  0.1830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[4]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  0.2886751345948129+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[5]*((psix*psix+nu*psiy*psiy)* -0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[6]*((psix*psix+nu*psiy*psiy)* -0.1830127018922193+(etax*etax+nu*etay*etay)* -0.1830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[7]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  0.0773502691896258+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[8]*((psix*psix+nu*psiy*psiy)* -0.0773502691896258+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[9]*((psix*psix+nu*psiy*psiy)*  0.1830127018922193+(etax*etax+nu*etay*etay)* -0.6830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[10]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  0.2886751345948129+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[11]*((psix*psix+nu*psiy*psiy)* -0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129);
stress[1] = dis[0]*((nu*psix*psix+psiy*psiy)*  0.6830127018922193+(nu*etax*etax+etay*etay)*  0.6830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[1]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  1.0773502691896257+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[2]*((nu*psix*psix+psiy*psiy)* -1.0773502691896257+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[3]*((nu*psix*psix+psiy*psiy)* -0.6830127018922193+(nu*etax*etax+etay*etay)*  0.1830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[4]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  0.2886751345948129+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[5]*((nu*psix*psix+psiy*psiy)* -0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[6]*((nu*psix*psix+psiy*psiy)* -0.1830127018922193+(nu*etax*etax+etay*etay)* -0.1830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[7]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  0.0773502691896258+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[8]*((nu*psix*psix+psiy*psiy)* -0.0773502691896258+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[9]*((nu*psix*psix+psiy*psiy)*  0.1830127018922193+(nu*etax*etax+etay*etay)* -0.6830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[10]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  0.2886751345948129+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[11]*((nu*psix*psix+psiy*psiy)* -0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129);
stress[2] = dis[0]*(((1-nu)*psix*psiy)*  0.6830127018922193+((1-nu)*etax*etay)*  0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  1.0773502691896257+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[2]*(((1-nu)*psix*psiy)* -1.0773502691896257+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[3]*(((1-nu)*psix*psiy)* -0.6830127018922193+((1-nu)*etax*etay)*  0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[5]*(((1-nu)*psix*psiy)* -0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[6]*(((1-nu)*psix*psiy)* -0.1830127018922193+((1-nu)*etax*etay)* -0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  0.0773502691896258+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[8]*(((1-nu)*psix*psiy)* -0.0773502691896258+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[9]*(((1-nu)*psix*psiy)*  0.1830127018922193+((1-nu)*etax*etay)* -0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[11]*(((1-nu)*psix*psiy)* -0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129);
stress[3] = dis[0]*((psix*psix+nu*psiy*psiy)* -0.6830127018922193+(etax*etax+nu*etay*etay)*  0.1830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[1]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  0.2886751345948129+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[2]*((psix*psix+nu*psiy*psiy)*  0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[3]*((psix*psix+nu*psiy*psiy)*  0.6830127018922193+(etax*etax+nu*etay*etay)*  0.6830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[4]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  1.0773502691896257+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[5]*((psix*psix+nu*psiy*psiy)*  1.0773502691896257+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[6]*((psix*psix+nu*psiy*psiy)*  0.1830127018922193+(etax*etax+nu*etay*etay)* -0.6830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[7]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  0.2886751345948129+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[8]*((psix*psix+nu*psiy*psiy)*  0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[9]*((psix*psix+nu*psiy*psiy)* -0.1830127018922193+(etax*etax+nu*etay*etay)* -0.1830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[10]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)*  0.0773502691896258+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[11]*((psix*psix+nu*psiy*psiy)*  0.0773502691896258+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129);
stress[4] = dis[0]*((nu*psix*psix+psiy*psiy)* -0.6830127018922193+(nu*etax*etax+etay*etay)*  0.1830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[1]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  0.2886751345948129+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[2]*((nu*psix*psix+psiy*psiy)*  0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[3]*((nu*psix*psix+psiy*psiy)*  0.6830127018922193+(nu*etax*etax+etay*etay)*  0.6830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[4]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  1.0773502691896257+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[5]*((nu*psix*psix+psiy*psiy)*  1.0773502691896257+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[6]*((nu*psix*psix+psiy*psiy)*  0.1830127018922193+(nu*etax*etax+etay*etay)* -0.6830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[7]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  0.2886751345948129+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[8]*((nu*psix*psix+psiy*psiy)*  0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[9]*((nu*psix*psix+psiy*psiy)* -0.1830127018922193+(nu*etax*etax+etay*etay)* -0.1830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[10]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)*  0.0773502691896258+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[11]*((nu*psix*psix+psiy*psiy)*  0.0773502691896258+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129);
stress[5] = dis[0]*(((1-nu)*psix*psiy)* -0.6830127018922193+((1-nu)*etax*etay)*  0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[2]*(((1-nu)*psix*psiy)*  0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[3]*(((1-nu)*psix*psiy)*  0.6830127018922193+((1-nu)*etax*etay)*  0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  1.0773502691896257+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[5]*(((1-nu)*psix*psiy)*  1.0773502691896257+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[6]*(((1-nu)*psix*psiy)*  0.1830127018922193+((1-nu)*etax*etay)* -0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[8]*(((1-nu)*psix*psiy)*  0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[9]*(((1-nu)*psix*psiy)* -0.1830127018922193+((1-nu)*etax*etay)* -0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)*  0.0773502691896258+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[11]*(((1-nu)*psix*psiy)*  0.0773502691896258+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129);
stress[6] = dis[0]*((psix*psix+nu*psiy*psiy)* -0.1830127018922193+(etax*etax+nu*etay*etay)* -0.1830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[1]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -0.0773502691896258+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[2]*((psix*psix+nu*psiy*psiy)*  0.0773502691896258+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[3]*((psix*psix+nu*psiy*psiy)*  0.1830127018922193+(etax*etax+nu*etay*etay)* -0.6830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[4]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -0.2886751345948129+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[5]*((psix*psix+nu*psiy*psiy)*  0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[6]*((psix*psix+nu*psiy*psiy)*  0.6830127018922193+(etax*etax+nu*etay*etay)*  0.6830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[7]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -1.0773502691896257+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[8]*((psix*psix+nu*psiy*psiy)*  1.0773502691896257+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[9]*((psix*psix+nu*psiy*psiy)* -0.6830127018922193+(etax*etax+nu*etay*etay)*  0.1830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[10]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -0.2886751345948129+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[11]*((psix*psix+nu*psiy*psiy)*  0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129);
stress[7] = dis[0]*((nu*psix*psix+psiy*psiy)* -0.1830127018922193+(nu*etax*etax+etay*etay)* -0.1830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[1]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -0.0773502691896258+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[2]*((nu*psix*psix+psiy*psiy)*  0.0773502691896258+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[3]*((nu*psix*psix+psiy*psiy)*  0.1830127018922193+(nu*etax*etax+etay*etay)* -0.6830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[4]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -0.2886751345948129+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[5]*((nu*psix*psix+psiy*psiy)*  0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[6]*((nu*psix*psix+psiy*psiy)*  0.6830127018922193+(nu*etax*etax+etay*etay)*  0.6830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[7]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -1.0773502691896257+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[8]*((nu*psix*psix+psiy*psiy)*  1.0773502691896257+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[9]*((nu*psix*psix+psiy*psiy)* -0.6830127018922193+(nu*etax*etax+etay*etay)*  0.1830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[10]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -0.2886751345948129+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[11]*((nu*psix*psix+psiy*psiy)*  0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129);
stress[8] = dis[0]*(((1-nu)*psix*psiy)* -0.1830127018922193+((1-nu)*etax*etay)* -0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -0.0773502691896258+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[2]*(((1-nu)*psix*psiy)*  0.0773502691896258+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[3]*(((1-nu)*psix*psiy)*  0.1830127018922193+((1-nu)*etax*etay)* -0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[5]*(((1-nu)*psix*psiy)*  0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[6]*(((1-nu)*psix*psiy)*  0.6830127018922193+((1-nu)*etax*etay)*  0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -1.0773502691896257+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[8]*(((1-nu)*psix*psiy)*  1.0773502691896257+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[9]*(((1-nu)*psix*psiy)* -0.6830127018922193+((1-nu)*etax*etay)*  0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[11]*(((1-nu)*psix*psiy)*  0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129);
stress[9] = dis[0]*((psix*psix+nu*psiy*psiy)*  0.1830127018922193+(etax*etax+nu*etay*etay)* -0.6830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[1]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -0.2886751345948129+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[2]*((psix*psix+nu*psiy*psiy)* -0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[3]*((psix*psix+nu*psiy*psiy)* -0.1830127018922193+(etax*etax+nu*etay*etay)* -0.1830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[4]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -0.0773502691896258+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[5]*((psix*psix+nu*psiy*psiy)* -0.0773502691896258+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[6]*((psix*psix+nu*psiy*psiy)* -0.6830127018922193+(etax*etax+nu*etay*etay)*  0.1830127018922193+(psix*etax+nu*psiy*etay)* -0.5000000000000000)+dis[7]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -0.2886751345948129+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[8]*((psix*psix+nu*psiy*psiy)* -0.2886751345948129+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)* -0.2886751345948129)+dis[9]*((psix*psix+nu*psiy*psiy)*  0.6830127018922193+(etax*etax+nu*etay*etay)*  0.6830127018922193+(psix*etax+nu*psiy*etay)*  0.5000000000000000)+dis[10]*((psix*psix+nu*psiy*psiy)*  0.0000000000000000+(etax*etax+nu*etay*etay)* -1.0773502691896257+(psix*etax+nu*psiy*etay)*  0.2886751345948129)+dis[11]*((psix*psix+nu*psiy*psiy)* -1.0773502691896257+(etax*etax+nu*etay*etay)*  0.0000000000000000+(psix*etax+nu*psiy*etay)*  0.2886751345948129);
stress[10] = dis[0]*((nu*psix*psix+psiy*psiy)*  0.1830127018922193+(nu*etax*etax+etay*etay)* -0.6830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[1]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -0.2886751345948129+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[2]*((nu*psix*psix+psiy*psiy)* -0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[3]*((nu*psix*psix+psiy*psiy)* -0.1830127018922193+(nu*etax*etax+etay*etay)* -0.1830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[4]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -0.0773502691896258+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[5]*((nu*psix*psix+psiy*psiy)* -0.0773502691896258+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[6]*((nu*psix*psix+psiy*psiy)* -0.6830127018922193+(nu*etax*etax+etay*etay)*  0.1830127018922193+(nu*psix*etax+psiy*etay)* -0.5000000000000000)+dis[7]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -0.2886751345948129+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[8]*((nu*psix*psix+psiy*psiy)* -0.2886751345948129+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)* -0.2886751345948129)+dis[9]*((nu*psix*psix+psiy*psiy)*  0.6830127018922193+(nu*etax*etax+etay*etay)*  0.6830127018922193+(nu*psix*etax+psiy*etay)*  0.5000000000000000)+dis[10]*((nu*psix*psix+psiy*psiy)*  0.0000000000000000+(nu*etax*etax+etay*etay)* -1.0773502691896257+(nu*psix*etax+psiy*etay)*  0.2886751345948129)+dis[11]*((nu*psix*psix+psiy*psiy)* -1.0773502691896257+(nu*etax*etax+etay*etay)*  0.0000000000000000+(nu*psix*etax+psiy*etay)*  0.2886751345948129);
stress[11] = dis[0]*(((1-nu)*psix*psiy)*  0.1830127018922193+((1-nu)*etax*etay)* -0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[1]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[2]*(((1-nu)*psix*psiy)* -0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[3]*(((1-nu)*psix*psiy)* -0.1830127018922193+((1-nu)*etax*etay)* -0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[4]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -0.0773502691896258+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[5]*(((1-nu)*psix*psiy)* -0.0773502691896258+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[6]*(((1-nu)*psix*psiy)* -0.6830127018922193+((1-nu)*etax*etay)*  0.1830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.5000000000000000)+dis[7]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -0.2886751345948129+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[8]*(((1-nu)*psix*psiy)* -0.2886751345948129+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5* -0.2886751345948129)+dis[9]*(((1-nu)*psix*psiy)*  0.6830127018922193+((1-nu)*etax*etay)*  0.6830127018922193+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.5000000000000000)+dis[10]*(((1-nu)*psix*psiy)*  0.0000000000000000+((1-nu)*etax*etay)* -1.0773502691896257+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129)+dis[11]*(((1-nu)*psix*psiy)* -1.0773502691896257+((1-nu)*etax*etay)*  0.0000000000000000+(psix*etax+psiy*etay)*(1-nu)*0.5*  0.2886751345948129);

	double xmid=(nodes[0]->XYZ[0]+nodes[1]->XYZ[0]+nodes[2]->XYZ[0]+nodes[3]->XYZ[0])*0.25;
	double ymid=(nodes[0]->XYZ[1]+nodes[1]->XYZ[1]+nodes[2]->XYZ[1]+nodes[3]->XYZ[1])*0.25;
	double zmid=(nodes[0]->XYZ[2]+nodes[1]->XYZ[2]+nodes[2]->XYZ[2]+nodes[3]->XYZ[2])*0.25;
	position[0]=xmid+sqrt(1.0/3)*(-xeta-xpsi);
	position[1]=ymid+sqrt(1.0/3)*(-yeta-ypsi);
	position[2]=zmid;
	position[3]=xmid+sqrt(1.0/3)*(-xeta+xpsi);
	position[4]=ymid+sqrt(1.0/3)*(-yeta+ypsi);
	position[5]=zmid;
	position[6]=xmid+sqrt(1.0/3)*(xeta+xpsi);
	position[7]=ymid+sqrt(1.0/3)*(yeta+ypsi);
	position[8]=zmid;
	position[9]=xmid+sqrt(1.0/3)*(xeta-xpsi);
	position[10]=ymid+sqrt(1.0/3)*(yeta-ypsi);
	position[11]=zmid;

}
