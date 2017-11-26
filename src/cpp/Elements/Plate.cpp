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
	ElementMaterial = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	nodes[3] = &NodeList[N4 - 1];
	/*thinking about how to add rotating degrees of freedom*/
	return true;
}

void CPlate::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber 
		   << setw(9) << nodes[1]->NodeNumber << setw(9) << nodes[2]->NodeNumber 
		   << setw(9) << nodes[3]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}


/*to be continued*/