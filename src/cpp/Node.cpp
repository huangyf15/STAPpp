/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <string>

#include "Node.h"

CNode::CNode(double X, double Y, double Z):NodeNumber(0)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
    bcode[3] = 1;
    bcode[4] = 1;
    bcode[5] = 1;

    RotationDOFManuallyInputFlag = 0;	// Boundary code flag
};

// return total count of non-blank args in string
int getArgsCount(std::string buff)
{
    int count = 0;
    bool inNumFlag = false;
    for (unsigned i=0; i<buff.length(); ++i)
    {
        if (buff[i] != ' ' && buff[i] != '\t')
        {
            if (!inNumFlag)
            {
                inNumFlag = true;
                count++;
            }
        }
        else
        {
            inNumFlag = false;
        }
    }
    return count;
}

//	Read element data from stream Input
bool CNode::Read(ifstream& Input, unsigned int np)
{
	unsigned int N;

	Input >> N;	// node number
	if (N != np + 1) 
	{
		cerr << "*** Error *** Nodes must be inputted in order !" << endl 
			 << "   Expected node number : " << np + 1 << endl
			 << "   Provided node number : " << N << endl;

		return false;
	}

	NodeNumber = N;

	// Read the dataline
	string NodeInfo, ModiNodeInfo;
	getline(Input, NodeInfo);

	int tabBlockNum = getArgsCount(NodeInfo);

	// Determine the input format:
	//     While the last 3 bcodes are manually input, tabBlockNum is 9;
	//     While the default values are chosen, tabBlockNum is 6.
	// Default values of the last 3 bcodes (related to the rotation):
	//     Structure elements: active,     value = 0;
	//     Solid elements:     not active, value = 1.
	
	// Rewrite the flag marking whether the rotation DOF are manually input
	RotationDOFManuallyInputFlag = (tabBlockNum == 9);
	// Save the nodal infos to bcode[0:5] and XYZ[]
	if (tabBlockNum == 9)
	{
		sscanf(NodeInfo.c_str(), "%d%d%d%d%d%d%lf%lf%lf",
			bcode, bcode+1, bcode+2,
			bcode+3, bcode+4, bcode+5,
			XYZ, XYZ+1, XYZ+2);
	}
	else if (tabBlockNum == 6)
	{
		sscanf(NodeInfo.c_str(), "%d%d%d%lf%lf%lf",
			bcode, bcode+1, bcode+2,
			XYZ, XYZ+1, XYZ+2);
	}
	else
	{
		cerr << "*** Error *** NodeInfos must be inputted in the correct format! " << endl
			<< "  Present Number of Nodeinfos: " << tabBlockNum << endl
			<< "  Correct Number of Nodeinfos: 6 or 9 !" << endl;
		return false;
	}

	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output, unsigned int np)
{
	output << setw(9) << np + 1 
		<< setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		<< setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output, unsigned int np)
{
	output << setw(9) << np+1 << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, unsigned int np, double* Displacement)
{
	output << setw(5) << np + 1 << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}

#ifdef _VIB_
void CNode::WriteVibrationDisplacement(COutputter& output, unsigned int np, double* Displacement)
{
	output << setw(5) << np + 1 << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}

#endif