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
	getline(Input, NodeInfo, '\n');

	// Replace all spaces with tabs
	int spaceTempPos = 0;
	int spaceNextPos = NodeInfo.find(' ', spaceTempPos);
	while(spaceNextPos >= 0) {
		NodeInfo.replace(spaceNextPos, 1, "\t");
		spaceTempPos = spaceNextPos;
		spaceNextPos = NodeInfo.find(' ', spaceTempPos);
	}

	// Number of NodalInfos, which is computed from the number of tab blocks
	int tabBlockNum = 1;
	// Delete blocked tabs
	int tabTempPos = 0;
	int tabNextPos = NodeInfo.find("\t", tabTempPos + 1);
	while (tabNextPos >= 0) {
		if (tabNextPos - tabTempPos == 1) {
			NodeInfo.erase(tabNextPos, 1);
		}
		else {
			tabTempPos = tabNextPos;
			tabBlockNum++;
		}
		tabNextPos = NodeInfo.find("\t", tabTempPos + 1);
	}

	// Determine the input format:
	//     While the last 3 bcodes are manually input, tabBlockNum is 9;
	//     While the default values are chosen, tabBlockNum is 6.
	// Default values of the last 3 bcodes (related to the rotation):
	//     Structure elements: active,     value = 0;
	//     Solid elements:     not active, value = 1.
	if (tabBlockNum == 9) {
		// Rewrite the flag marking whether the rotation DOF are manually input
		RotationDOFManuallyInputFlag = 1;
		// Save the nodal infos to bcode[1:6] and XYZ[]
		int posTab = 0;
		string tempStr;
		for (unsigned int i = 0; i < 6; i++) {
			int posNextTab = NodeInfo.find('\t', posTab + 1);
			tempStr.assign(NodeInfo, posTab, posNextTab - posTab);
			posTab += posNextTab - posTab;
			cout << tempStr << endl;
			bcode[i] = atoi(tempStr.c_str());
		}
		for (unsigned int i = 0; i < 2; i++) {
			int posNextTab = NodeInfo.find('\t',posTab+1);
			tempStr.assign(NodeInfo, posTab, posNextTab - posTab);
			posTab += posNextTab - posTab;
			XYZ[i] = atof(tempStr.c_str());
		}
		int LastStrLength = NodeInfo.at(NodeInfo.size() - 1) - posTab;
		tempStr.assign(NodeInfo, posTab, LastStrLength);
		XYZ[2] = atof(tempStr.c_str());
	}
	else if (tabBlockNum == 6) {
		// Save the nodal infos to bcode[1:3] and XYZ[]
		int posTab = 0;
		string tempStr;
		for (unsigned int i = 0; i < 3; i++) {
			int posNextTab = NodeInfo.find('\t', posTab + 1);
			tempStr.assign(NodeInfo, posTab, posNextTab - posTab);
			posTab += posNextTab - posTab;
			cout << tempStr << endl;
			bcode[i] = atoi(tempStr.c_str());
		}
		for (unsigned int i = 0; i < 2; i++) {
			int posNextTab = NodeInfo.find('\t', posTab + 1);
			tempStr.assign(NodeInfo, posTab, posNextTab - posTab);
			posTab += posNextTab - posTab;
			XYZ[i] = atof(tempStr.c_str());
		}
		int LastStrLength = NodeInfo.at(NodeInfo.size() - 1) - posTab;
		tempStr.assign(NodeInfo, posTab, LastStrLength);
		XYZ[2] = atof(tempStr.c_str());
	}
	else {
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
