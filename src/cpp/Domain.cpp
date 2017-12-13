/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;
	NWK = 0;
	MK = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
	CSRStiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
	delete CSRStiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::Instance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::Instance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

    //	Read load data
    if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

    //	Read element data
    if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

    //	Update equation number
    CalculateEquationNumber();
    Output->OutputEquationNumber();

    return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{
	//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

	//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}

void CDomain::GenerateLocationMatrix()
{
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			ElementGrp.GetElement(Ele).GenerateLocationMatrix();
    }
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	// Adjust bcode for nodes of structure elements
	// Default values of the last 3 bcodes (related to the rotation):
	//     Structure elements: active,     value = 0;
	//     Solid elements:     not active, value = 1.
    for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
    {
        ElementTypes eleType = EleGrpList[EleGrp].GetElementType();
        if (eleType == ElementTypes::Beam ||
            eleType == ElementTypes::TimoshenkoEBMOD ||
            eleType == ElementTypes::TimoshenkoSRINT)
        {
            const unsigned int NumE = EleGrpList[EleGrp].GetNUME();
            for (unsigned int NumEle = 0; NumEle < NumE; NumEle++)
            {
                const unsigned int NEN = EleGrpList[EleGrp].GetElement(NumEle).GetNEN();
                CNode** ElementNode = EleGrpList[EleGrp].GetElement(NumEle).GetNodes();
                for (unsigned int NumNode = 0; NumNode < NEN; NumNode++)
                {
                    if (!ElementNode[NumNode]->RotationDOFManuallyInputFlag)
                    {
                        const unsigned int N = ElementNode[NumNode]->NodeNumber;
                        NodeList[N - 1].bcode[3] = 0;
                        NodeList[N - 1].bcode[4] = 0;
                        NodeList[N - 1].bcode[5] = 0;
                    }   
                }
            }
        }
        else if (eleType == ElementTypes::Shell || 
                 eleType == ElementTypes::Plate)
        {
            const unsigned int NumE = EleGrpList[EleGrp].GetNUME();
            for (unsigned int NumEle = 0; NumEle < NumE; NumEle++)
            {
                const unsigned int NEN = EleGrpList[EleGrp].GetElement(NumEle).GetNEN();
                CNode** ElementNode = EleGrpList[EleGrp].GetElement(NumEle).GetNodes();
                bool x_inpl = true; //judge if all points have the same x 
                bool y_inpl = true;
                bool z_inpl = true;
                for (unsigned int NumNode = 1; NumNode < NEN; NumNode++)
                {
                    if (x_inpl && abs(ElementNode[NumNode]->XYZ[0] - ElementNode[0]->XYZ[0]) > DBL_EPSILON)
                    {
                        x_inpl = false;
                    }
                    if (y_inpl && abs(ElementNode[NumNode]->XYZ[1] - ElementNode[0]->XYZ[1]) > DBL_EPSILON)
                    {
                        y_inpl = false;
                    }
                    if (z_inpl && abs(ElementNode[NumNode]->XYZ[2] - ElementNode[0]->XYZ[2]) > DBL_EPSILON)
                    {
                        z_inpl = false;
                    }
                }
                for (unsigned int NumNode = 0; NumNode < NEN; NumNode++)
                {
                    if (!ElementNode[NumNode]->RotationDOFManuallyInputFlag)
                    {
                        const unsigned int N = ElementNode[NumNode]->NodeNumber;
                        if (x_inpl)
                        {
                            NodeList[N - 1].bcode[4] = 0;
                            NodeList[N - 1].bcode[5] = 0;
                        }
                        else if (y_inpl)
                        {
                            NodeList[N - 1].bcode[3] = 0;
                            NodeList[N - 1].bcode[5] = 0;
                        }
                        else if (z_inpl)
                        {
                            NodeList[N - 1].bcode[3] = 0;
                            NodeList[N - 1].bcode[4] = 0;
                        }
                        else
                        {
                            NodeList[N - 1].bcode[3] = 0;
                            NodeList[N - 1].bcode[4] = 0;
                            NodeList[N - 1].bcode[5] = 0;
                        }
                    }   

                }
            }
        }
    }

	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
		if (!LoadCases[lcase].Read(Input, lcase))
			return false;

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
    unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			ElementGrp.GetElement(Ele).CalculateColumnHeight(ColumnHeights);
    }

	//	Maximum half bandwidth ( = max(ColumnHeights) + 1 )
	MK = ColumnHeights[0];

	for (unsigned int i=1; i<NEQ; i++)
		if (MK < ColumnHeights[i])
			MK = ColumnHeights[i];

	MK = MK + 1;

#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintColumnHeights();
#endif

}

//	Calculate address of diagonal elements in banded matrix
//	Caution: Address is numbered from 1 !
void CDomain::CalculateDiagnoalAddress()
{
    unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();
    unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

//	Calculate the address of diagonal elements
//	M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
	DiagonalAddress[0] = 1;
	for (unsigned int col = 1; col <= NEQ; col++)
		DiagonalAddress[col] = DiagonalAddress[col - 1] + ColumnHeights[col-1] + 1;

//	Number of elements in banded global stiffness matrix
	NWK = DiagonalAddress[NEQ] - DiagonalAddress[0];

#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintDiagonalAddress();
#endif

}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp.GetElement(0).SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			ElementGrp.GetElement(Ele).assembly(Matrix, StiffnessMatrix, CSRStiffnessMatrix);

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter::Instance()->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	const CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

	//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
		{
			#ifdef MKL
				Force[dof - 1 + (NEQ*NLCASE*(LoadCase-1))] += LoadData->load[lnum];
			#else
				Force[dof - 1] += LoadData->load[lnum];
			#endif
		}
            
	}

	return true;
}

//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
//	Allocate for global force/displacement vector
#ifdef MKL
	Force = new double[NEQ * NLCASE];
	clear(Force, NEQ * NLCASE);
#else
	Force = new double[NEQ];
	clear(Force, NEQ);
#endif

	GenerateLocationMatrix();

//  Create the banded stiffness matrix
#ifdef MKL
	CSRStiffnessMatrix = new CSRMatrix<double>(NEQ);
	CalculateCSRColumns();
	GetCSRStiffnessMatrix().allocate();
#else
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
	//	Calculate column heights
	CalculateColumnHeights();
	//	Calculate address of diagonal elements in banded matrix
	CalculateDiagnoalAddress();
	StiffnessMatrix->Allocate();
#endif
	
	COutputter* Output = COutputter::Instance();
	Output->OutputTotalSystemData();
}

void CDomain::CalculateCSRColumns()
{
    auto& matrix = GetCSRStiffnessMatrix();
    matrix.beginPostionMark();

    for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
        for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            auto& element = ElementGrp.GetElement(Ele);
            unsigned LMSize = element.GetLMSize();
            unsigned* LM = element.GetLocationMatrix();
            for (unsigned i=0; i<LMSize; ++i)
            {
                unsigned index1 = LM[i];
                if (!index1) continue;
                for (unsigned j=i; j<LMSize; ++j)
                {
                    unsigned index2 = LM[j];
                    if (!index2) continue;
                    unsigned row = std::min(index1, index2);
                    unsigned column = std::max(index1, index2);
                    matrix.markPosition(row, column);
                }
            }
        }
    }
}
