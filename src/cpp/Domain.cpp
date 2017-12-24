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
#include <Eigen/Dense>
#include "mkl.h"

using namespace std;
using namespace Eigen;

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
#ifdef _VIB_
	MassMatrix = nullptr;
#endif
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
#ifdef _VIB_
	delete MassMatrix;
#endif
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

#ifdef _VIB_
	if (ReadVibNum())
		Output->PrintVibModNum();
	else
		return false;
#endif
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

    GenerateLocationMatrix();

#ifdef MKL
    CSRStiffnessMatrix = new CSRMatrix<double>(NEQ);
    CalculateCSRColumns();
#else
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
#ifdef _VIB_
	MassMatrix = new CSkylineMatrix<double>(NEQ);
#endif
    //	Calculate column heights
	CalculateColumnHeights();
	//	Calculate address of diagonal elements in banded matrix
	CalculateDiagnoalAddress();
#endif

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

#ifdef _VIB_
	unsigned int* ColumnHeights_mass = MassMatrix->GetColumnHeights();
	for (unsigned int i=0; i<NEQ; ++i){
		ColumnHeights_mass[i]=ColumnHeights[i];
	}
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

#ifdef _VIB_
	unsigned int* DiagonalAddress_mass = MassMatrix->GetDiagonalAddress();
	for (unsigned int i=0; i<=NEQ; ++i){
		DiagonalAddress_mass[i]=DiagonalAddress[i];
	}
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

//  Create the banded stiffness matrix
#ifdef MKL
	GetCSRStiffnessMatrix().allocate();
#else
	StiffnessMatrix->Allocate();
#endif

#ifdef _VIB_
	MassMatrix->Allocate();
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
                    if (index1 < index2) {
                        matrix.markPosition(index1, index2);
                    }
                    else 
                    {
                        matrix.markPosition(index2, index1);
                    }
                }
            }
        }
    }
}

#ifdef _VIB_

void CDomain::AssembleMassMatrix()
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
			ElementGrp.GetElement(Ele).assembly_mass(Matrix, MassMatrix);

		delete[] Matrix;
		Matrix = nullptr;
	}
/*
#ifdef _DEBUG_
	COutputter::Instance()->PrintStiffnessMatrix();
#endif
*/
}

bool CDomain::VibSolver(unsigned int NVibModes){
    unsigned int numEq = GetNEQ();
    //an inadvanced way to get the origin vectors
	double* vector_o = new double[numEq*NVibModes];
	for (unsigned int i=0; i<numEq; ++i){
		for (unsigned j=0; j<NVibModes; ++j){
			if (i==0 || j==i-1) {
				vector_o[j*numEq+i]=1.0;
			}
			else {
				vector_o[j*numEq+i]=0.0;
			}
		}
	}
    double* K_reduced = new double[NVibModes*(NVibModes+1)/2];

	double* vector_oY= new double[numEq*NVibModes];
	double* M_reduced = new double[NVibModes*(NVibModes+1)/2];
	CLDLTSolver* KSolver= new CLDLTSolver(*StiffnessMatrix);
	CLDLTSolver* MSolver= new CLDLTSolver(*MassMatrix);
	MSolver->Multiple(vector_o,vector_oY,numEq,NVibModes);
	KSolver->LDLT();
	double* vec_x_n = new double[NVibModes*numEq];
	double* vec_y_n = new double[NVibModes*numEq];
	double* lambdas_this = new double[NVibModes];
	double* lambdas_n = new double[NVibModes];
	double* eig_vecs = new double[NVibModes*NVibModes];
	double l_diffs = 1.0;
    while (l_diffs > 1.19e-7){
		for (unsigned int i=0; i<NVibModes; ++i){
			lambdas_this[i] = lambdas_n[i];
		}
		for (unsigned int i=0; i<NVibModes*numEq; ++i){
			vec_x_n[i] = vector_oY[i];
		}
		for (unsigned int i=0; i<NVibModes; ++i){
			KSolver->BackSubstitution(vec_x_n+i*numEq);
		}
		MSolver->Multiple(vec_x_n,vec_y_n,numEq,NVibModes);
		int ord=NVibModes;
		for (unsigned int i=0; i<NVibModes; ++i){
			for (unsigned int j=i; j<NVibModes; ++j){
				K_reduced[i+j*(j+1)/2] = 0.0;
				M_reduced[i+j*(j+1)/2] = 0.0;
				for (unsigned int k=0; k<numEq; ++k){
					K_reduced[i+j*(j+1)/2] += vec_x_n[i*numEq+k] * vector_oY[j*numEq+k];
					M_reduced[i+j*(j+1)/2] += vec_x_n[i*numEq+k] * vec_y_n[i*numEq+k];
				}
			}
		}
		//LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', ord,M_reduced, ord);
		if (LAPACKE_dspgvd(LAPACK_COL_MAJOR, 1,'V','U',ord,K_reduced,M_reduced,lambdas_n,eig_vecs,ord))
			return false;
		for (unsigned int i=0; i<numEq; ++i){
			for (unsigned int j=0; j<NVibModes; ++j){
				vec_x_n[i+j*numEq] = 0.0;
				for (unsigned int k=0; k<NVibModes; ++k){
					vec_x_n[i+j*numEq] += vector_o[i+k*numEq] * eig_vecs[k+j*numEq];
				}
			}
		}
	    for (unsigned int i=0; i<NVibModes*numEq; ++i){
			vector_o[i] = vec_x_n[i];
		}
		MSolver->Multiple(vec_x_n,vec_y_n,numEq,NVibModes);

	}
	VibDisp = new double[NVibModes*numEq];
	for (unsigned int i=0; i<NVibModes*numEq; ++i){
		VibDisp[i] = vector_o[i];
	}
	for (unsigned int i=0; i<NVibModes; ++i){
		EigenValues[i] = lambdas_n[i];
	}
	delete vector_o;
	delete vector_oY;
	delete vec_x_n;
	delete vec_y_n;
	delete lambdas_n;
	delete lambdas_this;
	delete K_reduced;
	delete M_reduced;
	delete eig_vecs;
	delete KSolver;
	delete MSolver;
	return true;
}

bool CDomain::ReadVibNum() {
	Input >> numEig;
	return true;
}
#endif