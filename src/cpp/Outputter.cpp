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
#include "Outputter.h"
#include "SkylineMatrix.h"

#include <iostream>
#include <iomanip>
#include <ctime>


using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;
COutputter2* COutputter2::_instance2 = nullptr;


//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

COutputter2::COutputter2(string FileName2)
{
	OutputFile2.open(FileName2);

	if (!OutputFile2)
	{
		cerr << "*** Error *** File " << FileName2 << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::Instance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
	return _instance;
}

//	Return the single instance of the class
COutputter2* COutputter2::Instance2(string FileName2)
{
	if (!_instance2)
		_instance2 = new COutputter2(FileName2);
	return _instance2;
}


//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::Instance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this, np);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::Instance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this, np);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::Instance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, ELEMENTS CURRENTLY" << endl
			  << "     EQ.3, NOT AVAILABLE" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				PrintBarElementData(EleGrp);
				break;
			case ElementTypes::Quadrilateral:
				PrintQuadrilateralElementData(EleGrp);
				break;
			case ElementTypes::Triangle:
				PrintTriangleElementData(EleGrp);
				break;
			case ElementTypes::Hexahedron:
				PrintHexElementData(EleGrp);
				break;
			case ElementTypes::Beam:
				PrintBeamElementData(EleGrp);
				break;
			case ElementTypes::TimoshenkoSRINT:
				PrintTimoshenkoSRINTElementData(EleGrp);
				break;
			case ElementTypes::TimoshenkoEBMOD:
				PrintTimoshenkoSRINTElementData(EleGrp);
				break;
			case ElementTypes::Plate:
				PrintPlateElementData(EleGrp);
				break;
			case ElementTypes::Shell:
				PrintShellElementData(EleGrp);
				break;

			default:
				std::cerr << "unknown ElementType " << ElementType << std::endl;
				exit(2);
				break;
		}
	}
}

//	Output bar element data
void COutputter::PrintBarElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	const unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}

//	Output quadrilateral element data
void COutputter::PrintQuadrilateralElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL CONSTANTS . . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S        POISSON'S" << endl
		  << " NUMBER     MODULUS          RATIO" << endl
		  << "               E              nu" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE      MATERIAL" << endl
		  << " NUMBER-N      I        J        K        L      SET NUMBER" << endl;

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}

// Output Beam element data
void COutputter::PrintBeamElementData(unsigned int EleGrp)
{
    CDomain* FEMData = CDomain::Instance();

    CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
    unsigned int NUMMAT = ElementGroup.GetNUMMAT();

    *this << " M A T E R I A L   D E F I N I T I O N" << endl << endl;
    *this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
    *this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
          << endl
          << endl;

    *this << "  SET       YOUNG'S       POISSON'S                               CROSS-SECTIONAL CONSTANTS [ONLY VALID FOR BOX]" << endl
          << " NUMBER     MODULUS         RATIO                 " << endl
          << "               E             nu               a               b               t1               t2               t3               t4" << endl;

    *this << setiosflags(ios::scientific) << setprecision(5);

    //	Loop over for all property sets
    for (unsigned int mset = 0; mset < NUMMAT; mset++)
        ElementGroup.GetMaterial(mset).Write(*this, mset);

    *this << endl << endl << " E L E M E N T   I N F O R M A T I O N" << endl;
    *this << " ELEMENT     NODE     NODE       MATERIAL" << endl
          << " NUMBER-N      I        J       SET NUMBER" << endl;

    const unsigned int NUME = ElementGroup.GetNUME();

    //	Loop over for all elements in group EleGrp
    for (unsigned int Ele = 0; Ele < NUME; Ele++)
        ElementGroup.GetElement(Ele).Write(*this, Ele);

    *this << endl;
}

// Output Triangle element data
void COutputter::PrintTriangleElementData(unsigned int EleGrp)

{
	
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND POISSON'S RATIO  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S        POISSON'S" << endl
		<< " NUMBER     MODULUS          RATIO" << endl
		<< "               E              nu" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		<< endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE        MATERIAL" << endl
		<< " NUMBER-N      I        J        K        SET NUMBER" << endl;

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}

//	Output TimoshenkoSRINT Beam element data
void COutputter::PrintTimoshenkoSRINTElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET      YOUNG'S      POISSON'S  CROSS-SECTIONAL MOMENT-INER  MOMENT-INER   1ST COORD     2ND COORD     3RD COORD" << endl
		<< " NUMBER    MODULUS        RATIO         AREA      LOCAL Y-AXIS LOCAL Z-AXIS  LOCAL Y-AXIS  LOCAL Y-AXIS  LOCAL Y-AXIS" << endl
		<< "              E             nu          Area          Iyy          Izz         THETAY1       THETAY2       THETAY3" << endl << endl;

	*this << setiosflags(ios::scientific) << setprecision(4);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		<< endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     MATERIAL" << endl
		<< " NUMBER-N      I        J        SET NUMBER" << endl;

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}
	  



//	Output hexahedron element data
void COutputter::PrintHexElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND POISSON'S RATIO  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S        POISSON'S" << endl
		  << " NUMBER     MODULUS          RATIO" << endl
		  << "               E              nu" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE      NODE     NODE     NODE     NODE     NODE        MATERIAL" << endl
		  << " NUMBER-N      1        2        3        4        5        6        7        8          SET NUMBER" << endl;

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}

//	Output TimoshenkoEBMOD Beam element data
void COutputter::PrintTimoshenkoEBMODElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET      YOUNG'S      POISSON'S  CROSS-SECTIONAL MOMENT-INER  MOMENT-INER   1ST COORD     2ND COORD     3RD COORD" << endl
		<< " NUMBER    MODULUS        RATIO         AREA      LOCAL Y-AXIS LOCAL Z-AXIS  LOCAL Y-AXIS  LOCAL Y-AXIS  LOCAL Y-AXIS" << endl
		<< "              E             nu          Area          Iyy          Izz         THETAY1       THETAY2       THETAY3" << endl << endl;

	*this << setiosflags(ios::scientific) << setprecision(4);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		<< endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     MATERIAL" << endl
		<< " NUMBER-N      I        J        SET NUMBER" << endl;

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::Instance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this, lcase);

		*this << endl;
	}
}

void COutputter::PrintPlateElementData(unsigned int EleGrp)
{
    CDomain* FEMData = CDomain::Instance();

    CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
    unsigned int NUMMAT = ElementGroup.GetNUMMAT();

    *this << " M A T E R I A L   D E F I N I T I O N" << endl << endl;
    *this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
    *this << " AND POISSON'S RATIO  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
          << endl
          << endl;

    *this << "  SET       YOUNG'S        THICKNESS        POISSON'S" << endl
          << " NUMBER     MODULUS         (HEIGHT)          RATIO" << endl
          << "               E               h                nu" << endl;

    *this << setiosflags(ios::scientific) << setprecision(5);

    //	Loop over for all property sets
    for (unsigned int mset = 0; mset < NUMMAT; mset++)
        ElementGroup.GetMaterial(mset).Write(*this, mset);

    *this << endl << endl << " E L E M E N T   I N F O R M A T I O N" << endl;
    *this << " ELEMENT     NODE     NODE     NODE     NODE      MATERIAL" << endl
          << " NUMBER-N      I        J        K        L      SET NUMBER" << endl;

    //	Loop over for all elements in group EleGrp
    for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
        ElementGroup.GetElement(Ele).Write(*this, Ele);

    *this << endl;
}

void COutputter::PrintShellElementData(unsigned int EleGrp)
{
    CDomain* FEMData = CDomain::Instance();

    CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
    unsigned int NUMMAT = ElementGroup.GetNUMMAT();

    *this << " M A T E R I A L   D E F I N I T I O N" << endl << endl;
    *this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
    *this << " AND POISSON'S RATIO  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
          << endl
          << endl;

    *this << "  SET       YOUNG'S        THICKNESS        POISSON'S" << endl
          << " NUMBER     MODULUS         (HEIGHT)          RATIO" << endl
          << "               E               h                nu" << endl;

    *this << setiosflags(ios::scientific) << setprecision(5);

    //	Loop over for all property sets
    for (unsigned int mset = 0; mset < NUMMAT; mset++)
        ElementGroup.GetMaterial(mset).Write(*this, mset);

    *this << endl << endl << " E L E M E N T   I N F O R M A T I O N" << endl;
    *this << " ELEMENT     NODE     NODE     NODE     NODE      MATERIAL" << endl
          << " NUMBER-N      I        J        K        L      SET NUMBER" << endl;

    //	Loop over for all elements in group EleGrp
    for (unsigned int Ele = 0; Ele < ElementGroup.GetNUME(); Ele++)
        ElementGroup.GetElement(Ele).Write(*this, Ele);

    *this << endl;
}


//	Print nodal displacement
void COutputter::OutputNodalDisplacement(unsigned int lcase)
{
	CDomain* FEMData = CDomain::Instance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << " LOAD CASE" << setw(5) << lcase + 1 << endl
		  << endl
		  << endl;

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT      X-ROTATION        Y-ROTATION        Z-ROTATION" << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, np, Displacement);

	*this << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::Instance();

	double* Displacement = FEMData->GetDisplacement();

	const unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress(&stress, Displacement);

					CBarMaterial& material = *static_cast<CBarMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
						<< stress << endl;
				}

				*this << endl;
				break;
			
			case ElementTypes::Quadrilateral: // Quadrilateral element
				*this << "    ELEMENT   GAUSS P           GUASS POINTS POSITIONS"
					<< "                       GUASS POINTS STRESSES"
					#ifdef _TEST_
					<< "                      GUASS POINTS DISPLACEMENTS            INTEGRATE"
					#endif
					<< endl;
				*this << "     NUMBER    INDEX        X             Y             Z" 
					<< "               SX'X'         SY'Y'        SX'Y'"
					#ifdef _TEST_
					<< "              UX            UY           UZ            WEIGHTS"
					#endif
					<< endl;
				double stresses[12];
				double Positions[12];
				#ifdef _TEST_
				double GaussDisplacements[12];
				double weights[4];
				#endif

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					#ifndef _TEST_
					static_cast<CQuadrilateral&>(
						EleGrp.GetElement(Ele)).ElementStress(stresses, Displacement, Positions);
					#else
					static_cast<CQuadrilateral&>(
						EleGrp.GetElement(Ele)).ElementStress(
							stresses, Displacement, Positions, GaussDisplacements, weights);
					#endif

					for (unsigned i=0; i<4; ++i) { // four gauss points
						*this << setw(8) << Ele + 1;
						*this << setw(10) << i+1;
						*this << setw(17) << Positions[i*3] << setw(14) << Positions[i*3+1] << setw(14) << Positions[i*3+2];
						*this << setw(17) << stresses[i*3] << setw(14) << stresses[i*3+1] << setw(14) << stresses[i*3+2];
						// *this << setw(32) << stresses[i] << std::endl;
						#ifdef _TEST_
						*this << setw(17) << GaussDisplacements[i*3] 
							<< setw(14) << GaussDisplacements[i*3+1] 
							<< setw(14) << GaussDisplacements[i*3+2];
						*this << setw(15) << weights[i];
						#endif
						*this << std::endl;
					}
				}
				*this << endl;

				break;

			case ElementTypes::Beam: // Beam element
				*this << "  ELEMENT          SXX                 SYY                   SZZ" << endl
					<< "  NUMBER" << endl;

				double beamstress[3];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress(beamstress, Displacement);

					CBeamMaterial& material =
						*static_cast<CBeamMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << beamstress[0] << setw(22)
						<< beamstress[1] << setw(22) << beamstress[2] << endl;
				}

				*this << endl;
				break;	

			case ElementTypes::Triangle: // 3T element
				double stress3T[3];
				#ifndef _TEST_
				*this << "  ELEMENT            LOCAL    ELEMENT    STRESS" << endl
					  << "  NUMBER         SXX            SYY            SXY" << endl;
				#else
				double GPPosition[9];
				double GPDisplacement[9];
				double weights3T[3];
				*this << "  ELEMENT    GP               GAUSS POINTS POSITION    "
					  << "                GAUSS POINTS DISPLACEMENTS       " 
					  << "               GAUSS POINTS STRESSES              INTEGRATE"
					  << std::endl
					  << "   INDEX   INDEX          X            Y              Z"
					  << "                DX           DY           DZ     "
					  << "          SXX           SYY           SXY          WEIGHTS"
					  << std::endl;
				#endif

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					#ifndef _TEST_
					static_cast<CTriangle&>(Element).ElementStress(stress3T, Displacement);
					#else
					static_cast<CTriangle&>(Element).ElementStress(stress3T, Displacement, GPPosition, GPDisplacement, weights3T);
					#endif
					CTriangleMaterial material = *static_cast<CTriangleMaterial*>(Element.GetElementMaterial());
					
					#ifndef _TEST_
					*this << setw(5) << Ele + 1 << setw(20) << stress3T[0] 
					      << setw(15) << stress3T[1] << setw(15) << stress3T[2] << endl;
					#else
					for (unsigned GPIndex=0; GPIndex<3; GPIndex++)
					{
						*this << setw(6) << Ele+1 << setw(8) << GPIndex+1 
							  << setw(18) << GPPosition[3*GPIndex] 
							  << setw(14) << GPPosition[3*GPIndex + 1] 
							  << setw(14) << GPPosition[3*GPIndex + 2]
							  << setw(17) << GPDisplacement[3*GPIndex]
							  << setw(14) << GPDisplacement[3*GPIndex + 1]
							  << setw(14) << GPDisplacement[3*GPIndex + 2]
							  << setw(17) << stress3T[0]
							  << setw(14) << stress3T[1]
							  << setw(14) << stress3T[2]
							  << setw(14) << weights3T[GPIndex]
							  << std::endl;
					}
					#endif
				}

				*this << endl;
				break;


			case ElementTypes::Hexahedron: // 8H element
				*this << "node      X              Y              Z              XY              YZ              XZ" << endl
					<< "NUMBER" << endl;

				double stressHex[48];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress(stressHex, Displacement);

					CHexMaterial& material = *static_cast<CHexMaterial*>(Element.GetElementMaterial());
					*this  << Ele + 1 << setw(15) << stressHex[0]<< setw(16) << stressHex[1] << setw(16)<< stressHex[2]<< setw(16)<< stressHex[3]<< setw(16)<< stressHex[4]<< setw(16)<< stressHex[5]<< endl
						<< setw(16) << stressHex[6]<< setw(16) << stressHex[7] << setw(16)<< stressHex[8]<< setw(16)<< stressHex[9]<< setw(16)<< stressHex[10]<< setw(16)<< stressHex[11]<< endl
						<< setw(16) << stressHex[12]<< setw(16) << stressHex[13] << setw(16)<< stressHex[14]<< setw(16)<< stressHex[15]<< setw(16)<< stressHex[16]<< setw(16)<< stressHex[17]<< endl
						<< setw(16) << stressHex[18]<< setw(16) << stressHex[19] << setw(16)<< stressHex[20]<< setw(16)<< stressHex[21]<< setw(16)<< stressHex[22]<< setw(16)<< stressHex[23]<< endl
						<< setw(16) << stressHex[24]<< setw(16) << stressHex[25] << setw(16)<< stressHex[26]<< setw(16)<< stressHex[27]<< setw(16)<< stressHex[28]<< setw(16)<< stressHex[29]<< endl
						<< setw(16) << stressHex[30]<< setw(16) << stressHex[31] << setw(16)<< stressHex[32]<< setw(16)<< stressHex[33]<< setw(16)<< stressHex[34]<< setw(16)<< stressHex[35]<< endl
						<< setw(16) << stressHex[36]<< setw(16) << stressHex[37] << setw(16)<< stressHex[38]<< setw(16)<< stressHex[39]<< setw(16)<< stressHex[40]<< setw(16)<< stressHex[41]<< endl
						<< setw(16) << stressHex[42]<< setw(16) << stressHex[43] << setw(16)<< stressHex[44]<< setw(16)<< stressHex[45]<< setw(16)<< stressHex[46]<< setw(16)<< stressHex[47]<< endl;
				}

				*this << endl;
				break;

			case ElementTypes::TimoshenkoSRINT: // TimoshenkoSRINT beam element
				double TimoshenkoStresses[3];
				double TimoshenkoForces[12];

				*this << "  ELEMENT        FORCE_X1    FORCE_X2    FORCE_Y1    FORCE_Y2    FORCE_Z1    FORCE_Z2   MOMENT_X1   MOMENT_X2   MOMENT_Y1   MOMENT_Y2   MOMENT_Z1   MOMENT_Z2  STRESS_XX  STRESS_YY  STRESS_XZ" << endl
					<< "  NUMBER" << endl;
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					static_cast<CTimoshenkoSRINT&>(
						EleGrp.GetElement(Ele)).ElementStress(TimoshenkoStresses, TimoshenkoForces, Displacement);
					*this << setw(5) << Ele + 1 << setw(18) << TimoshenkoForces[0] << setw(13) << TimoshenkoForces[1] 
						<< setw(13) << TimoshenkoForces[2] << setw(13) << TimoshenkoForces[3] << setw(13) << TimoshenkoForces[4]
						<< setw(13) << TimoshenkoForces[5] << setw(13) << TimoshenkoForces[6] << setw(13) << TimoshenkoForces[7]
						<< setw(13) << TimoshenkoForces[8] << setw(13) << TimoshenkoForces[9] << setw(13) << TimoshenkoForces[10]
						<< setw(13) << TimoshenkoForces[11] << setw(13) << TimoshenkoStresses[0] << setw(13) << TimoshenkoStresses[1] 
						<< setw(13) << TimoshenkoStresses[2] << endl;
					*this << std::endl;
				}
				break;

			case ElementTypes::TimoshenkoEBMOD: // TimoshenkoEBMOD beam element
				double TimoshenkoEBStresses[3];
				double TimoshenkoEBForces[12];

				*this << "  ELEMENT        FORCE_X1    FORCE_X2    FORCE_Y1    FORCE_Y2    FORCE_Z1    FORCE_Z2   MOMENT_X1   MOMENT_X2   MOMENT_Y1   MOMENT_Y2   MOMENT_Z1   MOMENT_Z2  STRESS_XX  STRESS_YY  STRESS_XZ" << endl
					<< "  NUMBER" << endl;
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					static_cast<CTimoshenkoEBMOD&>(
						EleGrp.GetElement(Ele)).ElementStress(TimoshenkoEBStresses, TimoshenkoEBForces, Displacement);
					*this << setw(5) << Ele + 1 << setw(18) << TimoshenkoEBForces[0] << setw(13) << TimoshenkoEBForces[1]
						<< setw(13) << TimoshenkoEBForces[2] << setw(13) << TimoshenkoEBForces[3] << setw(13) << TimoshenkoEBForces[4]
						<< setw(13) << TimoshenkoEBForces[5] << setw(13) << TimoshenkoEBForces[6] << setw(13) << TimoshenkoEBForces[7]
						<< setw(13) << TimoshenkoEBForces[8] << setw(13) << TimoshenkoEBForces[9] << setw(13) << TimoshenkoEBForces[10]
						<< setw(13) << TimoshenkoEBForces[11] << setw(13) << TimoshenkoEBStresses[0] << setw(13) << TimoshenkoEBStresses[1]
						<< setw(13) << TimoshenkoEBStresses[2] << endl;
					*this << std::endl;
				}
				break;
			case ElementTypes::Plate:
				*this << "    ELEMENT   GAUSS P           GUASS POINTS POSITIONS"
					<< "                       GUASS POINTS STRESSES"
					<< endl;
				*this << "     NUMBER    INDEX        X             Y             Z" 
					<< "               SX'X'_MAX     SY'Y'_MAX    SX'Y'_MAX"
					<< endl;
				double stresses4PE[12];
				double Positions4PE[12];
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					static_cast<CPlate&>(
						EleGrp.GetElement(Ele)).ElementStress(stresses4PE, Displacement, Positions4PE);
					
					for (unsigned i=0; i<4; ++i) { // four gauss points
						*this << setw(8) << Ele + 1;
						*this << setw(10) << i+1;
						*this << setw(17) << Positions4PE[i*3] << setw(14) << Positions4PE[i*3+1] << setw(14) << Positions4PE[i*3+2];
						*this << setw(17) << stresses4PE[i*3] << setw(14) << stresses4PE[i*3+1] << setw(14) << stresses4PE[i*3+2];
						// *this << setw(32) << stresses[i] << std::endl;
						
						*this << std::endl;
					}
				}
				*this << endl;

				break;
                        case ElementTypes::Shell:
                                *this << "    ELEMENT   GAUSS P           GUASS POINTS POSITIONS"
                                        << "                       GUASS POINTS STRESSES" << endl;
                                *this << "     NUMBER    INDEX        X             Y             Z"
                                        << "               SX'X'_MAX     SY'Y'_MAX    SX'Y'_MAX" << endl;
                                double stresses4SE[15];
                                double Positions4SE[15];
                                for (unsigned int Ele = 0; Ele < NUME; Ele++)
                                {

                                        dynamic_cast<CShell&>(EleGrp.GetElement(Ele))
                                                .ElementStress2(stresses4SE, Displacement, Positions4SE);

                                        for (unsigned i = 0; i < 5; ++i){
                                          // four gauss points;
                                          //THE FIFTH POINT IS THE CENTRE POINT FOR IN-PLANE STRESSES
                                                  *this << setw(8) << Ele + 1;
                                                  *this << setw(10) << i + 1;
                                                  *this << setw(17) << Positions4SE[i * 3] << setw(14) << Positions4SE[i * 3 + 1]
                                                        << setw(14) << Positions4SE[i * 3 + 2];
                                                  *this << setw(17) << stresses4SE[i * 3] << setw(14) << stresses4SE[i * 3 + 1]
                                                        << setw(14) << stresses4SE[i * 3 + 2];
                                               // *this << setw(32) << stresses[i] << std::endl;

                                                  *this << std::endl;
                                         }
                                }
                                break;

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}



// postprocess 
void COutputter2::OutputElementStress2()
{
	CDomain* FEMData = CDomain::Instance();

	double* Displacement = FEMData->GetDisplacement();

	const unsigned int NUMEG = FEMData->GetNUMEG(); // number of element group

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << "TITLE = \" 3DFEA \" "<< endl
			  << "VARIABLES = VARIABLES = \"X\",\"Y\",\"Z\", \"STRESS\"  " << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();  // NUMBER OF ELEMENT
		unsigned int NUMNP = FEMData->GetNUMNP();// NUMBER OF NODE POINTS

		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				{
					*this << "ZONE T=\"SCENE1\", N=" <<NUMNP << ",E=" << NUME << " F=FEPOINT , ET= BRICK, C= RED" << endl;

					double stress;
					double Positions[8];
					unsigned int* ELENUM= new unsigned int[NUME*2]; 
					//********************88double Positions[6];
					//loop for every element 
					for (unsigned int Ele = 0; Ele < NUME; Ele++) 
					{
						CElement& Element = EleGrp.GetElement(Ele);
						Element.ElementStress2(&stress, Displacement, Positions);

						CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
				
						*this <<Positions[0]<< setw(4) <<Positions[1]<<setw(4) <<Positions[2]<<setw(4)  << endl
							  <<Positions[3]<< setw(4) <<Positions[4]<<setw(4) <<Positions[5]<<setw(4)  << endl;
						ELENUM[Ele*2] = (unsigned int)Positions[6];
						ELENUM[Ele*2+1] = (unsigned int)Positions[7];
					}
					for (unsigned int Ele = 0; Ele < NUME; Ele++) 
					{
						*this <<ELENUM[Ele*2]<<setw(4)<<ELENUM[Ele*2+1]<<endl;
					}

					//*this << endl;
					break;
				}
			
			case ElementTypes::Quadrilateral: // Quadrilateral element
				*this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= QUADRILATERAL, C= RED" << endl;
				
				double stresses[12];
				double Position4Q[12];
				

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					
					dynamic_cast<CQuadrilateral&>(
						EleGrp.GetElement(Ele)).ElementStress2(stresses, Displacement, Position4Q);
					
					for (unsigned i=0; i<4; ++i) { // four gauss points
						
						*this << Position4Q[i*3] << setw(4) << Position4Q[i*3+1] << setw(4) << Position4Q[i*3+2]<< setw(4) ;
						
						
						*this << std::endl;
					}
				}
				*this << endl;

				break;

			case ElementTypes::Beam: // Beam element
				*this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= BRICK, C= RED" << endl;

				double beamstress[3];
				double PositionBeam[24];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress2(beamstress, Displacement, PositionBeam);

					CBeamMaterial& material =
						*dynamic_cast<CBeamMaterial*>(Element.GetElementMaterial());

					*this << PositionBeam[0] <<setw(4) << PositionBeam[1] <<setw(4) << PositionBeam[2] <<setw(4) <<endl
						  << PositionBeam[3] <<setw(4) << PositionBeam[4] <<setw(4) << PositionBeam[5] <<setw(4) <<endl
						  << PositionBeam[6] <<setw(4) << PositionBeam[7] <<setw(4) << PositionBeam[8] <<setw(4) <<endl
						  << PositionBeam[9] <<setw(4) << PositionBeam[10] <<setw(4) << PositionBeam[11] <<setw(4) <<endl
						  << PositionBeam[12] <<setw(4) << PositionBeam[13] <<setw(4) << PositionBeam[14] <<setw(4) <<endl
						  << PositionBeam[15] <<setw(4) << PositionBeam[16] <<setw(4) << PositionBeam[17] <<setw(4) <<endl
						  << PositionBeam[18] <<setw(4) << PositionBeam[19] <<setw(4) << PositionBeam[20] <<setw(4) <<endl
						  << PositionBeam[21] <<setw(4) << PositionBeam[22] <<setw(4) << PositionBeam[23] <<setw(4) <<endl;
				}

				*this << endl;
				break;	

			case ElementTypes::Triangle: // 3T element
				*this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= TRIANGLE, C= RED" << endl;
				double stress3T[3];
				double Position3T[9];
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress2(stress3T, Displacement, Position3T);
					CTriangleMaterial material = *dynamic_cast<CTriangleMaterial*>(Element.GetElementMaterial());
					
					*this << Position3T[0] <<setw(4) << Position3T[1] <<setw(4) << Position3T[2] <<setw(4) <<endl
						  << Position3T[3] <<setw(4) << Position3T[4] <<setw(4) << Position3T[5] <<setw(4) <<endl
						  << Position3T[6] <<setw(4) << Position3T[7] <<setw(4) << Position3T[8] <<setw(4) <<endl;
						
				}

				*this << endl;
				break;


			case ElementTypes::Hexahedron: // 8H element
				*this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= BRICK, C= RED" << endl;
				
				double stressHex[48];
				double Position8H[24];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress2(stressHex, Displacement, Position8H);

					CHexMaterial& material = *dynamic_cast<CHexMaterial*>(Element.GetElementMaterial());
					*this << Position8H[0] <<setw(4) << Position8H[1] <<setw(4) << Position8H[2] <<setw(4) <<endl
						  << Position8H[3] <<setw(4) << Position8H[4] <<setw(4) << Position8H[5] <<setw(4) <<endl
						  << Position8H[6] <<setw(4) << Position8H[7] <<setw(4) << Position8H[8] <<setw(4) <<endl
						  << Position8H[9] <<setw(4) << Position8H[10] <<setw(4) << Position8H[11] <<setw(4) <<endl
						  << Position8H[12] <<setw(4) << Position8H[13] <<setw(4) << Position8H[14] <<setw(4) <<endl
						  << Position8H[15] <<setw(4) << Position8H[16] <<setw(4) << Position8H[17] <<setw(4) <<endl
						  << Position8H[18] <<setw(4) << Position8H[19] <<setw(4) << Position8H[20] <<setw(4) <<endl
						  << Position8H[21] <<setw(4) << Position8H[22] <<setw(4) << Position8H[23] <<setw(4) <<endl;
				}

				*this << endl;
				break;

			case ElementTypes::TimoshenkoSRINT: // TimoshenkoSRINT beam element
				double TimoshenkoStresses[3];
				double TimoshenkoForces[12];
				double PositionTB[24];
				*this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= BRICK, C= RED" << endl;
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					dynamic_cast<CTimoshenkoSRINT&>(
						EleGrp.GetElement(Ele)).ElementStress2(TimoshenkoStresses, TimoshenkoForces, Displacement, PositionTB);

					*this << PositionTB[0] <<setw(4) << PositionTB[1] <<setw(4) << PositionTB[2] <<setw(4) <<endl
						  << PositionTB[3] <<setw(4) << PositionTB[4] <<setw(4) << PositionTB[5] <<setw(4) <<endl
						  << PositionTB[6] <<setw(4) << PositionTB[7] <<setw(4) << PositionTB[8] <<setw(4) <<endl
						  << PositionTB[9] <<setw(4) << PositionTB[10] <<setw(4) << PositionTB[11] <<setw(4) <<endl
						  << PositionTB[12] <<setw(4) << PositionTB[13] <<setw(4) << PositionTB[14] <<setw(4) <<endl
						  << PositionTB[15] <<setw(4) << PositionTB[16] <<setw(4) << PositionTB[17] <<setw(4) <<endl
						  << PositionTB[18] <<setw(4) << PositionTB[19] <<setw(4) << PositionTB[20] <<setw(4) <<endl
						  << PositionTB[21] <<setw(4) << PositionTB[22] <<setw(4) << PositionTB[23] <<setw(4) <<endl;
					*this << std::endl;
				}
				break;

			case ElementTypes::TimoshenkoEBMOD: // TimoshenkoEBMOD beam element
				double TimoshenkoEBStresses[3];
				double TimoshenkoEBForces[12];
				double PositionTEB[24];
				*this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= BRICK, C= RED" << endl;
				
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					dynamic_cast<CTimoshenkoEBMOD&>(
						EleGrp.GetElement(Ele)).ElementStress2(TimoshenkoEBStresses, TimoshenkoEBForces, Displacement, PositionTEB);

					*this << PositionTEB[0] <<setw(4) << PositionTEB[1] <<setw(4) << PositionTEB[2] <<setw(4) <<endl
						  << PositionTEB[3] <<setw(4) << PositionTEB[4] <<setw(4) << PositionTEB[5] <<setw(4) <<endl
						  << PositionTEB[6] <<setw(4) << PositionTEB[7] <<setw(4) << PositionTEB[8] <<setw(4) <<endl
						  << PositionTEB[9] <<setw(4) << PositionTEB[10] <<setw(4) << PositionTEB[11] <<setw(4) <<endl
						  << PositionTEB[12] <<setw(4) << PositionTEB[13] <<setw(4) << PositionTEB[14] <<setw(4) <<endl
						  << PositionTEB[15] <<setw(4) << PositionTEB[16] <<setw(4) << PositionTEB[17] <<setw(4) <<endl
						  << PositionTEB[18] <<setw(4) << PositionTEB[19] <<setw(4) << PositionTEB[20] <<setw(4) <<endl
						  << PositionTEB[21] <<setw(4) << PositionTEB[22] <<setw(4) << PositionTEB[23] <<setw(4) <<endl;
					*this << std::endl;
				}
				break;

			case ElementTypes::Plate:
				*this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= QUADRILATERAL, C= RED" << endl;
				
					
				double stresses4PE[12];
				double Positions4PE[12];
				

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					dynamic_cast<CPlate&>(
						EleGrp.GetElement(Ele)).ElementStress2(stresses4PE, Displacement, Positions4PE);
					
					for (unsigned i=0; i<4; ++i) { // four gauss points
						*this << Positions4PE[0] <<setw(4) << Positions4PE[1] <<setw(4) << Positions4PE[2] <<setw(4) <<endl
						  << Positions4PE[3] <<setw(4) << Positions4PE[4] <<setw(4) << Positions4PE[5] <<setw(4) <<endl
						  << Positions4PE[6] <<setw(4) << Positions4PE[7] <<setw(4) << Positions4PE[8] <<setw(4) <<endl
						  << Positions4PE[9] <<setw(4) << Positions4PE[10] <<setw(4) << Positions4PE[11] <<setw(4) <<endl;
						
						*this << std::endl;
					}
				}
				*this << endl;

				break;


                case ElementTypes::Shell:
                 *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME << " F=FEPOINT , ET= QUADRILATERAL, C= RED" << endl;               
                                double stresses4SE[15];
                                double Positions4SE[15];
                                for (unsigned int Ele = 0; Ele < NUME; Ele++)
                                {
                                        dynamic_cast<CShell&>(EleGrp.GetElement(Ele))
                                                .ElementStress2(stresses4SE, Displacement, Positions4SE);

                                        for (unsigned i = 0; i < 5; ++i){
                                          // four gauss points;
                                          //THE FIFTH POINT IS THE CENTRE POINT FOR IN-PLANE STRESSES
                                                  *this << Positions4SE[0] <<setw(4) << Positions4SE[1] <<setw(4) << Positions4SE[2] <<setw(4) <<endl
														<< Positions4SE[3] <<setw(4) << Positions4SE[4] <<setw(4) << Positions4SE[5] <<setw(4) <<endl
														<< Positions4SE[6] <<setw(4) << Positions4SE[7] <<setw(4) << Positions4SE[8] <<setw(4) <<endl
														<< Positions4SE[9] <<setw(4) << Positions4SE[10] <<setw(4) << Positions4SE[11] <<setw(4) <<endl;
                                               // *this << setw(32) << stresses[i] << std::endl;

                                                  *this << std::endl;
                                         }
                                }
                                break;

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetNWK()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetMK()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = "
		  << FEMData->GetNWK() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}



#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	unsigned int* ColumnHeights = FEMData->GetStiffnessMatrix().GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
 		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	unsigned int* DiagonalAddress = FEMData->GetStiffnessMatrix().GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
#ifdef MKL
	*this << "*** _Debug_ *** CSR stiffness matrix" << std::endl;
	*this << CDomain::Instance()->GetCSRStiffnessMatrix() << std::endl;
#else
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double>& StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix.GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << StiffnessMatrix(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (unsigned I = 1; I <= NEQ; I++)
	{
		for (unsigned J = 1; J <= NEQ; J++)
		{
            unsigned i, j;
            i = std::min(I, J);
            j = std::max(I, J);
			unsigned H = DiagonalAddress[j] - DiagonalAddress[j - 1];
			if (j - i - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << StiffnessMatrix(i, j);
			}
		}

		*this << endl;
	}
#endif
	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement(unsigned int loadcase)
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << "  Load case = " << loadcase << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
