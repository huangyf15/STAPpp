/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

//! Outputer class is used to output results
class COutputter
{
private:

//!	File stream for output
	ofstream OutputFile;

protected:

//!	Constructor
	COutputter(string FileName);

//!	Designed as a single instance class
	static COutputter* _instance;

public:

//!	Return pointer to the output file stream
	inline ofstream* GetOutputFile() { return &OutputFile; }

//!	Return the single instance of the class
	static COutputter* Instance(string FileName = " ");

//!	Output current time and date
	void PrintTime(const struct tm * ptm, COutputter& output);

//!	Output logo and heading 
	void OutputHeading();

//!	Output nodal point data
	void OutputNodeInfo();

//!	Output equation numbers
	void OutputEquationNumber();

//!	Output element data
	void OutputElementInfo();

//!	Output bar element data
	void PrintBarElementData(unsigned int EleGrp);

//!	Output Quadrilateral element data
	void PrintQuadrilateralElementData(unsigned int EleGrp);

//!	Output Triangle element data
	void PrintTriangleElementData(unsigned int EleGrp);

//!	Output Quadrilateral element data
	void PrintHexElementData(unsigned int EleGrp);

//!	Output bar element data
	void PrintBeamElementData(unsigned int EleGrp);

//!	Output TimoshenkoSRINT Beam element data
	void PrintTimoshenkoSRINTElementData(unsigned int EleGrp);

// !Output TimoshenkoEBMOD Beam element data
	void PrintTimoshenkoEBMODElementData(unsigned int EleGrp);

//!	Output Plate element data
    void PrintPlateElementData(unsigned int EleGrp);

//!	Output Shell element data
    void PrintShellElementData(unsigned int EleGrp);

//!	Output 9Q element data
	void Print9QElementData(unsigned int EleGrp);

//!	Output Frustum Shell element data
	void PrintFrustumElementData(unsigned int EleGrp);

//!	Output Infinite element data
	void PrintInfiniteElementData(unsigned int EleGrp);

//!	Output 5Q element data
	void Print5QElementData(unsigned int EleGrp);

//!	Output load data 
	void OutputLoadInfo(); 

//!	Output displacement data
	void OutputNodalDisplacement(unsigned int lcase);

//!	Output element stresses 
	void OutputElementStress();

//!	Print total system data
	void OutputTotalSystemData();

//! Overload the operator <<
	template <typename T>
	COutputter& operator<<(const T& item) 
	{
		#ifndef _RUN_
		std::cout << item;
		#endif
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	COutputter& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		#ifndef _RUN_
		op(std::cout);
		#endif
		op(OutputFile);
		return *this;
	}

#ifdef _DEBUG_

//!	Print banded and full stiffness matrix for debuging
	void PrintStiffnessMatrix();
#ifdef _VIB_
//! Print matrix for mass
	void PrintMassMatrix();
#endif
//!	Print address of diagonal elements for debuging
	void PrintDiagonalAddress();

//!	Print column heights for debuging
	void PrintColumnHeights();

//!	Print displacement vector for debuging
	void PrintDisplacement(unsigned int loadcase);

#endif

#ifdef _VIB_
	void PrintVibModNum();

	void OutputVibDisps();
#endif


};
