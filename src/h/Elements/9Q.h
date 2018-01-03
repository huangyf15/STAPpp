
#pragma once

#include "Element.h"
#include <cmath>

using namespace std;

//! Bar element class
class C9Q : public CElement
{
public:
    //!	Constructor
    C9Q();

    //!	Desconstructor
    ~C9Q();

    //!	Read element data from stream Input
    virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

    //!	Write element data to stream
    virtual void Write(COutputter& output, unsigned int Ele);

    //! Generate location matrix: the global equation number that corresponding to each DOF of the
    //! element
    //	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

    //!	Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //!	Calculate element stress
    virtual void ElementStress(double* stress, double* Displacement, double* Positions = nullptr);

    //!	Return the size of the element stiffness matrix (stored as an array column by column)
    virtual unsigned int SizeOfStiffnessMatrix();

#ifdef _VIB_
//!	Calculate element mass matrix (Upper triangular matrix, stored as an array column by colum)
	virtual void ElementMass(double* mass); 
#endif

};
