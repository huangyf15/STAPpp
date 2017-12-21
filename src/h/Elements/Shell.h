#pragma once

#include "Element.h"

using namespace std;

//! Shell element class
class CShell : public CElement
{
public:
    //!	Constructor
    CShell();

    //!	Desconstructor
    ~CShell();

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
    virtual void ElementStress(double* stress, double* Displacement,double* position);

    //!	Calculate element stress
    virtual void ElementStress2(double* stress, double* Displacement, double* PrePositions, double* PostPositions);

    //!	Return the size of the element stiffness matrix (stored as an array column by column)
    virtual unsigned int SizeOfStiffnessMatrix();
};
