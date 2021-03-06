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

#include "Element.h"
#include "Elements/8H.h"
#include "Elements/9Q.h"
#include "Elements/Bar.h"
#include "Elements/Beam.h"
#include "Elements/Frustum.h"
#include "Elements/Plate.h"
#include "Elements/Quadrilateral.h"
#include "Elements/Shell.h"
#include "Elements/TimoshenkoEBMOD.h"
#include "Elements/TimoshenkoSRINT.h"
#include "Elements/Triangle.h"
#include "Elements/Infinite_4Q.h"
#include "Elements/5Q.h"
#include "Material.h"
#include "Node.h"

using namespace std;

enum ElementTypes
{
    UNDEFINED = 0,
    Bar = 1,
    Quadrilateral = 2,
    Triangle = 3,
    Hexahedron = 4,
    Beam = 5,
    Plate = 6,
    Shell = 7,
    TimoshenkoSRINT = 8,
    TimoshenkoEBMOD = 9,
    T9Q = 10,
    Infinite = 11,
    T5Q = 12,
    Frustum = 13
};

//! Element group class
class CElementGroup
{
private:
    //! List of all nodes in the domain, obtained from CDomain object
    static CNode* NodeList_;

    //! Element type of this group
    ElementTypes ElementType_;

    //! Element size of this group
    std::size_t ElementSize_;

    std::size_t MaterialSize_;

    //! Number of elements in this group
    unsigned int NUME_;

    //! Element List in this group
    CElement* ElementList_;

    //! Number of material/section property sets in this group
    unsigned int NUMMAT_;

    //! Material list in this group
    CMaterial* MaterialList_;

public:
    //! Constructor
    CElementGroup();

    //! Destructor
    ~CElementGroup();

    //! Read element group data from stream Input
    bool Read(ifstream& Input);

    void CalculateMemberSize();

    void AllocateElement(std::size_t size);

    void AllocateMaterial(std::size_t size);

    //! Read element data from the input data file
    bool ReadElementData(ifstream& Input);

    //! Return element type of this group
    ElementTypes GetElementType() { return ElementType_; }

    //! Return the number of elements in the group
    unsigned int GetNUME() { return NUME_; }

    CElement& GetElement(unsigned int index);

    CMaterial& GetMaterial(unsigned int index);

    //! Return the number of material/section property setss in this element group
    unsigned int GetNUMMAT() { return NUMMAT_; }
};
