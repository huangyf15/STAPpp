/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Domain.h"

CNode* CElementGroup::NodeList_ = nullptr;

//! Constructor
CElementGroup::CElementGroup()
{
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::Instance();
        NodeList_ = FEMData->GetNodeList();
    }

    ElementType_ = ElementTypes::UNDEFINED;

    NUME_ = 0;
    ElementList_ = nullptr;

    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! Deconstructor
CElementGroup::~CElementGroup()
{
    if (ElementList_)
        delete[] ElementList_;

    if (MaterialList_)
        delete[] MaterialList_;
}

CElement& CElementGroup::GetElement(unsigned int index)
{
    return *(CElement*)((std::size_t)(ElementList_) + index * ElementSize_);
}

CMaterial& CElementGroup::GetMaterial(unsigned int index)
{
    return *(CMaterial*)((std::size_t)(MaterialList_) + index * MaterialSize_);
}

void CElementGroup::CalculateMemberSize()
{
    switch (ElementType_)
    {
        case ElementTypes::UNDEFINED:
            std::cerr << "Setting element type to UNDEFINED." << std::endl;
            exit(5);
        case ElementTypes::Bar:
            ElementSize_ = sizeof(CBar);
            MaterialSize_ = sizeof(CBarMaterial);
            break;
        case ElementTypes::Quadrilateral:
            ElementSize_ = sizeof(CQuadrilateral);
            MaterialSize_ = sizeof(CQuadrilateralMaterial);
            break;
        case ElementTypes::Triangle:
            ElementSize_ = sizeof(CTriangle);
            MaterialSize_ = sizeof(CTriangleMaterial);
            break;
        case ElementTypes::Hexahedron:
            ElementSize_ = sizeof(CHex);
            MaterialSize_ = sizeof(CHexMaterial);
            break;
        case ElementTypes::Beam:
            ElementSize_ = sizeof(CBeam);
            MaterialSize_ = sizeof(CBeamMaterial);
            break;
        case ElementTypes::TimoshenkoSRINT:
            ElementSize_ = sizeof(CTimoshenkoSRINT);
            MaterialSize_ = sizeof(CTimoshenkoMaterial);
            break;
        case ElementTypes::TimoshenkoEBMOD:
            ElementSize_ = sizeof(CTimoshenkoEBMOD);
            MaterialSize_ = sizeof(CTimoshenkoMaterial);
            break;
        case ElementTypes::Plate:
            ElementSize_ = sizeof(CPlate);
            MaterialSize_ = sizeof(CPlateMaterial);
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not finished yet. See CElementGroup::CalculateMemberSize." << std::endl;
            exit(5);
            break;
    }
}

void CElementGroup::AllocateElement(std::size_t size)
{
    switch (ElementType_)
    {
        case ElementTypes::Bar:
            ElementList_ = new CBar[size];
            break;
        case ElementTypes::Quadrilateral:
            ElementList_ = new CQuadrilateral[size];
            break; 
        case ElementTypes::Triangle:
            ElementList_ = new CTriangle[size];
            break;
        case ElementTypes::Hexahedron:
            ElementList_ = new CHex[size];
            break;
        case ElementTypes::Beam:
            ElementList_ = new CBeam[size];
            break;
        case ElementTypes::TimoshenkoSRINT:
            ElementList_ = new CTimoshenkoSRINT[size];
            break;
        case ElementTypes::TimoshenkoEBMOD:
            ElementList_ = new CTimoshenkoEBMOD[size];
            break;
        case ElementTypes::Plate:
            ElementList_ = new CPlate[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not finished yet. See CElementGroup::AllocateElement." << std::endl;
            exit(5);
    }
}

void CElementGroup::AllocateMaterial(std::size_t size)
{
    switch (ElementType_)
    {
        case ElementTypes::Bar:
            MaterialList_ = new CBarMaterial[size];
            break;
        case ElementTypes::Quadrilateral:
            MaterialList_ = new CQuadrilateralMaterial[size];
            break;
        case ElementTypes::Triangle:
            MaterialList_ = new CTriangleMaterial[size];
            break;
        case ElementTypes::Hexahedron:
            MaterialList_ = new CHexMaterial[size];
            break;
        case ElementTypes::Beam:
            MaterialList_ = new CBeamMaterial[size];
            break;
        case ElementTypes::TimoshenkoSRINT:
            MaterialList_ = new CTimoshenkoMaterial[size];
            break;
        case ElementTypes::TimoshenkoEBMOD:
            MaterialList_ = new CTimoshenkoMaterial[size];
            break;
        case ElementTypes::Plate:
            MaterialList_ = new CPlateMaterial[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not finished yet. See CElementGroup::AllocateMaterial." << std::endl;
            exit(5);
    }
}

//! Read element group data from stream Input
bool CElementGroup::Read(ifstream& Input)
{
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;

    CalculateMemberSize();

    if (!ReadElementData(Input))
        return false;

    return true;
}

//  Read bar element data from the input data file
bool CElementGroup::ReadElementData(ifstream& Input)
{
    //  Read material/section property lines
    AllocateMaterial(NUMMAT_);

    //  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
        if (!GetMaterial(mset).Read(Input, mset))
            return false;

    //  Read element data lines
    AllocateElement(NUME_);

    //  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
        if (!GetElement(Ele).Read(Input, Ele, MaterialList_, NodeList_))
            return false;

    return true;
}

//  Read quadrilateral element data from the input data file
bool CElementGroup::ReadQuadrilateralElementData(ifstream& Input)
{
    //  Read material/section property lines
    MaterialList_ = new CQuadrilateralMaterial[NUMMAT_]; // Materials for group EleGrp
    CQuadrilateralMaterial* mlist = dynamic_cast<CQuadrilateralMaterial*>(MaterialList_);

    //  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
        if (!mlist[mset].Read(Input, mset))
            return false;

    //  Read element data lines
    ElementList_ = new CQuadrilateral[NUME_]; // Elements of group EleGrp
    CQuadrilateral* elist = dynamic_cast<CQuadrilateral*>(ElementList_);

    //  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
        if (!elist[Ele].Read(Input, Ele, MaterialList_, NodeList_))
            return false;

    return true;
}

bool CElementGroup::ReadPlateElementData(ifstream& Input)
{
    //  Read material/section property lines
    MaterialList_ = new CPlateMaterial[NUMMAT_]; // Materials for group EleGrp
    CPlateMaterial* mlist = dynamic_cast<CPlateMaterial*>(MaterialList_);

    //  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
        if (!mlist[mset].Read(Input, mset))
            return false;

    //  Read element data lines
    ElementList_ = new CPlate[NUME_]; // Elements of group EleGrp
    CPlate* elist = dynamic_cast<CPlate*>(ElementList_);

    //  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
        if (!elist[Ele].Read(Input, Ele, MaterialList_, NodeList_))
            return false;

    return true;
}
