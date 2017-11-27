/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/Triangle.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//  Constructor
CTriangle::CTriangle()
{
    NEN = 3; // Each element has 3 nodes
    nodes = new CNode*[NEN];

    ND = 9;
    LocationMatrix = new unsigned int[ND];

    ElementMaterial = NULL;
}

//  Desconstructor
CTriangle::~CTriangle() {}

//  Read element data from stream Input
bool CTriangle::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N;

    Input >> N; // element number

    if (N != Ele + 1)
    {
        cerr << "*** Error *** Elements must be inputted in order !" << endl
             << "    Expected element : " << Ele + 1 << endl
             << "    Provided element : " << N << endl;

        return false;
    }

    unsigned int MSet;       // Material property set number
    unsigned int N1, N2, N3; // node number

    Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial = dynamic_cast<CTriangleMaterial*>(MaterialSets) + MSet - 1;
    nodes[0] = NodeList + N1 - 1;
    nodes[1] = NodeList + N2 - 1;
    nodes[2] = NodeList + N3 - 1;

    return true;
}

//  Write element data to stream
void CTriangle::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele + 1 << setw(11) << nodes[0]->NodeNumber << setw(9)
           << nodes[1]->NodeNumber << setw(12) << nodes[2]->NodeNumber << setw(12)
           << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//  Caution:  Equation number is numbered from 1 !
void CTriangle::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//  Return the size of the element stiffness matrix (stored as an array column by column)
unsigned int CTriangle::SizeOfStiffnessMatrix()
{
    return 45;
    // 9 * 10 / 2
}

inline double dot(const double* p1, const double* p2)
{
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

void Convert2d23d3T(const double* ke, double* Matrix, const double i[3], const double j[3])
{
    Matrix[0] = i[0] * (i[0] * ke[0] + j[0] * ke[2]) + j[0] * j[0] * ke[1];
    Matrix[1] = i[1] * (i[1] * ke[0] + j[1] * ke[2]) + j[1] * j[1] * ke[1];
    Matrix[2] = j[0] * j[1] * ke[1] + i[0] * (i[1] * ke[0] + j[1] * ke[2]);
    Matrix[3] = i[2] * i[2] * ke[0] + j[2] * j[2] * ke[1] + i[2] * j[2] * ke[2];
    Matrix[4] = j[1] * j[2] * ke[1] + i[1] * (i[2] * ke[0] + j[2] * ke[2]);
    Matrix[5] = j[0] * j[2] * ke[1] + i[0] * (i[2] * ke[0] + j[2] * ke[2]);
    Matrix[6] = i[0] * i[0] * ke[3] + j[0] * j[0] * ke[6] + i[0] * j[0] * ke[7];
    Matrix[7] = i[0] * (j[2] * ke[4] + i[2] * ke[5]) + j[0] * (j[2] * ke[8] + i[2] * ke[9]);
    Matrix[8] = i[0] * (j[1] * ke[4] + i[1] * ke[5]) + j[0] * (j[1] * ke[8] + i[1] * ke[9]);
    Matrix[9] = i[0] * i[0] * ke[5] + j[0] * j[0] * ke[8] + i[0] * j[0] * (ke[4] + ke[9]);
    Matrix[10] = i[1] * i[1] * ke[3] + j[1] * j[1] * ke[6] + i[1] * j[1] * ke[7];
    Matrix[11] = j[0] * j[1] * ke[6] + i[0] * (i[1] * ke[3] + j[1] * ke[7]);
    Matrix[12] = i[1] * (j[2] * ke[4] + i[2] * ke[5]) + j[1] * (j[2] * ke[8] + i[2] * ke[9]);
    Matrix[13] = i[1] * i[1] * ke[5] + j[1] * j[1] * ke[8] + i[1] * j[1] * (ke[4] + ke[9]);
    Matrix[14] = j[0] * (i[1] * ke[4] + j[1] * ke[8]) + i[0] * (i[1] * ke[5] + j[1] * ke[9]);
    Matrix[15] = i[2] * i[2] * ke[3] + j[2] * j[2] * ke[6] + i[2] * j[2] * ke[7];
    Matrix[16] = j[1] * j[2] * ke[6] + i[1] * (i[2] * ke[3] + j[2] * ke[7]);
    Matrix[17] = j[0] * j[2] * ke[6] + i[0] * (i[2] * ke[3] + j[2] * ke[7]);
    Matrix[18] = i[2] * i[2] * ke[5] + j[2] * j[2] * ke[8] + i[2] * j[2] * (ke[4] + ke[9]);
    Matrix[19] = j[1] * (i[2] * ke[4] + j[2] * ke[8]) + i[1] * (i[2] * ke[5] + j[2] * ke[9]);
    Matrix[20] = j[0] * (i[2] * ke[4] + j[2] * ke[8]) + i[0] * (i[2] * ke[5] + j[2] * ke[9]);
    Matrix[21] = i[0] * i[0] * ke[10] + j[0] * j[0] * ke[15] + i[0] * j[0] * ke[16];
    Matrix[22] = i[0] * (j[2] * ke[11] + i[2] * ke[12]) + j[0] * (j[2] * ke[17] + i[2] * ke[18]);
    Matrix[23] = i[0] * (j[1] * ke[11] + i[1] * ke[12]) + j[0] * (j[1] * ke[17] + i[1] * ke[18]);
    Matrix[24] = i[0] * i[0] * ke[12] + j[0] * j[0] * ke[17] + i[0] * j[0] * (ke[11] + ke[18]);
    Matrix[25] = i[0] * (j[2] * ke[13] + i[2] * ke[14]) + j[0] * (j[2] * ke[19] + i[2] * ke[20]);
    Matrix[26] = i[0] * (j[1] * ke[13] + i[1] * ke[14]) + j[0] * (j[1] * ke[19] + i[1] * ke[20]);
    Matrix[27] = i[0] * i[0] * ke[14] + j[0] * j[0] * ke[19] + i[0] * j[0] * (ke[13] + ke[20]);
    Matrix[28] = i[1] * i[1] * ke[10] + j[1] * j[1] * ke[15] + i[1] * j[1] * ke[16];
    Matrix[29] = j[0] * j[1] * ke[15] + i[0] * (i[1] * ke[10] + j[1] * ke[16]);
    Matrix[30] = i[1] * (j[2] * ke[11] + i[2] * ke[12]) + j[1] * (j[2] * ke[17] + i[2] * ke[18]);
    Matrix[31] = i[1] * i[1] * ke[12] + j[1] * j[1] * ke[17] + i[1] * j[1] * (ke[11] + ke[18]);
    Matrix[32] = j[0] * (i[1] * ke[11] + j[1] * ke[17]) + i[0] * (i[1] * ke[12] + j[1] * ke[18]);
    Matrix[33] = i[1] * (j[2] * ke[13] + i[2] * ke[14]) + j[1] * (j[2] * ke[19] + i[2] * ke[20]);
    Matrix[34] = i[1] * i[1] * ke[14] + j[1] * j[1] * ke[19] + i[1] * j[1] * (ke[13] + ke[20]);
    Matrix[35] = j[0] * (i[1] * ke[13] + j[1] * ke[19]) + i[0] * (i[1] * ke[14] + j[1] * ke[20]);
    Matrix[36] = i[2] * i[2] * ke[10] + j[2] * j[2] * ke[15] + i[2] * j[2] * ke[16];
    Matrix[37] = j[1] * j[2] * ke[15] + i[1] * (i[2] * ke[10] + j[2] * ke[16]);
    Matrix[38] = j[0] * j[2] * ke[15] + i[0] * (i[2] * ke[10] + j[2] * ke[16]);
    Matrix[39] = i[2] * i[2] * ke[12] + j[2] * j[2] * ke[17] + i[2] * j[2] * (ke[11] + ke[18]);
    Matrix[40] = j[1] * (i[2] * ke[11] + j[2] * ke[17]) + i[1] * (i[2] * ke[12] + j[2] * ke[18]);
    Matrix[41] = j[0] * (i[2] * ke[11] + j[2] * ke[17]) + i[0] * (i[2] * ke[12] + j[2] * ke[18]);
    Matrix[42] = i[2] * i[2] * ke[14] + j[2] * j[2] * ke[19] + i[2] * j[2] * (ke[13] + ke[20]);
    Matrix[43] = j[1] * (i[2] * ke[13] + j[2] * ke[19]) + i[1] * (i[2] * ke[14] + j[2] * ke[20]);
    Matrix[44] = j[0] * (i[2] * ke[13] + j[2] * ke[19]) + i[0] * (i[2] * ke[14] + j[2] * ke[20]);
}

//  Calculate element stiffness matrix
//  Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CTriangle::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];

    // make p31 p21
    double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};
    double p32[3] = {n3.XYZ[0] - n2.XYZ[0], n3.XYZ[1] - n2.XYZ[1], n3.XYZ[2] - n2.XYZ[2]};

    // n = p31 cross p21 (normalized)
    double n[3] = {p31[1] * p21[2] - p31[2] * p21[1], p31[2] * p21[0] - p31[0] * p21[2],
                   p31[0] * p21[1] - p31[1] * p21[0]};
    // generate area and normalize n at the same time
    double Area = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= Area;
    n[1] /= Area;
    n[2] /= Area;
    Area /= 2.0;
    // i = normalized p21, so that y21 = 0
    double i[3] = {p21[0], p21[1], p21[1]};
    normalize(i);
    // j = n cross i
    double const j[3] = {n[1] * i[2] - n[2] * i[1], n[2] * i[0] - n[0] * i[2],
                         n[0] * i[1] - n[1] * i[0]};

    // generate M here
    double x32 = dot(p32, i);
    double y23 = -dot(p32, j);
    double x13 = -dot(p31, i);
    double y31 = dot(p31, j);
    double x21 = dot(p21, i);

    CTriangleMaterial& material =
        *dynamic_cast<CTriangleMaterial*>(ElementMaterial); // Pointer to material of the element

    auto E = material.E;
    auto v = material.nu;

    double ke[21];
    double cof = E / (8 * Area * (1 - v * v));
    double va = 1 - v;

    // ke is stored down to up
    ke[0] = cof * (va * (x32 * x32) + 2 * (y23 * y23));
    ke[1] = cof * (2 * (x32 * x32) + va * (y23 * y23));
    ke[2] = cof * ((1 + v) * x32 * y23);
    ke[3] = cof * (va * (x13 * x13) + 2 * (y31 * y31));
    ke[4] = cof * (x13 * (y23 - v * y23) + 2 * v * x32 * y31);
    ke[5] = cof * (x13 * (x32 - v * x32) + 2 * y23 * y31);
    ke[6] = cof * (2 * (x13 * x13) + va * (y31 * y31));
    ke[7] = cof * ((1 + v) * x13 * y31);
    ke[8] = cof * (2 * x13 * x32 + va * y23 * y31);
    ke[9] = cof * (2 * v * x13 * y23 + va * x32 * y31);
    ke[10] = cof * (va * (x21 * x21));
    ke[11] = cof * (va * x21 * y31);
    ke[12] = cof * (va * x13 * x21);
    ke[13] = cof * (va * x21 * y23);
    ke[14] = cof * (va * x21 * x32);
    ke[15] = cof * (2 * (x21 * x21));
    ke[16] = 0;
    ke[17] = cof * (2 * x13 * x21);
    ke[18] = cof * (2 * v * x21 * y31);
    ke[19] = cof * (2 * x21 * x32);
    ke[20] = cof * (2 * v * x21 * y23);

    Convert2d23d3T(ke, Matrix, i, j);
}

//  Calculate element stress
void CTriangle::ElementStress(double stress[3], double* Displacement)
{
    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];

    // make p31 p21
    double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};
    double p32[3] = {n3.XYZ[0] - n2.XYZ[0], n3.XYZ[1] - n2.XYZ[1], n3.XYZ[2] - n2.XYZ[2]};

    // n = p31 cross p21 (normalized)
    double n[3] = {p31[1] * p21[2] - p31[2] * p21[1], p31[2] * p21[0] - p31[0] * p21[2],
                   p31[0] * p21[1] - p31[1] * p21[0]};
    // generate area and normalize n at the same time
    double Area = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= Area;
    n[1] /= Area;
    n[2] /= Area;
    Area /= 2.0;
    // i = normalized p21, so that y21 = 0
    double i[3] = {p21[0], p21[1], p21[1]};
    normalize(i);
    // j = n cross i
    double const j[3] = {n[1] * i[2] - n[2] * i[1], n[2] * i[0] - n[0] * i[2],
                         n[0] * i[1] - n[1] * i[0]};

    // generate M here
    double x32 = dot(p32, i);
    double y23 = -dot(p32, j);
    double x13 = -dot(p31, i);
    double y31 = dot(p31, j);
    double x21 = dot(p21, i);

    // form d first.
    // d represent 3d displacements at boundary nodes.
    double d[9];
    for (unsigned index = 0; index < 9; ++index)
    {
        if (LocationMatrix[index])
            d[index] = Displacement[LocationMatrix[index] - 1];
        else
            d[index] = 0;
    }

    double de[6] = {
        i[0] * d[0] + i[1] * d[3] + i[2] * d[6],
        j[0] * d[0] + j[1] * d[3] + j[2] * d[6], // node 1 (2d)
        i[0] * d[1] + i[1] * d[4] + i[2] * d[7],
        j[0] * d[1] + j[1] * d[4] + j[2] * d[7], // node 2 (2d)
        i[0] * d[2] + i[1] * d[5] + i[2] * d[8],
        j[0] * d[2] + j[1] * d[5] + j[2] * d[8] // node 3 (2d)
    };

    CTriangleMaterial& material = *dynamic_cast<CTriangleMaterial*>(ElementMaterial);

    double v = material.nu;
    double cof = material.E / (4 * Area * (1 - v * v));

    stress[0] =
        2 * cof * (y23 * de[0] + y31 * de[2] + v * (x32 * de[1] + x13 * de[3] + x21 * de[5]));
    stress[1] = 2 * (v * y23 * de[0] + x32 * de[1] + v * y31 * de[2] + x13 * de[3] + x21 * de[5]);
    stress[2] = (1 - v) * (x32 * de[0] + y23 * de[1] + x13 * de[2] + y31 * de[3] + x21 * de[4]);
}
