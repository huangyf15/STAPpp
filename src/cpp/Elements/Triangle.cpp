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
    NEN = 3; // Each element has 2 nodes
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
    nodes[0] = &NodeList[N1 - 1];
    nodes[1] = &NodeList[N2 - 1];
    nodes[2] = &NodeList[N3 - 1];

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
//  For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//  has 21 elements
unsigned int CTriangle::SizeOfStiffnessMatrix() { return 45; }

inline void normalize(double ptr[3])
{
    double sum = std::sqrt((ptr[0] * ptr[0]) + (ptr[1] * ptr[1]) + (ptr[2] + ptr[2]));
    ptr[0] /= sum;
    ptr[1] /= sum;
    ptr[2] /= sum;
}

inline double dot(double* p1, double* p2) { return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]; }

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
    double const n[3] = {p31[1] * p21[2] - p31[2] * p21[1], p31[2] * p21[0] - p31[0] * p21[2],
                         p31[0] * p21[1] - p31[1] * p21[0]};
    // generate area and normalize n at the same time
    double Aera = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= Aera;
    n[1] /= Aera;
    n[2] /= Aera;
    Aera /= 2.0;
    // i = normalized p21
    // i is manually set parallel to p21 so that y21 = 0
    double const i[3] = {p21[0], p21[1], p21[1]};
    normalize(i);
    // j = n cross i
    double const j[3] = {n[1] * i[2] - n[2] * i[1], n[2] * i[0] - n[0] * i[2],
                         n[0] * i[1] - n[1] * i[0]};

    // by here, a conversion matrix is formed,
    // as (x', y') = ((i0, i1, i2), (j0, j1, j2)) . (x, y, z)

    // form R = {
    //  i0, i1, i2, 0,0,0, 0,0,0
    //  j0, j1, j2, 0,0,0, 0,0,0
    //  0,  0,  0,  [ i ],   0
    //      0       [ j ],   0
    //   ....
    // }

    //  consider Ke = R^T ke' R = R^T (A B^T D B) R, B = BB / (2*A)
    // let M = (BB * R), then Ke = (1/4A) (M^T D M)
    // to see how M is generated, see ../../memo/tri.nb, a mathematica file.

    // generate M here
    double x32 = dot(p32, i);
    double y23 = -dot(p32, j) double x13 = -dot(p31, i);
    double y31 = dot(p31, j);
    double x21 = dot(p21, i);
    // double y12 = -dot(p21, j); but notice y12 = 0
    double M[27] =
        {
            i[0] * y23, i[1] * y23, i[2] * y23, i[0] * y31, i[1] * y31, i[2] * y31, 0, 0, 0,
            // first row
            j[0] * x32, j[1] * x32, j[2] * x32, j[0] * x13, j[1] * x13, j[2] * x13, j[0] * x21,
            j[1] * x21, j[2] * x21,
            // second row
            i[0] * x32 + j[0] * y23, i[1] * x32 + j[1] * y23, i[2] * x32 + j[2] * y23,
            i[0] * x13 + j[0] * y31, i[1] * x13 + j[1] * y31, i[2] * x13 + j[2] * y31, i[0] * x21,
            i[1] * x21, i[2]* x21
            // last row
        }

    CTriangleMaterial* material =
        dynamic_cast<CTriangleMaterial*>(ElementMaterial); // Pointer to material of the element

    const double& E = material->E;
    const double& v = material->nu;
    const double d_33 = (1 - v) / 2.0;

    double ke[21];
#ifdef m
#error "macro m is predefined"
#else
#define m(ii, jj) (M[9 * (ii - 1) + (jj - 1)])
    for (unsigned int i = 1; i <= 9; ++i)
    {
        for (unsigned int j = 1; j <= i; ++j)
        {
            // k_ij = m_ki d_kl m_lj
            ke[i - 1 + (j * (j - 1) / 2)] = (
                m(1, i) * (m(1, j) + v * m(2, j)) + 
                m(2, i) * (m(2, j) + v * m(1, j)) + 
                m(3, i) * m(3, j);
        }
    }
#undef m
}

//  Calculate element stress
void CTriangle::ElementStress(double* stress, double* Displacement)
{
    // CTriangleMaterial* material =
    //     dynamic_cast<CTriangleMaterial*>(ElementMaterial); // Pointer to material of the element

    // double DX[3];  //  dx = x2-x1, dy = y2-y1, dz = z2-z1
    // double L2 = 0; //  Square of bar length (L^2)

    // for (unsigned int i = 0; i < 3; i++)
    // {
    //     DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
    //     L2 = L2 + DX[i] * DX[i];
    // }

    // double S[6];
    // for (unsigned int i = 0; i < 3; i++)
    // {
    //     S[i] = -DX[i] * material->E / L2;
    //     S[i + 3] = -S[i];
    // }

    // *stress = 0.0;
    // for (unsigned int i = 0; i < 6; i++)
    // {
    //     if (LocationMatrix[i])
    //         *stress += S[i] * Displacement[LocationMatrix[i] - 1];
    // }
}
