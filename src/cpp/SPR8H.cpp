#include "Elements/8H.h"
#include "Eigen/dense"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;
using namespace std;


void  CHex::ElementPostSPR(double* stressG, double* Displacement , double* PrePositions, double* PostPositions, double* PositionG, unsigned int* mingtianzaishuo)
{
	
	//construct ideal displacement vector
		//LM can be used here

//ideal displacement matrix de
	VectorXd idealdisp(24); 
		for (unsigned int i=0;i<24;i++ )
		{ if (LocationMatrix[i])
			idealdisp(i) = Displacement[LocationMatrix[i]-1];
		  else
		    idealdisp(i) = 0;
		}
	
	CHexMaterial* material = dynamic_cast<CHexMaterial*>(ElementMaterial);	// Pointer to material of the element
	
	//construct De
	double v=material->nu;
	double k = material->E * (1-v)/(1+v)/(1-2*v);

	MatrixXd D(6,6); // constitutive matrix
	D << 1,			v/(1-v),	v/(1-v),		0,					0,					0,
		v/(1-v),	1,			v/(1-v),		0,					0,					0,
		v/(1-v),	v/(1-v),    1,				0,					0,					0,
		0,			0,		    0,				(1-2*v)/2/(1-v),    0,					0,
		0,			0,			0,				0,					(1-2*v)/2/(1-v),	0,
		0,			0,			0,				0,					0,					(1-2*v)/2/(1-v);
	D=D*k;

	//construct coordinate matrix
	MatrixXd coorxyz(8,3);
	for (unsigned int i=0;i<8;i++)
	{ 
		for (unsigned int j=0;j<3;j++)
	   {
			coorxyz(i,j)=nodes[i]->XYZ[j];
	    }
	}


// construct jacobi matrix
	double xi;
	double eta;
	double zeta;
	double detj;

	double xi8[8];
	double eta8[8];
	double zeta8[8];

	xi8[0]=0.577350269189626;
	xi8[1]=0.577350269189626;
	xi8[2]=-0.577350269189626;
	xi8[3]=-0.577350269189626;
	xi8[4]=0.577350269189626;
	xi8[5]=0.577350269189626;
	xi8[6]=-0.577350269189626;
	xi8[7]=-0.577350269189626;
	
	eta8[0]=-0.577350269189626;
	eta8[1]=0.577350269189626;
	eta8[2]=0.577350269189626;
	eta8[3]=-0.577350269189626;
	eta8[4]=-0.577350269189626;
	eta8[5]=0.577350269189626;
	eta8[6]=0.577350269189626;
	eta8[7]=-0.577350269189626;

	zeta8[0]=-0.577350269189626;
	zeta8[1]=-0.577350269189626;
	zeta8[2]=-0.577350269189626;
	zeta8[3]=-0.577350269189626;
	zeta8[4]=0.577350269189626;
	zeta8[5]=0.577350269189626;
	zeta8[6]=0.577350269189626;
	zeta8[7]=0.577350269189626;

	
	MatrixXd GN(3,8);
	MatrixXd J;
	MatrixXd Jni;
	MatrixXd Bele; // elements in Be
	MatrixXd Be(6,24);
	MatrixXd stressXYZ(6,8); //8 gauss point, 6 stress 
	MatrixXd interpo(8,8);  //interpolation for stress recovery 
	MatrixXd recovery; //

	for (unsigned p=0;p<8;p++)
	{
		
		xi   = xi8[p];
		eta  = eta8[p];
		zeta = zeta8[p];
		

	GN << (1-eta)*(1-zeta), (1+eta)*(1-zeta),  -(1+eta)*(1-zeta),  -(1-eta)*(1-zeta),  (1-eta)*(1+zeta),  (1+eta)*(1+zeta),  -(1+eta)*(1+zeta),  -(1-eta)*(1+zeta),  
		  -(1+xi)*(1-zeta),  (1+xi)*(1-zeta),    (1-xi)*(1-zeta),  -(1-xi)*(1-zeta),   -(1+xi)*(1+zeta),   (1+xi)*(1+zeta),   (1-xi)*(1+zeta),   -(1-xi)*(1+zeta),  
		  -(1+xi)*(1-eta),   -(1+xi)*(1+eta),    -(1-xi)*(1+eta),   -(1-xi)*(1-eta),    (1+xi)*(1-eta),    (1+xi)*(1+eta),    (1-xi)*(1+eta),    (1-xi)*(1-eta);    
	GN=GN/8; // coefficient

	J=GN*coorxyz;
	Jni=J.inverse();
	detj=J.determinant();
	Bele=Jni*GN;
	

	// assign the value of Be	
	Be << Bele(0,0),0,0,				Bele(0,1), 0,0,				Bele(0,2), 0,0,				Bele(0,3), 0,0,				Bele(0,4), 0,0,				Bele(0,5), 0,0,				Bele(0,6), 0,0,				Bele(0,7), 0,0,
		  0,Bele(1,0),0,				0,Bele(1,1),0,				0,Bele(1,2),0,				0,Bele(1,3),0,				0,Bele(1,4),0,				0,Bele(1,5),0,				0,Bele(1,6),0,				0,Bele(1,7),0,
		  0,0,Bele(2,0),				0,0,Bele(2,1),				0,0,Bele(2,2),				0,0,Bele(2,3),				0,0,Bele(2,4),				0,0,Bele(2,5),				0,0,Bele(2,6),				0,0,Bele(2,7),
		  Bele(1,0),Bele(0,0),0,		Bele(1,1),Bele(0,1),0,		Bele(1,2),Bele(0,2),0,		Bele(1,3),Bele(0,3),0,		Bele(1,4),Bele(0,4),0,		Bele(1,5),Bele(0,5),0,		Bele(1,6),Bele(0,6),0,		Bele(1,7),Bele(0,7),0,
		  0,Bele(2,0),Bele(1,0),		0,Bele(2,1),Bele(1,1),		0,Bele(2,2),Bele(1,2),		0,Bele(2,3),Bele(1,3),		0,Bele(2,4),Bele(1,4),		0,Bele(2,5),Bele(1,5),		0,Bele(2,6),Bele(1,6),		0,Bele(2,7),Bele(1,7),
		  Bele(2,0),0,Bele(0,0),		Bele(2,1),0,Bele(0,1),		Bele(2,2),0,Bele(0,2),		Bele(2,3),0,Bele(0,3),		Bele(2,4),0,Bele(0,4),		Bele(2,5),0,Bele(0,5),		Bele(2,6),0,Bele(0,6),		Bele(2,7),0,Bele(0,7);
	

	// 6,1 for each gauss point stress 
	//a column of stressXYZ
	stressXYZ.col(p)=D*Be*idealdisp;  
	//loop for every gauss point 

	}
	
	// stress recovery for tress on nodes
	double a=2.549038105676658;
	double b=-0.683012701892219;
	double c=0.183012701892219;
	double d=-0.049038105676658;

	interpo << a, b, c, b, b, c, d, c,
			   b, a, b, c, c, b, c, d,
			   c, b, a, b, d, c, b, c,
			   b, c, b, a, c, d, c, b,
			   b, c, d, c, a, b, c, b,
			   c, b, c, d, b, a, b, c,
			   d, c, b, c, c, b, a, b,
			   c, d, c, b, b, c, b ,a;


// interpolate stress on nodes in order of xx,yy,zz,xy,yz,xz
	for (unsigned i=0;i<6;i++)
	{
		recovery = interpo*(stressXYZ.row(i).transpose());
		for (unsigned j=0;j<8;j++)
		{
			stressHex[6*j+i]=recovery(i);
		}
	}
	
}


