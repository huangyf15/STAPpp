
#include "PostOutputter.h"
#include "Domain.h"

PostOutputter* PostOutputter::_instance = nullptr;

PostOutputter::PostOutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}


//	Return the single instance of the class
PostOutputter* PostOutputter::Instance(string FileName)
{
	if (!_instance)
		_instance = new PostOutputter(FileName);
	return _instance;
}


// postprocess 
void PostOutputter::OutputElementStress()
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
						*this << ELENUM[Ele*2] << setw(4) << ELENUM[Ele*2+1] << endl;
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
