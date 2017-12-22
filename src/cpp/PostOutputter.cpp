
#include "PostOutputter.h"
#include "Domain.h"

#define POINTS_DATA_WIDTH 14

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

// Postprocess
void PostOutputter::OutputElementStress()
{
	// the amplification factor 
	double coeff = 1E4;

	CDomain* FEMData = CDomain::Instance();

	double* Displacement = FEMData->GetDisplacement();

	const unsigned int NUMEG = FEMData->GetNUMEG(); // number of element group

	// loop for each element group
	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{

		*this << "TITLE = \" STAPpp FEM \" " << endl
			<< "VARIABLES = \"X_POST\",\"Y_POST\",\"Z_POST\", \"STRESS_XX\",\"STRESS_YY\",\"STRESS_ZZ\",\"STRESS_XY\",\"STRESS_YZ\",\"STRESS_ZX\",  " << endl;

		// Get the ElementGroup and related infos
		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		ElementTypes ElementType = EleGrp.GetElementType();	// ElementType
		unsigned int NUME = EleGrp.GetNUME();				// Number of elements
		unsigned int NUMNP = FEMData->GetNUMNP();			// Number of node points

		switch (ElementType)
		{
		case ElementTypes::Bar: // Bar element

			*this << "ZONE T = \"SCENE1\", N = " << NUME * 8 << ", E = " << NUME
				<< ", F = FEPOINT , ET = BRICK, C = RED" << endl;

			double PrePositionBar[24];
			double PostPositionBar[24];
			double stressBar[48];

			// Loop for each element
			// Node infos in the present ZONE
			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				CElement& Element = EleGrp.GetElement(Ele);
				Element.ElementPostInfo(stressBar, Displacement, PrePositionBar, PostPositionBar);

				for (unsigned nodeIndex = 0; nodeIndex < 8; nodeIndex++)
				{
					for (unsigned DegOF = 0; DegOF < 3; DegOF++)
					{
						*this << setw(POINTS_DATA_WIDTH) << (1 - coeff) * PrePositionBar[3 * nodeIndex + DegOF] + coeff * PostPositionBar[3 * nodeIndex + DegOF];
					}
					for (unsigned DegOF = 0; DegOF < 6; DegOF++)
					{
						*this << setw(POINTS_DATA_WIDTH) << stressBar[6 * nodeIndex + DegOF];
					}
					*this << endl;
				}
			}

			// Node numbers corresponding to each element
			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				for (unsigned NumEleNode = 0; NumEleNode < 8; NumEleNode++)
				{
					*this << setw(POINTS_DATA_WIDTH) << Ele * 8 + NumEleNode + 1;
				}
				*this << endl;
			}
			break;

        case ElementTypes::Quadrilateral: // Quadrilateral element

            *this << "ZONE T = \"SCENE1\", N=" << NUME*4 << ",E =" << NUME
                  << " ,F = FEPOINT , ET = QUADRILATERAL, C = RED" << endl;

            double stress4Q[24];
            double PrePosition4Q[12];
            double PostPosition4Q[12];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {

                dynamic_cast<CQuadrilateral&>(
                    EleGrp.GetElement(Ele)).ElementPostInfo(stress4Q, Displacement, PrePosition4Q, PostPosition4Q);

                for (unsigned ni = 0; ni < 4; ++ni)
                {
                    for (unsigned dof = 0; dof < 3; ++dof)
                        *this << setw(POINTS_DATA_WIDTH) << PostPosition4Q[ni * 3 + dof];
                    for (unsigned dof = 0; dof < 6; ++dof)
                        *this << setw(POINTS_DATA_WIDTH) << stress4Q[ni * 6 + dof];
                    *this << std::endl;
                }
            }
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned NumEleNode = 0; NumEleNode < 4; NumEleNode++)
                {
                    *this << setw(POINTS_DATA_WIDTH) << Ele * 4 + NumEleNode + 1;
                }
                *this << endl;
            }
            *this << endl;

            break;

        case ElementTypes::Beam: // Beam element
            *this << "ZONE T = \"SCENE1\", N=" <<NUME*8 << ",E =" << NUME
                  << " ,F = FEPOINT , ET = BRICK, C = RED" << endl;

            double beamstress[48];
			double prePositionBeam[24];
            double postPositionBeam[24];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementPostInfo(beamstress, Displacement, prePositionBeam, postPositionBeam);

                CBeamMaterial& material =
                    *dynamic_cast<CBeamMaterial*>(Element.GetElementMaterial());

                for (unsigned i = 0; i < 8; i++){

					for (unsigned DegOF = 0; DegOF < 3; DegOF++){
						*this << setw(POINTS_DATA_WIDTH) << postPositionBeam[3 * i + DegOF];
					}

					for (unsigned DegOF = 0; DegOF < 6; DegOF++){
						*this << setw(POINTS_DATA_WIDTH) << beamstress[6 * i + DegOF];
					}
                    
					*this << endl;
				}
            }

            // Node numbers corresponding to each element
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned NumEleNode = 0; NumEleNode < 8; NumEleNode++)
                {
                    *this << setw(POINTS_DATA_WIDTH) << Ele * 8 + NumEleNode + 1;
                }
                *this << endl;
            }
            

            *this << endl;

            break;

        case ElementTypes::Triangle: // 3T element

            *this << "ZONE T = \"SCENE1\", N =" <<NUME * 3 << ",E =" << NUME
                  << " ,F = FEPOINT , ET = TRIANGLE, C = RED" << endl;

            double stress3T[3];
            double PrePosition3T[9];
            double PostPosition3T[9];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementPostInfo(stress3T, Displacement, PrePosition3T, PostPosition3T);
                CTriangleMaterial material =
                    *dynamic_cast<CTriangleMaterial*>(Element.GetElementMaterial());

                for (unsigned nodeIndex = 0; nodeIndex < 3; nodeIndex++)
                {
					for (unsigned dof = 0; dof < 3; dof++)
                        *this << setw(POINTS_DATA_WIDTH) << PostPosition3T[nodeIndex * 3 + dof];
                    *this << setw(POINTS_DATA_WIDTH) << stress3T[0]
                          << setw(POINTS_DATA_WIDTH) << stress3T[1]
                          << setw(POINTS_DATA_WIDTH) << 0.0
                          << setw(POINTS_DATA_WIDTH) << stress3T[2]
                          << setw(POINTS_DATA_WIDTH) << 0.0
                          << setw(POINTS_DATA_WIDTH) << 0.0 << std::endl;
                }
            }
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned _ = 0; _ < 3; _++)
                    *this << setw(POINTS_DATA_WIDTH) << Ele * 3 + _ + 1;
                *this << std::endl;
            }

            *this << endl;

            break;

        case ElementTypes::Hexahedron: // 8H element
        {
            *this << "ZONE T=\"SCENE1\", N=" << NUME*8 << ",E=" << NUME
                  << " ,F=FEPOINT , ET= BRICK, C= RED" << endl;

            double stressHex[48];
            double PrePosition8H[24];
            double Position8H[24];
            // for SPR
            unsigned int* glo2ET= new unsigned int[NUMNP]; 

            //call the SPR function 

            // for (unsigned int Ele = 0; Ele < NUME; Ele++)
            //{
            //  dynamic_cast<CHex&>(EleGrp.GetElement(Ele)).ElementPostInfoSPR(stressHex, Displacement, PrePosition8H, Position8H);
            //}

            // normal output way 

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementPostInfo(stressHex, Displacement, PrePosition8H, Position8H);

                CHexMaterial& material = *dynamic_cast<CHexMaterial*>(Element.GetElementMaterial());
                for (unsigned _ = 0; _ < 8; _++)
                {
                    *this << PrePosition8H[_*3 + 0]+coeff*(Position8H[_ * 3 + 0]-PrePosition8H[_*3 + 0]) << setw(POINTS_DATA_WIDTH)
                          << PrePosition8H[_*3 + 1]+coeff*(Position8H[_ * 3 + 1]-PrePosition8H[_*3 + 1]) << setw(POINTS_DATA_WIDTH)
                          << PrePosition8H[_*3 + 2]+coeff*(Position8H[_ * 3 + 2]-PrePosition8H[_*3 + 2])<< setw(POINTS_DATA_WIDTH) 
                          << stressHex[_*6 + 0] << setw(16)<< stressHex[_*6 + 1] << setw(16)<< stressHex[_*6 +2]<< setw(16)<< stressHex[_*6 +3]<< setw(16)<< stressHex[_*6 +4]<< setw(16)<< stressHex[_*6 +5]
                          << endl;
                }
            }
            for (unsigned int Ele= 0; Ele < NUME; Ele++) 
            {
                for (unsigned int i=1;i<9;i++)
                    *this << Ele*8+i<<setw(16);
                    *this<<endl;
            }
            *this << endl;

            break;
        }

        case ElementTypes::TimoshenkoSRINT: // TimoshenkoSRINT beam element
            *this << "ZONE T=\"SCENE1\", N=" <<NUME * 8 << ",E=" << NUME
                 << " ,F=FEPOINT , ET= BRICK, C= RED" << endl;

            double stressTimoSRINT[3];
            double PrePositionTimoSRINT[24];
            double PositionTimoSRINT[24];
            
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                dynamic_cast<CTimoshenkoSRINT&>(EleGrp.GetElement(Ele))
                    .ElementPostInfo(stressTimoSRINT, Displacement, PrePositionTimoSRINT, PositionTimoSRINT);

                for (unsigned _ = 0; _ < 8; _++)
                    *this << PositionTimoSRINT[_ * 3 + 0] << setw(POINTS_DATA_WIDTH)
                          << PositionTimoSRINT[_ * 3 + 1] << setw(POINTS_DATA_WIDTH)
                          << PositionTimoSRINT[_ * 3 + 2] << setw(POINTS_DATA_WIDTH) << endl;
                *this << std::endl;
            }

            break;

		case ElementTypes::TimoshenkoEBMOD: // TimoshenkoEBMOD beam element
			*this << "ZONE T=\"SCENE1\", N=" << NUME*8 << ",E=" << NUME
				<< " ,F=FEPOINT , ET= BRICK, C= RED" << endl;

			double stressTimoEBMOD[3];
			double PrePositionTimoEBMOD[24];
			double PositionTimoEBMOD[24];

			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				dynamic_cast<CTimoshenkoEBMOD&>(EleGrp.GetElement(Ele))
					.ElementPostInfo(stressTimoEBMOD, Displacement, PrePositionTimoEBMOD, PositionTimoEBMOD);

				for (unsigned _ = 0; _ < 8; _++)
					*this << PositionTimoEBMOD[_ * 3 + 0] << setw(POINTS_DATA_WIDTH)
					<< PositionTimoEBMOD[_ * 3 + 1] << setw(POINTS_DATA_WIDTH)
					<< PositionTimoEBMOD[_ * 3 + 2] << setw(POINTS_DATA_WIDTH) << endl;
				*this << std::endl;
			}

			break;

        case ElementTypes::Plate:
            *this << "ZONE T=\"SCENE1\", N=" << 8 * NUME << " E=" << NUME
                  << " F=FEPOINT , ET= BRICK, C= RED" << endl;

            double stresses4PE[48];
			double PrePositions4PE[24];
            double Positions4PE[24];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                dynamic_cast<CPlate&>(EleGrp.GetElement(Ele))
                    .ElementPostInfo(stresses4PE, Displacement, PrePositions4PE, Positions4PE);

                for (unsigned i = 0; i < 4; ++i)
                { // four gauss points
                    *this << setw(POINTS_DATA_WIDTH) << Positions4PE[3*i] << setw(POINTS_DATA_WIDTH) << Positions4PE[3 * i + 1]
                          << setw(POINTS_DATA_WIDTH) << Positions4PE[3 * i + 2] << setw(POINTS_DATA_WIDTH)
                          << stresses4PE[6 * i] << setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 1]
                          << setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 2] << setw(POINTS_DATA_WIDTH)
                          << stresses4PE[6 * i + 3] << setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 4]
                          << setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 5]
                          << endl;

                }
            }
            *this << endl;
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned int i = 0; i < 8; ++i){
                    *this << setw(10) << 8 * Ele + i + 1;
                }
                *this << std::endl;
            }
            break;

        case ElementTypes::Shell:
            *this << "ZONE T=\"SCENE1\", N=" << 8 * NUME << " E=" << NUME
                  << " F=FEPOINT , ET= BRICK, C= RED" << endl;

            double stresses4SE[48];
			double PrePostions4SE[24];
            double Positions4SE[24];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                dynamic_cast<CShell&>(EleGrp.GetElement(Ele))
                    .ElementPostInfo(stresses4SE, Displacement, PrePostions4SE, Positions4SE);

                for (unsigned i = 0; i < 8; ++i)
                {
                    // four gauss points;
                    *this << setw(POINTS_DATA_WIDTH) << Positions4SE[3*i] << setw(POINTS_DATA_WIDTH) << Positions4SE[3 * i + 1]
                          << setw(POINTS_DATA_WIDTH) << Positions4SE[3 * i + 2] << setw(POINTS_DATA_WIDTH)
                          << stresses4SE[6 * i] << setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 1]
                          << setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 2] << setw(POINTS_DATA_WIDTH)
                          << stresses4SE[6 * i + 3] << setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 4]
                          << setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 5]
                          << endl;

                }
            }
                        for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
            for (unsigned int i = 0; i < 8; ++i){
                    *this << setw(10) << 8 * Ele + i + 1;
                }
                *this << std::endl;
            }


            break;

        default: // Invalid element type
            cerr << "*** Error *** Element type " << ElementType << " has not been implemented.\n\n";
        }
    }
}
