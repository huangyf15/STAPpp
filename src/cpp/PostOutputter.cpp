
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

// Return the single instance of the class
PostOutputter* PostOutputter::Instance(string FileName)
{
    if (!_instance)
        _instance = new PostOutputter(FileName);
    return _instance;
}

// Postprocess
void PostOutputter::OutputElementStress()
{
    // The amplification factor coeff
    // Output = InitPosition + coeff * Displacement
    //        = InitPosition + coeff * (FinalPosition - InitPosition);
    double coeff = 1.0;

    CDomain* FEMData = CDomain::Instance();

    double* Displacement = FEMData->GetDisplacement();

    const unsigned int NUMEG = FEMData->GetNUMEG(); // Number of element groups

    *this << "TITLE = \" STAPpp FEM \" " << endl
          << "VARIABLES = \"X_POST\", \"Y_POST\", \"Z_POST\", "
             "\"STRESS_I\", \"STRESS_II\", \"STRESS_III\", \"STRESS_VONMISES\", "
             "\"STRESS_XX\", \"STRESS_YY\", \"STRESS_ZZ\", \"STRESS_XY\", \"STRESS_YZ\", \"STRESS_ZX\""
          << endl;

    // loop for each element group
    for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
    {

        // Get the ElementGroup and related infos
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
        ElementTypes ElementType = EleGrp.GetElementType(); // ElementType
        unsigned int NUME = EleGrp.GetNUME();               // Number of elements

        switch (ElementType)
        {
        case ElementTypes::Bar: // Bar element

            *this << "ZONE T = \"Bridge\", N = " << NUME * 8 << ", E = " << NUME
                  << ", F = FEPOINT , ET = BRICK, C = RED" << endl;

            double PrePositionBar[24];
            double PostPositionBar[24];
            double stressBar[48];
            double cmptStressBar[4];  // cmptStressBar = {stressI, stressII, stressIII, stress_vonMises};

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
                        *this << setw(POINTS_DATA_WIDTH)
                              << (1 - coeff) * PrePositionBar[3 * nodeIndex + DegOF] +
                                     coeff * PostPositionBar[3 * nodeIndex + DegOF];
                    }

                    cmptStressBar[0] = stressBar[6 * nodeIndex] + stressBar[6 * nodeIndex + 1] + stressBar[6 * nodeIndex + 2];
                    cmptStressBar[1] = stressBar[6 * nodeIndex]*stressBar[6 * nodeIndex + 1] - stressBar[6 * nodeIndex + 3]*stressBar[6 * nodeIndex + 3]
                                     + stressBar[6 * nodeIndex]*stressBar[6 * nodeIndex + 2] - stressBar[6 * nodeIndex + 5]*stressBar[6 * nodeIndex + 5]
                                     + stressBar[6 * nodeIndex + 1]*stressBar[6 * nodeIndex + 2] - stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex + 4];
                    cmptStressBar[2] = stressBar[6 * nodeIndex]*stressBar[6 * nodeIndex + 1]*stressBar[6 * nodeIndex + 2]
                                     + stressBar[6 * nodeIndex + 3]*stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex + 5]*2
                                     - stressBar[6 * nodeIndex + 1]*stressBar[6 * nodeIndex + 5]*stressBar[6 * nodeIndex + 5]
                                     - stressBar[6 * nodeIndex + 2]*stressBar[6 * nodeIndex + 3]*stressBar[6 * nodeIndex + 3]
                                     - stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex];
                    cmptStressBar[3] = sqrt(cmptStressBar[0]*cmptStressBar[0] - cmptStressBar[1]);
                    *this << setw(POINTS_DATA_WIDTH) << cmptStressBar[0]
                          << setw(POINTS_DATA_WIDTH) << cmptStressBar[1]
                          << setw(POINTS_DATA_WIDTH) << cmptStressBar[2]
                          << setw(POINTS_DATA_WIDTH) << cmptStressBar[3];

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

            *this << "ZONE T = \"Bridge\", N = " << NUME * 4 << ",E = " << NUME
                  << " ,F = FEPOINT , ET = QUADRILATERAL, C = RED" << endl;

            double stress4Q[24];
            double PrePosition4Q[12];
            double PostPosition4Q[12];
            double cmptStress4Q[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {

                EleGrp.GetElement(Ele).ElementPostInfo(stress4Q, Displacement, PrePosition4Q,
                                                       PostPosition4Q);

                for (unsigned ni = 0; ni < 4; ++ni)
                {
                    for (unsigned dof = 0; dof < 3; ++dof)
                        *this << setw(POINTS_DATA_WIDTH) << PostPosition4Q[ni * 3 + dof];

                    cmptStress4Q[0] = stress4Q[6 * ni] + stress4Q[6 * ni + 1] + stress4Q[6 * ni + 2];
                    cmptStress4Q[1] = stress4Q[6 * ni]*stress4Q[6 * ni + 1] - stress4Q[6 * ni + 3]*stress4Q[6 * ni + 3]
                                     + stress4Q[6 * ni]*stress4Q[6 * ni + 2] - stress4Q[6 * ni + 5]*stress4Q[6 * ni + 5]
                                     + stress4Q[6 * ni + 1]*stress4Q[6 * ni + 2] - stress4Q[6 * ni + 4]*stress4Q[6 * ni + 4];
                    cmptStress4Q[2] = stress4Q[6 * ni]*stress4Q[6 * ni + 1]*stress4Q[6 * ni + 2]
                                     + stress4Q[6 * ni + 3]*stress4Q[6 * ni + 4]*stress4Q[6 * ni + 5]*2
                                     - stress4Q[6 * ni + 1]*stress4Q[6 * ni + 5]*stress4Q[6 * ni + 5]
                                     - stress4Q[6 * ni + 2]*stress4Q[6 * ni + 3]*stress4Q[6 * ni + 3]
                                     - stress4Q[6 * ni + 4]*stress4Q[6 * ni + 4]*stress4Q[6 * ni];
                    cmptStress4Q[3] = sqrt(cmptStress4Q[0]*cmptStress4Q[0] - cmptStress4Q[1]);
                    *this << setw(POINTS_DATA_WIDTH) << cmptStress4Q[0]
                          << setw(POINTS_DATA_WIDTH) << cmptStress4Q[1]
                          << setw(POINTS_DATA_WIDTH) << cmptStress4Q[2]
                          << setw(POINTS_DATA_WIDTH) << cmptStress4Q[3];

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
            *this << "ZONE T = \"Bridge\", N = " << NUME * 8 << ",E = " << NUME
                  << " ,F = FEPOINT , ET = BRICK, C = RED" << endl;

            double beamstress[48];
            double prePositionBeam[24];
            double postPositionBeam[24];
            double cmptStressBeam[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementPostInfo(beamstress, Displacement, prePositionBeam,
                                        postPositionBeam);

                CBeamMaterial& material =
                    *dynamic_cast<CBeamMaterial*>(Element.GetElementMaterial());

                for (unsigned i = 0; i < 8; i++)
                {

                    for (unsigned DegOF = 0; DegOF < 3; DegOF++)
                    {
                        *this << setw(POINTS_DATA_WIDTH)
                              << (1 - coeff) * prePositionBeam[3 * i + DegOF] +
                                     coeff * postPositionBeam[3 * i + DegOF];
                    }

                    cmptStressBeam[0] = beamstress[6 * i] + beamstress[6 * i + 1] + beamstress[6 * i + 2];
                    cmptStressBeam[1] = beamstress[6 * i]*beamstress[6 * i + 1] - beamstress[6 * i + 3]*beamstress[6 * i + 3]
                                     + beamstress[6 * i]*beamstress[6 * i + 2] - beamstress[6 * i + 5]*beamstress[6 * i + 5]
                                     + beamstress[6 * i + 1]*beamstress[6 * i + 2] - beamstress[6 * i + 4]*beamstress[6 * i + 4];
                    cmptStressBeam[2] = beamstress[6 * i]*beamstress[6 * i + 1]*beamstress[6 * i + 2]
                                     + beamstress[6 * i + 3]*beamstress[6 * i + 4]*beamstress[6 * i + 5]*2
                                     - beamstress[6 * i + 1]*beamstress[6 * i + 5]*beamstress[6 * i + 5]
                                     - beamstress[6 * i + 2]*beamstress[6 * i + 3]*beamstress[6 * i + 3]
                                     - beamstress[6 * i + 4]*beamstress[6 * i + 4]*beamstress[6 * i];
                    cmptStressBeam[3] = sqrt(cmptStressBeam[0]*cmptStressBeam[0] - cmptStressBeam[1]);
                    *this << setw(POINTS_DATA_WIDTH) << cmptStressBeam[0]
                          << setw(POINTS_DATA_WIDTH) << cmptStressBeam[1]
                          << setw(POINTS_DATA_WIDTH) << cmptStressBeam[2]
                          << setw(POINTS_DATA_WIDTH) << cmptStressBeam[3];

                    for (unsigned DegOF = 0; DegOF < 6; DegOF++)
                    {
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

            *this << "ZONE T = \"Bridge\", N = " << NUME * 3 << ",E = " << NUME
                  << " ,F = FEPOINT , ET = TRIANGLE, C = RED" << endl;

            double stress3T[3];
            double PrePosition3T[9];
            double PostPosition3T[9];
            double cmptStress3T[4];

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

                    cmptStress3T[0] = stress3T[6 * nodeIndex] + stress3T[6 * nodeIndex + 1];
                    cmptStress3T[1] = stress3T[6 * nodeIndex]*stress3T[6 * nodeIndex + 1] - stress3T[6 * nodeIndex + 2]*stress3T[6 * nodeIndex + 2];
					cmptStress3T[2] = sqrt(cmptStress3T[0] * cmptStress3T[0] - cmptStress3T[1]);
                    cmptStress3T[3] = 0.0;
                    *this << setw(POINTS_DATA_WIDTH) << cmptStress3T[0]
                          << setw(POINTS_DATA_WIDTH) << cmptStress3T[1]
                          << setw(POINTS_DATA_WIDTH) << cmptStress3T[2]
                          << setw(POINTS_DATA_WIDTH) << cmptStress3T[3];

                    *this << setw(POINTS_DATA_WIDTH) << stress3T[0] << setw(POINTS_DATA_WIDTH)
                          << stress3T[1] << setw(POINTS_DATA_WIDTH) << 0.0
                          << setw(POINTS_DATA_WIDTH) << stress3T[2] << setw(POINTS_DATA_WIDTH)
                          << 0.0 << setw(POINTS_DATA_WIDTH) << 0.0 << std::endl;
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
            *this << "ZONE T= \"Bridge\", N = " << NUME * 8 << " ,E = " << NUME
                  << " ,F = FEPOINT , ET = BRICK, C = RED" << endl;

            double stressHex[48];
            double PrePosition8H[24];
            double Position8H[24];
            double cmptStress8H[4];
            // for SPR
            unsigned int* glo2ET = new unsigned int[NUME * 8];

            // call the SPR function

            // for (unsigned int Ele = 0; Ele < NUME; Ele++)
            //{
            //  dynamic_cast<CHex&>(EleGrp.GetElement(Ele)).ElementPostInfoSPR(stressHex,
            //  Displacement, PrePosition8H, Position8H);
            //}

            // normal output way

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementPostInfo(stressHex, Displacement, PrePosition8H, Position8H);

                CHexMaterial& material = *dynamic_cast<CHexMaterial*>(Element.GetElementMaterial());
                for (unsigned _ = 0; _ < 8; _++)
                {
                    *this << setw(POINTS_DATA_WIDTH)
                          << PrePosition8H[_ * 3 + 0] +
                                 coeff * (Position8H[_ * 3 + 0] - PrePosition8H[_ * 3 + 0])
                          << setw(POINTS_DATA_WIDTH)
                          << PrePosition8H[_ * 3 + 1] +
                                 coeff * (Position8H[_ * 3 + 1] - PrePosition8H[_ * 3 + 1])
                          << setw(POINTS_DATA_WIDTH)
                          << PrePosition8H[_ * 3 + 2] +
                                 coeff * (Position8H[_ * 3 + 2] - PrePosition8H[_ * 3 + 2]);

                    cmptStress8H[0] = stressHex[6 * _] + stressHex[6 * _ + 1] + stressHex[6 * _ + 2];
                    cmptStress8H[1] = stressHex[6 * _]*stressHex[6 * _ + 1] - stressHex[6 * _ + 3]*stressHex[6 * _ + 3]
                                     + stressHex[6 * _]*stressHex[6 * _ + 2] - stressHex[6 * _ + 5]*stressHex[6 * _ + 5]
                                     + stressHex[6 * _ + 1]*stressHex[6 * _ + 2] - stressHex[6 * _ + 4]*stressHex[6 * _ + 4];
                    cmptStress8H[2] = stressHex[6 * _]*stressHex[6 * _ + 1]*stressHex[6 * _ + 2]
                                     + stressHex[6 * _ + 3]*stressHex[6 * _ + 4]*stressHex[6 * _ + 5]*2
                                     - stressHex[6 * _ + 1]*stressHex[6 * _ + 5]*stressHex[6 * _ + 5]
                                     - stressHex[6 * _ + 2]*stressHex[6 * _ + 3]*stressHex[6 * _ + 3]
                                     - stressHex[6 * _ + 4]*stressHex[6 * _ + 4]*stressHex[6 * _];
                    cmptStress8H[3] = sqrt(cmptStress8H[0]*cmptStress8H[0] - cmptStress8H[1]);
                    *this << setw(POINTS_DATA_WIDTH) << cmptStress8H[0]
                          << setw(POINTS_DATA_WIDTH) << cmptStress8H[1]
                          << setw(POINTS_DATA_WIDTH) << cmptStress8H[2]
                          << setw(POINTS_DATA_WIDTH) << cmptStress8H[3];

                    *this << setw(POINTS_DATA_WIDTH) << stressHex[_ * 6 + 0]
                          << setw(POINTS_DATA_WIDTH) << stressHex[_ * 6 + 1]
                          << setw(POINTS_DATA_WIDTH) << stressHex[_ * 6 + 2]
                          << setw(POINTS_DATA_WIDTH) << stressHex[_ * 6 + 3]
                          << setw(POINTS_DATA_WIDTH) << stressHex[_ * 6 + 4]
                          << setw(POINTS_DATA_WIDTH) << stressHex[_ * 6 + 5] << endl;
                }
            }
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned int i = 1; i < 9; i++)
                {
                    *this << setw(POINTS_DATA_WIDTH) << Ele * 8 + i;
                }
                *this << endl;
            }
            *this << endl;

            break;
        }

        case ElementTypes::TimoshenkoSRINT: // TimoshenkoSRINT beam element

            break;

        case ElementTypes::TimoshenkoEBMOD: // TimoshenkoEBMOD beam element

            break;

        case ElementTypes::Plate:
            *this << "ZONE T = \"Bridge\", N = " << 8 * NUME << " E = " << NUME
                  << " F = FEPOINT , ET = BRICK, C = RED" << endl;

            double stresses4PE[48];
            double PrePositions4PE[24];
            double Positions4PE[24];
            double cmptStressPlate[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                EleGrp.GetElement(Ele).ElementPostInfo(stresses4PE, Displacement, PrePositions4PE,
                                                       Positions4PE);

                for (unsigned i = 0; i < 4; ++i)
                { // four gauss points
					*this << setw(POINTS_DATA_WIDTH)
						<< (1 - coeff) * PrePositions4PE[3 * i] + coeff * Positions4PE[3 * i]
						<< setw(POINTS_DATA_WIDTH)
						<< (1 - coeff) * PrePositions4PE[3 * i + 1] +
						coeff * Positions4PE[3 * i + 1]
						<< setw(POINTS_DATA_WIDTH)
						<< (1 - coeff) * PrePositions4PE[3 * i + 2] +
						coeff * Positions4PE[3 * i + 2];

					cmptStressPlate[0] = stresses4PE[6 * i] + stresses4PE[6 * i + 1] + stresses4PE[6 * i + 2];
					cmptStressPlate[1] = stresses4PE[6 * i] * stresses4PE[6 * i + 1] - stresses4PE[6 * i + 3] * stresses4PE[6 * i + 3]
						+ stresses4PE[6 * i] * stresses4PE[6 * i + 2] - stresses4PE[6 * i + 5] * stresses4PE[6 * i + 5]
						+ stresses4PE[6 * i + 1] * stresses4PE[6 * i + 2] - stresses4PE[6 * i + 4] * stresses4PE[6 * i + 4];
                    cmptStressPlate[2] = stresses4PE[6 * i]*stresses4PE[6 * i + 1]*stresses4PE[6 * i + 2]
                                     + stresses4PE[6 * i + 3]*stresses4PE[6 * i + 4]*stresses4PE[6 * i + 5]*2
                                     - stresses4PE[6 * i + 1]*stresses4PE[6 * i + 5]*stresses4PE[6 * i + 5]
                                     - stresses4PE[6 * i + 2]*stresses4PE[6 * i + 3]*stresses4PE[6 * i + 3]
                                     - stresses4PE[6 * i + 4]*stresses4PE[6 * i + 4]*stresses4PE[6 * i];
					cmptStressPlate[2] = sqrt(cmptStressPlate[0] * cmptStressPlate[0] - cmptStressPlate[1]);
					*this << setw(POINTS_DATA_WIDTH) << cmptStressPlate[0]
						<< setw(POINTS_DATA_WIDTH) << cmptStressPlate[1]
						<< setw(POINTS_DATA_WIDTH) << cmptStressPlate[2]
						<< setw(POINTS_DATA_WIDTH) << cmptStressPlate[3];

					*this << setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i]
						<< setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 1]
						<< setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 2]
						<< setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 3]
						<< setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 4]
						<< setw(POINTS_DATA_WIDTH) << stresses4PE[6 * i + 5] << endl;
                }
            }
            *this << endl;
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned int i = 0; i < 8; ++i)
                {
                    *this << setw(POINTS_DATA_WIDTH) << 8 * Ele + i + 1;
                }
                *this << std::endl;
            }
            break;

        case ElementTypes::Shell:
            *this << "ZONE T = \"Bridge\", N = " << 8 * NUME << " E = " << NUME
                  << " F = FEPOINT , ET = BRICK, C = RED" << endl;

            double stresses4SE[48];
            double PrePositions4SE[24];
            double Positions4SE[24];
            double cmptStressShell[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                EleGrp.GetElement(Ele).ElementPostInfo(stresses4SE, Displacement, PrePositions4SE,
                                                       Positions4SE);

                for (unsigned i = 0; i < 8; ++i)
                {
                    // four gauss points;
					*this << setw(POINTS_DATA_WIDTH)
						<< (1 - coeff) * PrePositions4SE[3 * i] + coeff * Positions4SE[3 * i]
						<< setw(POINTS_DATA_WIDTH)
						<< (1 - coeff) * PrePositions4SE[3 * i + 1] +
						coeff * Positions4SE[3 * i + 1]
						<< setw(POINTS_DATA_WIDTH)
						<< (1 - coeff) * PrePositions4SE[3 * i + 2] +
						coeff * Positions4SE[3 * i + 2];

					cmptStressShell[0] = stresses4SE[6 * i] + stresses4SE[6 * i + 1] + stresses4SE[6 * i + 2];
					cmptStressShell[1] = stresses4SE[6 * i] * stresses4SE[6 * i + 1] - stresses4SE[6 * i + 3] * stresses4SE[6 * i + 3]
						+ stresses4SE[6 * i] * stresses4SE[6 * i + 2] - stresses4SE[6 * i + 5] * stresses4SE[6 * i + 5]
						+ stresses4SE[6 * i + 1] * stresses4SE[6 * i + 2] - stresses4SE[6 * i + 4] * stresses4SE[6 * i + 4];
                    cmptStressShell[2] = stresses4SE[6 * i]*stresses4SE[6 * i + 1]*stresses4SE[6 * i + 2]
                                     + stresses4SE[6 * i + 3]*stresses4SE[6 * i + 4]*stresses4SE[6 * i + 5]*2
                                     - stresses4SE[6 * i + 1]*stresses4SE[6 * i + 5]*stresses4SE[6 * i + 5]
                                     - stresses4SE[6 * i + 2]*stresses4SE[6 * i + 3]*stresses4SE[6 * i + 3]
                                     - stresses4SE[6 * i + 4]*stresses4SE[6 * i + 4]*stresses4SE[6 * i];
					cmptStressShell[3] = sqrt(cmptStressShell[0] * cmptStressShell[0] - cmptStressShell[1]);
					*this << setw(POINTS_DATA_WIDTH) << cmptStressShell[0]
						<< setw(POINTS_DATA_WIDTH) << cmptStressShell[1]
						<< setw(POINTS_DATA_WIDTH) << cmptStressShell[2]
						<< setw(POINTS_DATA_WIDTH) << cmptStressShell[3];

					*this << setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i]
						<< setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 1]
						<< setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 2]
						<< setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 3]
						<< setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 4]
						<< setw(POINTS_DATA_WIDTH) << stresses4SE[6 * i + 5] << endl;
                }
            }
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned int i = 0; i < 8; ++i)
                {
                    *this << setw(POINTS_DATA_WIDTH) << 8 * Ele + i + 1;
                }
                *this << std::endl;
            }

            break;

        case ElementTypes::T9Q:

            break;

        default: // Invalid element type
            cerr << "*** Error *** Element type " << ElementType
                 << " has not been implemented.\n\n";
        }
    }
}
