
#include "PostOutputter.h"
#include "Domain.h"

#define POINTS_DATA_WIDTH 15

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
        *this << "TITLE = \" 3DFEA \" " << endl
              << "VARIABLES = VARIABLES = \"X\",\"Y\",\"Z\", \"STRESS\"  " << endl;

        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
        unsigned int NUME = EleGrp.GetNUME();     // NUMBER OF ELEMENT
        unsigned int NUMNP = FEMData->GetNUMNP(); // NUMBER OF NODE POINTS

        ElementTypes ElementType = EleGrp.GetElementType();

        switch (ElementType)
        {
        case ElementTypes::Bar: // Bar element
        {
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << ", E=" << NUME
                  << " F=FEPOINT , ET= BRICK, C= RED" << endl;

            double stress;
            double Positions[8];
            unsigned int* ELENUM = new unsigned int[NUME * 2];
            //********************88double Positions[6];
            // loop for every element
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementStress2(&stress, Displacement, Positions);

                CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());

                for (unsigned _ = 0; _ < 2; _++)
                    *this << setw(POINTS_DATA_WIDTH) << Positions[3 * _ + 0]
                          << setw(POINTS_DATA_WIDTH) << Positions[3 * _ + 1]
                          << setw(POINTS_DATA_WIDTH) << Positions[3 * _ + 2] << endl;

                ELENUM[Ele * 2] = (unsigned int)Positions[6];
                ELENUM[Ele * 2 + 1] = (unsigned int)Positions[7];
            }
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                *this << setw(POINTS_DATA_WIDTH) << ELENUM[Ele * 2]
                      << setw(POINTS_DATA_WIDTH) << ELENUM[Ele * 2 + 1] << endl;
            }

            //*this << endl;
            break;
        }

        case ElementTypes::Quadrilateral: // Quadrilateral element
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= QUADRILATERAL, C= RED" << endl;

            double stresses[12];
            double Position4Q[12];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {

                dynamic_cast<CQuadrilateral&>(
                    EleGrp.GetElement(Ele)).ElementStress2(stresses, Displacement, Position4Q);

                for (unsigned i = 0; i < 4; ++i)
                { // four gauss points

                    *this << setw(POINTS_DATA_WIDTH) << Position4Q[i * 3]
                          << setw(POINTS_DATA_WIDTH) << Position4Q[i * 3 + 1]
                          << setw(POINTS_DATA_WIDTH) << Position4Q[i * 3 + 2] << std::endl;
                }
            }
            *this << endl;

            break;

        case ElementTypes::Beam: // Beam element
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= BRICK, C= RED" << endl;

            double beamstress[3];
            double PositionBeam[24];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementStress2(beamstress, Displacement, PositionBeam);

                CBeamMaterial& material =
                    *dynamic_cast<CBeamMaterial*>(Element.GetElementMaterial());

                for (unsigned _ = 0; _ < 8; _++)
                    *this << setw(POINTS_DATA_WIDTH) << PositionBeam[3 * _ + 0]
                          << setw(POINTS_DATA_WIDTH) << PositionBeam[3 * _ + 1]
                          << setw(POINTS_DATA_WIDTH) << PositionBeam[3 * _ + 2]  << endl;
            }

            *this << endl;
            break;

        case ElementTypes::Triangle: // 3T element
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= TRIANGLE, C= RED" << endl;
            double stress3T[3];
            double Position3T[9];
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementStress2(stress3T, Displacement, Position3T);
                CTriangleMaterial material =
                    *dynamic_cast<CTriangleMaterial*>(Element.GetElementMaterial());

                for (unsigned _ = 0; _ < 3; _++)
                    *this << Position3T[_ * 3 + 0] << setw(POINTS_DATA_WIDTH)
                          << Position3T[_ * 3 + 1] << setw(POINTS_DATA_WIDTH)
                          << Position3T[_ * 3 + 2] << setw(POINTS_DATA_WIDTH) << endl;
            }

            *this << endl;
            break;

        case ElementTypes::Hexahedron: // 8H element
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= BRICK, C= RED" << endl;

            double stressHex[48];
            double Position8H[24];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp.GetElement(Ele);
                Element.ElementStress2(stressHex, Displacement, Position8H);

                CHexMaterial& material = *dynamic_cast<CHexMaterial*>(Element.GetElementMaterial());
                for (unsigned _ = 0; _ < 8; _++)
                    *this << Position8H[_ * 3 + 0] << setw(POINTS_DATA_WIDTH)
                          << Position8H[_ * 3 + 1] << setw(POINTS_DATA_WIDTH)
                          << Position8H[_ * 3 + 2] << setw(POINTS_DATA_WIDTH) << endl;
            }

            *this << endl;
            break;

        case ElementTypes::TimoshenkoSRINT: // TimoshenkoSRINT beam element
            double TimoshenkoStresses[3];
            double TimoshenkoForces[12];
            double PositionTB[24];
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= BRICK, C= RED" << endl;
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                dynamic_cast<CTimoshenkoSRINT&>(EleGrp.GetElement(Ele))
                    .ElementStress2(TimoshenkoStresses, TimoshenkoForces, Displacement, PositionTB);

                for (unsigned _ = 0; _ < 8; _++)
                    *this << PositionTB[_ * 3 + 0] << setw(POINTS_DATA_WIDTH)
                          << PositionTB[_ * 3 + 1] << setw(POINTS_DATA_WIDTH)
                          << PositionTB[_ * 3 + 2] << setw(POINTS_DATA_WIDTH) << endl;
                *this << std::endl;
            }
            break;

        case ElementTypes::TimoshenkoEBMOD: // TimoshenkoEBMOD beam element
            double TimoshenkoEBStresses[3];
            double TimoshenkoEBForces[12];
            double PositionTEB[24];
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= BRICK, C= RED" << endl;

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                dynamic_cast<CTimoshenkoEBMOD&>(EleGrp.GetElement(Ele))
                    .ElementStress2(TimoshenkoEBStresses, TimoshenkoEBForces, Displacement,
                                    PositionTEB);

                for (unsigned _ = 0; _ < 8; _++)
                    *this << PositionTEB[_ * 3 + 0] << setw(POINTS_DATA_WIDTH)
                          << PositionTEB[_ * 3 + 1] << setw(POINTS_DATA_WIDTH)
                          << PositionTEB[_ * 3 + 2] << setw(POINTS_DATA_WIDTH) << endl;

                *this << std::endl;
            }
            break;

        case ElementTypes::Plate:
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= QUADRILATERAL, C= RED" << endl;

            double stresses4PE[12];
            double Positions4PE[12];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                dynamic_cast<CPlate&>(EleGrp.GetElement(Ele))
                    .ElementStress2(stresses4PE, Displacement, Positions4PE);

                for (unsigned i = 0; i < 4; ++i)
                { // four gauss points
                    *this << Positions4PE[0] << setw(POINTS_DATA_WIDTH) << Positions4PE[1]
                          << setw(POINTS_DATA_WIDTH) << Positions4PE[2] << setw(POINTS_DATA_WIDTH)
                          << endl
                          << Positions4PE[3] << setw(POINTS_DATA_WIDTH) << Positions4PE[4]
                          << setw(POINTS_DATA_WIDTH) << Positions4PE[5] << setw(POINTS_DATA_WIDTH)
                          << endl
                          << Positions4PE[6] << setw(POINTS_DATA_WIDTH) << Positions4PE[7]
                          << setw(POINTS_DATA_WIDTH) << Positions4PE[8] << setw(POINTS_DATA_WIDTH)
                          << endl
                          << Positions4PE[9] << setw(POINTS_DATA_WIDTH) << Positions4PE[10]
                          << setw(POINTS_DATA_WIDTH) << Positions4PE[11] << setw(POINTS_DATA_WIDTH)
                          << endl;

                    *this << std::endl;
                }
            }
            *this << endl;

            break;

        case ElementTypes::Shell:
            *this << "ZONE T=\"SCENE1\", N=" << NUMNP << "E=" << NUME
                  << " F=FEPOINT , ET= QUADRILATERAL, C= RED" << endl;
            double stresses4SE[15];
            double Positions4SE[15];
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {

                dynamic_cast<CShell&>(EleGrp.GetElement(Ele))
                    .ElementStress2(stresses4SE, Displacement, Positions4SE);

                for (unsigned i = 0; i < 5; ++i)
                {
                    // four gauss points;
                    // THE FIFTH POINT IS THE CENTRE POINT FOR IN-PLANE STRESSES
                    *this << Positions4SE[0] << setw(POINTS_DATA_WIDTH) << Positions4SE[1]
                          << setw(POINTS_DATA_WIDTH) << Positions4SE[2] << setw(POINTS_DATA_WIDTH)
                          << endl
                          << Positions4SE[3] << setw(POINTS_DATA_WIDTH) << Positions4SE[4]
                          << setw(POINTS_DATA_WIDTH) << Positions4SE[5] << setw(POINTS_DATA_WIDTH)
                          << endl
                          << Positions4SE[6] << setw(POINTS_DATA_WIDTH) << Positions4SE[7]
                          << setw(POINTS_DATA_WIDTH) << Positions4SE[8] << setw(POINTS_DATA_WIDTH)
                          << endl
                          << Positions4SE[9] << setw(POINTS_DATA_WIDTH) << Positions4SE[10]
                          << setw(POINTS_DATA_WIDTH) << Positions4SE[11] << setw(POINTS_DATA_WIDTH)
                          << endl;
                    // *this << setw(32) << stresses[i] << std::endl;

                    *this << std::endl;
                }
            }
            break;

        default: // Invalid element type
            cerr << "*** Error *** Elment type " << ElementType << " has not been implemented.\n\n";
        }
    }
}
