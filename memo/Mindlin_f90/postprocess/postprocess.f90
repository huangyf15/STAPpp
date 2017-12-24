
    SUBROUTINE POSTPROCESS(DIM,NNUM,FNUM,ID,X,Y,Z,U,ELID,FORCE)
    
    USE GLOBALS
    
    IMPLICIT NONE
    ! NNUM --- Number of nodes in one element
    ! ELID --- ID number for nodes in elements
    INTEGER, INTENT(IN)         ::  DIM,NNUM,FNUM,ID(6,NUMNP),ELID(NNUM,NPAR(2))
    INTEGER                     ::  I, J
    REAL(8), INTENT(IN)         ::  X(NUMNP),Y(NUMNP),Z(NUMNP),U(NEQ)
    REAL(8)                     ::  DISP(NUMNP*DIM)
    REAL(8)                     ::  FORCE(FNUM,NUMNP)
    
    CHARACTER*20                :: NAME
    
    !IF(STATUS .EQ. 0) THEN
        !OPEN(IPOST,FILE = "STAP90.DAT", STATUS = "REPLACE")
        !WRITE(IPOST, '("TITLE = """,A,"""")') HED
    !    STATUS = STATUS + 1
    !END IF
    
    IF(IND == 1) THEN
        IF(DIM == 2) THEN
            ! 4Q
            IF(NPAR(1) == 2) THEN
                NAME = 'FEQUADRILATERAL'
                CALL OUTPUT_2
            ! TRUSS
            ELSE IF(NPAR(1)==1) THEN
                NAME = 'FELINESEG'
                CALL OUTPUT_2
            ! TRIANGLE
            ELSE IF(NPAR(1)==3) THEN
                NAME = 'FETRIANGLE'
                CALL OUTPUT_2
            ! BEAM
            ELSE IF(NPAR(1)==5) THEN
                NAME = 'FELINESEG'
                CALL OUTPUT_2
            ! 8Q
            ELSE IF(NPAR(1)==8) THEN
                NAME = 'FEQUADRILATERAL'
                CALL OUTPUT_2
            ! 9Q
            ELSE IF(NPAR(1)==9) THEN
                NAME = 'FEQUADRILATERAL'
                CALL OUTPUT_2
            ! 6T
            ELSE IF(NPAR(1)==10) THEN
                NAME = 'FETRIANGLE'
                CALL OUTPUT_2

            END IF
        ELSE IF(DIM ==3) THEN
            ! 8H
            IF(NPAR(1)==4) THEN
                NAME = 'FEBRICK'
                CALL OUTPUT_3
            ! TRUSS
            ELSE IF(NPAR(1)==1) THEN
                NAME = 'FELINESEG'
                CALL OUTPUT_3
            ! BEAM
            ELSE IF(NPAR(1)==5) THEN
                NAME = 'FELINESEG'
                CALL OUTPUT_3
            ! PLATE
            ELSE IF(NPAR(1)==6) THEN
                NAME = 'FEQUADRILATERAL'
                CALL OUTPUT_3
            ! SHELL
            ELSE IF(NPAR(1)==7) THEN
                NAME = 'FEQUADRILATERAL'
                CALL OUTPUT_3
            ! 4T
            ELSE IF(NPAR(1)==11) THEN
                NAME = 'FETETRAHEDRON'
                CALL OUTPUT_3

                
            END IF
        END IF
    ELSE IF(IND == 3) THEN
        STATUS(2) = STATUS(2) + 1
        
        ! Generate the displacement
        DISP = 0.
        DO I = 1, NUMNP
            DO J = 1,DIM
                IF(ID(J,NUMNP) .NE. 0) THEN
                    DISP((I-1)*DIM + J) = U(ID(J,NUMNP))
                END IF
            END DO
        END DO
        
        IF(DIM == 2) THEN
            WRITE(IPOST, '("ZONE T=""STRAINED"",DATAPACKING=POINT,",\)')
            ! 4Q
            IF(NPAR(1) == 2) THEN
                WRITE(IPOST,'("ZONETYPE=FEQUADRILATERAL,NODES=" , I5 , ", ELEMENTS=", I5 ,", C=GREEN, VARSHARELIST=([1,2] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(5(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(:,I)
                END DO
            ! TRUSS
            ELSE IF(NPAR(1) == 1) THEN
                WRITE(IPOST,'("ZONETYPE=FELINESEG,NODES=" , I5 , ", ELEMENTS=", I5 ,", C=GREEN, VARSHARELIST=([1,2] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(5(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(1,I),0.D0,0.D0
                END DO
            ! 3T
            ELSE IF(NPAR(1) == 3) THEN
                WRITE(IPOST,'("ZONETYPE=FETRIANGLE,NODES=" , I5 , ", ELEMENTS=", I5 ,", C=GREEN, VARSHARELIST=([1,2] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(5(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(:,I)
                END DO
            ! BEAM
            ELSE IF(NPAR(1) == 5) THEN
                WRITE(IPOST,'("ZONETYPE=FELINESEG,NODES=" , I5 , ", ELEMENTS=", I5 ,", C=GREEN, VARSHARELIST=([1,2] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(5(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(:,I)
                END DO
            ! 8Q
            ELSE IF(NPAR(1) == 8) THEN
                WRITE(IPOST,'("ZONETYPE=FEQUADRILATERAL,NODES=" , I5 , ", ELEMENTS=", I5 ,", C=GREEN, VARSHARELIST=([1,2] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(5(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(:,I)
                END DO
            ! 9Q
            ELSE IF(NPAR(1) == 9) THEN
                WRITE(IPOST,'("ZONETYPE=FEQUADRILATERAL,NODES=" , I5 , ", ELEMENTS=", I5 ,", C=GREEN, VARSHARELIST=([1,2] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(5(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(:,I)
                END DO
            ! 6T
            ELSE IF(NPAR(1) == 10) THEN
                WRITE(IPOST,'("ZONETYPE=FETRIANGLE,NODES=" , I5 , ", ELEMENTS=", I5 ,", C=GREEN, VARSHARELIST=([1,2] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(5(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(:,I)
                END DO

                
            END IF
        ELSE IF(DIM ==3) THEN
            WRITE(IPOST, '("ZONE T=""STRAINED"",DATAPACKING=POINT,",\)')
            ! 8H
            IF(NPAR(1) == 4) THEN
                WRITE(IPOST,'("ZONETYPE=FEBRICK,NODES=" , I5 , ", ELEMENTS=", I5, ", C=GREEN, VARSHARELIST=([1,2,3] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(9(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(:,I)
                END DO
            ! TRUSS
            ELSE IF(NPAR(1) == 1) THEN
                WRITE(IPOST,'("ZONETYPE=FELINESEG,NODES=" , I5 , ", ELEMENTS=", I5, ", C=GREEN, VARSHARELIST=([1,2,3] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(9(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(1,I), (0.D0, J = 1,5)
                END DO
            ! BEAM
            ELSE IF(NPAR(1) == 5) THEN
                WRITE(IPOST,'("ZONETYPE=FELINESEG,NODES=" , I5 , ", ELEMENTS=", I5, ", C=GREEN, VARSHARELIST=([1,2,3] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(9(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(1,I), 0.D0, 0.D0, FORCE(2:3,I), 0.D0
                END DO
            ! PLATE
            ELSE IF(NPAR(1) == 6) THEN
                WRITE(IPOST,'("ZONETYPE=FEQUADRILATERAL,NODES=" , I5 , ", ELEMENTS=", I5, ", C=GREEN, VARSHARELIST=([1,2,3] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(9(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(1,I), 0.D0, 0.D0, FORCE(2:3,I), 0.D0
                END DO
            ! SHELL
            ELSE IF(NPAR(1) == 7) THEN
                WRITE(IPOST,'("ZONETYPE=FEQUADRILATERAL,NODES=" , I5 , ", ELEMENTS=", I5, ", C=GREEN, VARSHARELIST=([1,2,3] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(9(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(1,I), 0.D0, 0.D0, FORCE(2:3,I), 0.D0
                END DO
            ! 4T
            ELSE IF(NPAR(1) == 11) THEN
                WRITE(IPOST,'("ZONETYPE=FETETRAHEDRON,NODES=" , I5 , ", ELEMENTS=", I5, ", C=GREEN, VARSHARELIST=([1,2,3] = 1), CONNECTIVITYSHAREZONE = ",I5)') NUMNP, NPAR(2), STATUS(2)
                DO I = 1, NUMNP
                    WRITE(IPOST, '(9(E14.6,2X))') (DISP((I-1)*DIM + J), J = 1,DIM), FORCE(1,I), 0.D0, 0.D0, FORCE(2:3,I), 0.D0
                END DO

                
            END IF
        END IF
    END IF
    

    CONTAINS
    
    SUBROUTINE OUTPUT_2
    
    USE GLOBALS
    IMPLICIT NONE
    
    
                IF(STATUS(1) .EQ. 0) THEN
  WRITE(IPOST, '("TITLE = """,A80,"""")') HED
                    WRITE(IPOST, '("VARIABLES = ""X"",""Y"",""U"",""V"",&
                        ""STRESS_X"",""STRESS_Y"",""SHEAR""",/,&
                        "ZONE T=""INITIAL"", N=",I5,",DATAPACKING=POINT, ",\)') NUMNP
                    STATUS(1) = STATUS(1) + 1
                    WRITE(IPOST,'("ZONETYPE=",A20,",NODES=",I5,", ELEMENTS=",I5,", C=BLUE")') NAME,NUMNP,NPAR(2)
                    DO I = 1, NUMNP
                        WRITE(IPOST, '(7(E14.6,2X))')   X(I),Y(I), (0.D0, J=1,5)
                    END DO
                ELSE
                    WRITE(IPOST,'("ZONE T=""INITIAL"", N=",I5,", DATAPACKING=POINT, ZONETYPE=",A20,",NODES=",I5,", ELEMENTS=",I5,", C=BLUE, VARSHARELIST = ([1,2,3,4,5,6,7]=1)")') NAME,NUMNP,NPAR(2)
                END IF
                DO I = 1, NPAR(2)
                    WRITE(IPOST, '(4I5,1X)') ELID(:,I)
                END DO

    
    END SUBROUTINE
    
    SUBROUTINE OUTPUT_3
    
    USE GLOBALS
    IMPLICIT NONE
    
    
                IF(STATUS(1) .EQ. 0) THEN
  WRITE(IPOST, '("TITLE = """,A80,"""")') HED
                    WRITE(IPOST, '("VARIABLES = ""X"",""Y"",""Z"",""U"",""V"",""W"",&
                        ""STRESS_X"",""STRESS_Y"",""STRESS_Z"",""SHEAR_YZ"",""SHEAR_ZX"",""SHEAR_XY""",/,&
                        "ZONE T=""INITIAL"", N=",I5,", DATAPACKING=POINT, ",\)') NUMNP
                    STATUS(1) = STATUS(1) + 1
                    WRITE(IPOST,'("ZONETYPE=",A20,",NODES=",I5,", ELEMENTS=",I5,", C=BLUE")') NAME,NUMNP,NPAR(2)
                    DO I = 1, NUMNP
                        WRITE(IPOST, '(12(E14.6,2X))')   X(I),Y(I),Z(I), (0.D0, J=1,9)
                    END DO
                ELSE
                    WRITE(IPOST,'("ZONE T=""INITIAL"", DATAPACKING=POINT, ZONETYPE=",A20,",NODES=",I5,", ELEMENTS=",I5,", C=BLUE, VARSHARELIST = ([1,2,3,4,5,6,7,8,9,10,11,12]=1)")') NAME,NUMNP,NPAR(2)
                END IF
                DO I = 1, NPAR(2)
                    WRITE(IPOST, '(8I5,1X)') ELID(:,I)
                END DO

    
    END SUBROUTINE

    END SUBROUTINE POSTPROCESS

    