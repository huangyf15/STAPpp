SUBROUTINE BEAMS
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the BEAM element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106,N107,N108,N109,N110,N111,N112,N113,N114,N115,N116,N117,N118

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 8*NUMMAT*ITWO + 2*NUMMAT + 17*NUME + 6*NUME*ITWO + NUMNP + 6*NUMNP*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)
				  ! THE POINTER SHOULD CHANGE FOR DIFFERENT ELEMENT GROUP!
				  
! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: AREA(NUMMAT)
! N103: LM(12,NUME)
! N104: XYZ(6,NUME)
! N105: MTAP(NUME)
! N106: INDX(NUMMAT)
! N107: G(NUMMAT)
! N108: EIY(NUMMAT)
! N109: EIZ(NUMMAT)
! N110: HGI(NUME)
! N111: HGJ(NUME)
! N112: DENSITY(NUMMAT)
! N113: YMAX(NUMMAT)
! N114: ZMAX(NUMMAT)
! N115: ELID(2,NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+12*NUME
  N105=N104+6*NUME*ITWO
  N106=N105+NUME
  N107=N106+NUMMAT
  N108=N107+NUMMAT*ITWO
  N109=N108+NUMMAT*ITWO
  N110=N109+NUMMAT*ITWO
  N111=N110+NUME
  N112=N111+NUME
  N113=N112+NUMMAT*ITWO
  N114=N113+NUMMAT*ITWO
  N115=N114+NUMMAT*ITWO
  N116=N115+2*NUME
  N117=N116+NUMNP
  N118=N117+6*NUMNP*ITWO
  NLAST=N118

  MIDEST=NLAST - NFIRST

  CALL BEAM (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105),A(N106),A(N107),A(N108),A(N109),A(N110),A(N111),A(N112),A(N113),A(N114),A(N115),A(N116),A(N117))

  RETURN

END SUBROUTINE BEAMS


SUBROUTINE BEAM (ID,X,Y,Z,U,MHT,E,AREA,LM,XYZ,MATP,INDX,G,EIY,EIZ,HGI,HGJ,DENSITY,YMAX,ZMAX,ELID,NNUM,FORCEA)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(12,NPAR(2)),MATP(NPAR(2)),MHT(NEQ),INDX(NPAR(3)),HGI(NPAR(2)),HGJ(NPAR(2)),ELID(2,NPAR(2)),NNUM(NUMNP)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),AREA(NPAR(3)),  &
             XYZ(6,NPAR(2)),U(NEQ),G(NPAR(3)),EIY(NPAR(3)),EIZ(NPAR(3)),DENSITY(NPAR(3)),YMAX(NPAR(3)),ZMAX(NPAR(3))
  REAL(8) :: S(12,12),SS(12,12),K,GJ,GAKL,COR(6,NPAR(2)),DISP(12),D(12),FORCE(12),WEIGHT(3), R(NEQ), VONMISES, FORCEA(6,NUMNP)
  INTEGER :: NPAR1, NUME, NUMMAT, ND, I, J, L, N, II
  INTEGER :: MTYPE, IPRINT, MODE
  REAL(8) :: XL3, XL2, XL, XX, YY, STR, P,KS(12,12),TRANS(12,12)
  REAL(8) :: SIGMAX,SIGMAY,SIGMAZ,TAUXY,TAUXZ,TAUYZ
  
  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=12

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
                   '     EQ.2, QUADS ELEMENTS',/,  &
                   '     EQ.3, BEAM ELEMENTS',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S     CROSS-SECTIONAL',/,  &
                   ' NUMBER     MODXULUS',10X,'AREA',/,  &
                   15 X,'E',14X,'A')")

     DO I=1,NUMMAT
        READ (IIN,'(I10,2F20.0,I10,6F20.0)') N,E(N),AREA(N),INDX(N),G(N),EIY(N),EIZ(N),DENSITY(N),YMAX(N),ZMAX(N) ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6,2X,I5,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6)") N,E(N),AREA(N),INDX(N),G(N),EIY(N),EIZ(N),DENSITY(N),YMAX(N),ZMAX(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE       MATERIAL    HINGED?    HINGED?',/,   &
                      ' NUMBER-N      I        J       SET NUMBER      I          J')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(6I10,6F20.0)') N,I,J,MTYPE,HGI(N),HGJ(N),(COR(II,N),II=1,6)  ! Read in element information

        ELID(:,N) = (/I,J/)

!       Save element information
        XYZ(1,N)=X(I)  ! Coordinates of the element's left node
        XYZ(2,N)=Y(I)
        XYZ(3,N)=Z(I)

        XYZ(4,N)=X(J)  ! Coordinates of the element's right node
        XYZ(5,N)=Y(J)
        XYZ(6,N)=Z(J)

        MATP(N)=MTYPE  ! Material type

        DO L=1,12
           LM(L,N)=0
        END DO

        DO L=1,6
           LM(L,N)=ID(L,I)     ! Connectivity matrix
           LM(L+6,N)=ID(L,J)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,7X,I5,7X,I5,7X,I5)") N,I,J,MTYPE,HGI(N),HGJ(N)

     END DO

     CALL POSTPROCESS(3,2,6,ID,X,Y,Z,U,ELID,FORCEA)

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

    ! Gravity
        REWIND ILOAD
        READ(ILOAD)R

     DO N=1,NUME
        MTYPE=MATP(N)

		
		GJ=G(MTYPE)*(EIY(MTYPE)+EIZ(MTYPE))/E(MTYPE)
	    
		CALL  TRANSPORT(TRANS,XL,XYZ(1,N),COR(1,N))! TRANSPORTION MATRIX
		XL2=XL**2
		XL3=XL**3
		
		!add gravity force
		WEIGHT = DENSITY(MTYPE)*XL*AREA(MTYPE)*GRAVITY/2.
		R(LM(1,N))=R(LM(1,N))+WEIGHT(1)
		R(LM(2,N))=R(LM(2,N))+WEIGHT(2)
		R(LM(3,N))=R(LM(3,N))+WEIGHT(3)
		R(LM(7,N))=R(LM(7,N))+WEIGHT(1)
		R(LM(8,N))=R(LM(8,N))+WEIGHT(2)
		R(LM(9,N))=R(LM(9,N))+WEIGHT(3)
		
		SS=0

			IF(HGI(N).EQ.1)THEN
				IF(HGJ(N).EQ.1)THEN
					SS(1,1) =  E(MTYPE)*AREA(MTYPE)/XL
					SS(1,7) = -SS(1,1)
					SS(7,1) = -SS(1,1)
					SS(7,7) =  SS(1,1)
				ELSE
					SS(1,1) = E(MTYPE)*AREA(MTYPE)/XL
					SS(2,2) = 3*EIZ(MTYPE)/XL3
					SS(3,3) = 3*EIY(MTYPE)/XL3
					
					SS(1,7) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(2,8) = -3*EIZ(MTYPE)/XL3
					SS(3,9) = -3*EIY(MTYPE)/XL3
					
					SS(7,1) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(8,2) = -3*EIZ(MTYPE)/XL3
					SS(9,3) = -3*EIY(MTYPE)/XL3
					
					SS(8,8)=   3*EIZ(MTYPE)/XL3
					SS(9,9)=   3*EIY(MTYPE)/XL3
					
					SS(2,12) =  3*EIZ(MTYPE)/XL2
					SS(3,11) = -3*EIY(MTYPE)/XL2
					
					SS(12,2) =  3*EIZ(MTYPE)/XL2
					SS(11,3) = -3*EIY(MTYPE)/XL2
					
					SS(11,11) =  3*EIY(MTYPE)/XL
					SS(12,12) =  3*EIZ(MTYPE)/XL
					
					SS(11,9) =  3*EIY(MTYPE)/XL2
					SS(12,8) = -3*EIZ(MTYPE)/XL2
					
					SS(9,11) =  3*EIY(MTYPE)/XL2
					SS(8,12) = -3*EIZ(MTYPE)/XL2
				
				ENDIF
			ELSE
				IF(HGJ(N).EQ.1)THEN
					SS(1,1) = E(MTYPE)*AREA(MTYPE)/XL
					SS(2,2) = 3*EIZ(MTYPE)/XL3
					SS(3,3) = 3*EIY(MTYPE)/XL3
					
					SS(1,7) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(2,8) = -3*EIZ(MTYPE)/XL3
					SS(3,9) = -3*EIY(MTYPE)/XL3
					
					SS(7,1) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(8,2) = -3*EIZ(MTYPE)/XL3
					SS(9,3) = -3*EIY(MTYPE)/XL3
					
					SS(8,8)=   3*EIZ(MTYPE)/XL3
					SS(9,9)=   3*EIY(MTYPE)/XL3
					
					SS(2,6) =  3*EIZ(MTYPE)/XL2
					SS(3,5) = -3*EIY(MTYPE)/XL2
					
					SS(6,2) =  3*EIZ(MTYPE)/XL2
					SS(5,3) = -3*EIY(MTYPE)/XL2
					
					SS(5,5) =  3*EIY(MTYPE)/XL
					SS(6,6) =  3*EIZ(MTYPE)/XL
					
					SS(5,9) =  3*EIY(MTYPE)/XL2
					SS(6,8) = -3*EIZ(MTYPE)/XL2
					
					SS(9,5) =  3*EIY(MTYPE)/XL2
					SS(8,6) = -3*EIZ(MTYPE)/XL2
				ELSE
					SS(1,1) = E(MTYPE)*AREA(MTYPE)/XL
					SS(1,7) = -SS(1,1)
					SS(2,2) = 12 * EIZ(MTYPE)/XL3
					SS(2,6) =  6 * EIZ(MTYPE)/XL2
					SS(2,8) =-12 * EIZ(MTYPE)/XL3
					SS(2,12)=  6 * EIZ(MTYPE)/XL2
					SS(3,3) = 12 * EIY(MTYPE)/XL3
					SS(3,5) =- 6 * EIY(MTYPE)/XL2
					SS(3,9) =-12 * EIY(MTYPE)/XL3
					SS(3,11)= -6 * EIY(MTYPE)/XL2
					SS(4,4) =              GJ/XL
					SS(4,10)= -            GJ/XL
					SS(5,5) =  4 * EIY(MTYPE)/XL
					SS(5,9) =  6 * EIY(MTYPE)/XL2
					SS(5,11)=  2 * EIY(MTYPE)/XL
					SS(6,6) =  4 * EIZ(MTYPE)/XL
					SS(6,8) = -6 * EIZ(MTYPE)/XL2
					SS(6,12)=  2 * EIZ(MTYPE)/XL
					SS(7,7)=SS(1,1)
					SS(8,8) = 12 * EIZ(MTYPE)/XL3
					SS(8,12)= -6 * EIZ(MTYPE)/XL2
					SS(9,9) = 12 * EIY(MTYPE)/XL3
					SS(9,11)=  6 * EIY(MTYPE)/XL2
					SS(10,10) =            GJ/XL
					SS(11,11) = 4* EIY(MTYPE)/XL
					SS(12,12) = 4* EIZ(MTYPE)/XL
					
					DO I=2,12
						DO J=1,I-1
							SS(I,J)=SS(J,I)
						END DO
					END DO
				ENDIF
            ENDIF
            
		S=MATMUL(TRANSPOSE(TRANS),MATMUL(SS,TRANS))
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO
		
     ! Gravity
        REWIND ILOAD
        WRITE(ILOAD)R

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

      NNUM = 0
      FORCEA = 0.D0

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
                                           '  ELEMENT',13X,'AXISF',12X,'SHEARY',12X,'SHEARZ',12X,'  TORQUE',12X,' BMY',12X,'  BMZ')") NG
        MTYPE=MATP(N)

		
		GJ=G(MTYPE)*(EIY(MTYPE)+EIZ(MTYPE))/E(MTYPE)
	    
		CALL  TRANSPORT(TRANS,XL,XYZ(1,N),COR(1,N))! TRANSPORTION MATRIX
		XL2=XL**2
		XL3=XL**3
	  
		SS=0
		
		
		!IF(MODE.EQ.1)THEN! EULER BEAM ELEMENT

			IF(HGI(N).EQ.1)THEN
				IF(HGJ(N).EQ.1)THEN
					SS(1,1) =  E(MTYPE)*AREA(MTYPE)/XL
					SS(1,7) = -SS(1,1)
					SS(7,1) = -SS(1,1)
					SS(7,7) =  SS(1,1)
				ELSE
					SS(1,1) = E(MTYPE)*AREA(MTYPE)/XL
					SS(2,2) = 3*EIZ(MTYPE)/XL3
					SS(3,3) = 3*EIY(MTYPE)/XL3
					
					SS(1,7) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(2,8) = -3*EIZ(MTYPE)/XL3
					SS(3,9) = -3*EIY(MTYPE)/XL3
					
					SS(7,1) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(8,2) = -3*EIZ(MTYPE)/XL3
					SS(9,3) = -3*EIY(MTYPE)/XL3
					
					SS(8,8)=   3*EIZ(MTYPE)/XL3
					SS(9,9)=   3*EIY(MTYPE)/XL3
					
					SS(2,12) =  3*EIZ(MTYPE)/XL2
					SS(3,11) = -3*EIY(MTYPE)/XL2
					
					SS(12,2) =  3*EIZ(MTYPE)/XL2
					SS(11,3) = -3*EIY(MTYPE)/XL2
					
					SS(11,11) =  3*EIY(MTYPE)/XL
					SS(12,12) =  3*EIZ(MTYPE)/XL
					
					SS(11,9) =  3*EIY(MTYPE)/XL2
					SS(12,8) = -3*EIZ(MTYPE)/XL2
					
					SS(9,11) =  3*EIY(MTYPE)/XL2
					SS(8,12) = -3*EIZ(MTYPE)/XL2
				
				ENDIF
			ELSE
				IF(HGJ(N).EQ.1)THEN
					SS(1,1) = E(MTYPE)*AREA(MTYPE)/XL
					SS(2,2) = 3*EIZ(MTYPE)/XL3
					SS(3,3) = 3*EIY(MTYPE)/XL3
					
					SS(1,7) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(2,8) = -3*EIZ(MTYPE)/XL3
					SS(3,9) = -3*EIY(MTYPE)/XL3
					
					SS(7,1) = -E(MTYPE)*AREA(MTYPE)/XL
					SS(8,2) = -3*EIZ(MTYPE)/XL3
					SS(9,3) = -3*EIY(MTYPE)/XL3
					
					SS(8,8)=   3*EIZ(MTYPE)/XL3
					SS(9,9)=   3*EIY(MTYPE)/XL3
					
					SS(2,6) =  3*EIZ(MTYPE)/XL2
					SS(3,5) = -3*EIY(MTYPE)/XL2
					
					SS(6,2) =  3*EIZ(MTYPE)/XL2
					SS(5,3) = -3*EIY(MTYPE)/XL2
					
					SS(5,5) =  3*EIY(MTYPE)/XL
					SS(6,6) =  3*EIZ(MTYPE)/XL
					
					SS(5,9) =  3*EIY(MTYPE)/XL2
					SS(6,8) = -3*EIZ(MTYPE)/XL2
					
					SS(9,5) =  3*EIY(MTYPE)/XL2
					SS(8,6) = -3*EIZ(MTYPE)/XL2
				ELSE
					SS(1,1) = E(MTYPE)*AREA(MTYPE)/XL
					SS(1,7) = -SS(1,1)
					SS(2,2) = 12 * EIZ(MTYPE)/XL3
					SS(2,6) =  6 * EIZ(MTYPE)/XL2
					SS(2,8) =-12 * EIZ(MTYPE)/XL3
					SS(2,12)=  6 * EIZ(MTYPE)/XL2
					SS(3,3) = 12 * EIY(MTYPE)/XL3
					SS(3,5) =- 6 * EIY(MTYPE)/XL2
					SS(3,9) =-12 * EIY(MTYPE)/XL3
					SS(3,11)= -6 * EIY(MTYPE)/XL2
					SS(4,4) =              GJ/XL
					SS(4,10)= -            GJ/XL
					SS(5,5) =  4 * EIY(MTYPE)/XL
					SS(5,9) =  6 * EIY(MTYPE)/XL2
					SS(5,11)=  2 * EIY(MTYPE)/XL
					SS(6,6) =  4 * EIZ(MTYPE)/XL
					SS(6,8) = -6 * EIZ(MTYPE)/XL2
					SS(6,12)=  2 * EIZ(MTYPE)/XL
					SS(7,7)=SS(1,1)
					SS(8,8) = 12 * EIZ(MTYPE)/XL3
					SS(8,12)= -6 * EIZ(MTYPE)/XL2
					SS(9,9) = 12 * EIY(MTYPE)/XL3
					SS(9,11)=  6 * EIY(MTYPE)/XL2
					SS(10,10) =            GJ/XL
					SS(11,11) = 4* EIY(MTYPE)/XL
					SS(12,12) = 4* EIZ(MTYPE)/XL
					
					DO I=2,12
						DO J=1,I-1
							SS(I,J)=SS(J,I)
						END DO
					END DO
				ENDIF
            ENDIF
            
            DO L=1,6
                I=LM(L,N)
                IF(I.GT.0) DISP( L )=U(I)
                J=LM(L+6,N)
                IF(J.GT.0) DISP(L+6)=U(J)
            ENDDO
            D=MATMUL(TRANS,DISP)
            FORCE=MATMUL(SS,DISP)
            WRITE(IOUT,"(E13.6,6X,E13.6,6X,E13.6,6X,E13.6,6X,E13.6,6X,E13.6)")(FORCE(I),I=1,6)
			WRITE(IOUT,"(E13.6,6X,E13.6,6X,E13.6,6X,E13.6,6X,E13.6,6X,E13.6)")(FORCE(I),I=6,12)
			
			SIGMAX=FORCE(1)/AREA(MTYPE)
			SIGMAY=FORCE(6)*YMAX(MTYPE)*E(MTYPE)/EIZ(MTYPE)
			SIGMAZ=FORCE(5)*ZMAX(MTYPE)*E(MTYPE)/EIY(MTYPE)
			TAUYZ= FORCE(4)*SQRT(YMAX(MTYPE)**2+ZMAX(MTYPE)**2)*G(MTYPE)/GJ
			TAUXZ= FORCE(2)/AREA(MTYPE)
			TAUXY= FORCE(3)/AREA(MTYPE)
			VONMISES=SQRT((SIGMAX-SIGMAY)**2+(SIGMAY-SIGMAZ)**2+(SIGMAZ-SIGMAX)**2+6*(TAUXY**2+TAUYZ**2+TAUXZ**2))/SQRT(2.D0)
		
            FORCEA(:,ELID(1,N)) = FORCEA(:,ELID(1,N)) + (/SIGMAX,SIGMAY,SIGMAZ,TAUYZ,TAUXZ,TAUXY/)
            FORCEA(:,ELID(2,N)) = FORCEA(:,ELID(2,N)) + (/SIGMAX,SIGMAY,SIGMAZ,TAUYZ,TAUXZ,TAUXY/)
        
        NNUM(ELID(:,N)) = NNUM(ELID(:,N)) + 1
	
     END DO
     
    ! For postprocess
     DO I = 1,NUMNP
         IF(NNUM(I) .NE. 0) FORCEA(:,I) = FORCEA(:,I)/NNUM(I)
     END DO
    CALL POSTPROCESS(3,2,6,ID,X,Y,Z,U,ELID,FORCEA)

  ELSE 
     STOP "ERROR"
  END IF

END SUBROUTINE BEAM

SUBROUTINE TRANSPORT(TRANS,XL,XYZ,COR)
  !give element the transport matrix and length
   USE GLOBALS,ONLY: IOUT
   IMPLICIT NONE
   REAL(8),INTENT(OUT)::TRANS(12,12),XL
   REAL(8),INTENT(IN)::XYZ(6),COR(6)
   REAL(8)::COSX0X,COSX0Y,COSX0Z,COSY0X,COSY0Y,COSY0Z,COSZ0X,COSZ0Y,COSZ0Z
   REAL(8)::TRANS0(3,3),L
   INTEGER::J
            XL=SQRT((XYZ(1)-XYZ(4))**2+(XYZ(2)-XYZ(5))**2+(XYZ(3)-XYZ(6))**2) !单元长度
     		COSX0X=(-XYZ(1)+XYZ(4))/XL
     		COSX0Y=(-XYZ(2)+XYZ(5))/XL
     		COSX0Z=(-XYZ(3)+XYZ(6))/XL
            L=SQRT(COR(1)**2+COR(2)**2+COR(3)**2)
     		COSY0X=COR(1)/L
     		COSY0Y=COR(2)/L
     		COSY0Z=COR(3)/L
            L=SQRT(COR(4)**2+COR(5)**2+COR(6)**2)
     		COSZ0X=COR(4)/L
     		COSZ0Y=COR(5)/L
     		COSZ0Z=COR(6)/L
     		TRANS0(1,:)=(/COSX0X,COSX0Y,COSX0Z/)
     		TRANS0(2,:)=(/COSY0X,COSY0Y,COSY0Z/)
     		TRANS0(3,:)=(/COSZ0X,COSZ0Y,COSZ0Z/)
     	! TRANSPORT MATRIX	
     		TRANS=0
     		TRANS(1:3,1:3)=TRANS0
     		TRANS(4:6,4:6)=TRANS0
     		TRANS(7:9,7:9)=TRANS0
     		TRANS(10:12,10:12)=TRANS0
         
END SUBROUTINE TRANSPORT   		