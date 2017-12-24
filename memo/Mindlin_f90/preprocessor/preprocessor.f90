PROGRAM PREPROCESS

   IMPLICIT NONE
   
   integer, parameter :: IIN=10		! Unit used for input
   integer, parameter :: IOUT=11	! Unit used for output
   integer, parameter :: ILOAD=12   ! Unit storing Concentrated force
   
   integer :: NUMNP		! Total number of nodal points
   integer :: REDUCED_NUMNP
   
   integer :: IND		! Solution phase indicator
						!   1 - Read and generate element information
						!   2 - Assemble structure stiffness matrix
						!   3 - Stress calculations
   integer :: NPAR(10)	! Element group control data
						!   NPAR(1) - Element type
						!             1 : Truss element
                        !             2 : 4Q element
						!   NPAR(2) - Number of elements
						!   NPAR(3) - Number of different sets of material and 
						!             cross-sectional  constants
                        !   NPAR(4) - 在QUAD单元中判断是平面应力（0）还是平面应变（1）
                        !   NPAR(5) - 在8Q\9Q单元中判断是等参元（0）还是亚参元（1）
   integer :: NUMEG		! Total number of element groups, > 0
   
   character*200 :: HED	! Master heading information for use in labeling the output
   
   integer :: NLCASE=1    ! 荷载工况数
   
   character*200 :: LINE
   character*200 :: PARTLINE
   integer :: NUMLINE
   
   integer :: stat
   integer :: I, J, K, L
   
   integer,allocatable :: ID(:,:), N(:)
   real(8),allocatable :: X(:), Y(:), Z(:)
   
   integer,allocatable :: EID(:,:)
   real(8),allocatable :: EX(:), EY(:), EZ(:)
   real(8),allocatable :: DX(:), DY(:), DZ(:)
   real(8),allocatable :: X1(:), Y1(:), Z1(:), X2(:), Y2(:), Z2(:), THETA(:)     
   integer :: NUMNPE
   
   integer,allocatable :: SET(:)
   
   integer :: ASSEMBLYLINE, SETLINE, SURFACELINE, MATERIALLINE, BOUNDARYLINE, LOADLINE
   
   integer,allocatable :: NODELINE(:), ELEMENTLINE(:), ENDPARTLINE(:), NUMELENP(:), NUMELE(:), INNUMNP(:)
   character*30,allocatable :: PARTNAME(:), INSTANCENAME(:)
   character*10,allocatable :: ELEMENTTYPE(:)
   logical,allocatable :: ISMOP(:) !判断是否mesh on part
   
   character*30,allocatable :: TIESETNAME(:)
   character*30,allocatable :: MASTER(:), SLAVE(:)
   integer,allocatable :: COMBINENODE(:,:), TIESETLIST(:), KILLNODE(:)
   integer :: NUMTIE, NUMSET, NUMSETNODE, MINLIST
   
   character*30 :: BOUNDARYSETNAME, BSETINSTANCENAME
   integer,allocatable :: FIXEDDIRACTION(:)
   integer :: NUMDIRACTION, NUMBSETNODE, BSETLIST
   
   integer :: NLOAD  !本工况中集中载荷的个数
   character*30 :: LOADSETNAME, LSETINSTANCENAME
   real(8),allocatable :: LOAD(:)
   integer,allocatable :: LDIRACTION(:) 
   integer :: NUMLDIR, NUMLSETNODE, LSETLIST
   real(8) :: GRAVITY, GX, GY, GZ
   logical :: ISGRAVITY
   
   integer :: NUMMAT
   character*30 :: SECTIONSETNAME, MATERIALNAME
   
   CALL OPENFILES()
   
   READ(IIN,*)
   READ(IIN,"(A200)",iostat=stat) LINE
   I=INDEX(LINE,"Model name:")
   IF (I>15) THEN
       HED=LINE(14:I-2)
   ELSE
       HED="Job-1"
   END IF
   WRITE (IOUT,"(A200)") HED
   
   ASSEMBLYLINE=2
   DO WHILE (LINE(1:9)/="*Assembly")
       ASSEMBLYLINE=ASSEMBLYLINE+1
       READ(IIN,"(A200)",iostat=stat) LINE
   END DO
   
   READ(IIN,*)
   READ(IIN,"(A200)",iostat=stat) LINE
   
   !计算单元组数
   NUMEG=0
   DO WHILE (LINE(1:9)=="*Instance")
       NUMEG=NUMEG+1
       DO WHILE (LINE(1:13)/="*End Instance")
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       READ(IIN,*)
       READ(IIN,"(A200)",iostat=stat) LINE
   END DO
   
   ALLOCATE(NODELINE(NUMEG), ELEMENTLINE(NUMEG), ENDPARTLINE(NUMEG), NUMELENP(NUMEG), NUMELE(NUMEG), &
            PARTNAME(NUMEG), INSTANCENAME(NUMEG), ELEMENTTYPE(NUMEG), INNUMNP(NUMEG+1), ISMOP(NUMEG), &
            DX(NUMEG), DY(NUMEG), DZ(NUMEG), X1(NUMEG), Y1(NUMEG), Z1(NUMEG), X2(NUMEG), Y2(NUMEG), Z2(NUMEG), THETA(NUMEG))
   
   !
   REWIND IIN
   DO I=1,ASSEMBLYLINE+1
       READ(IIN,*)
   END DO
   
   !
   ISMOP=.TRUE.
   DO I=1,NUMEG
       READ(IIN,"(A200)",iostat=stat) LINE
       J=INDEX(LINE,"name=")
       K=INDEX(LINE,"part=")
       INSTANCENAME(I)=LINE(J+5:K-3)
       PARTNAME(I)=LINE(K+5:)
       DO WHILE (LINE(1:13)/="*End Instance")
           IF (LINE(1:5)=="*Node") ISMOP(I)=.FALSE.
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       READ(IIN,*)
   END DO
   
   INNUMNP(1)=0
   DO I=1,NUMEG
       REWIND IIN
       READ(IIN,"(A200)",iostat=stat) LINE
       IF (ISMOP(I)) THEN
           PARTLINE="*Part, name="//PARTNAME(I)
       ELSE
           PARTLINE="*Instance, name="//TRIM(INSTANCENAME(I))//", part="//TRIM(PARTNAME(I))
       END IF   
       NODELINE(I)=2
       DO WHILE (LINE/=PARTLINE)
           NODELINE(I)=NODELINE(I)+1
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       ELEMENTLINE(I)=NODELINE(I)-1
       DO WHILE (LINE(1:8)/="*Element")
           ELEMENTLINE(I)=ELEMENTLINE(I)+1
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       NUMELENP(I)=ELEMENTLINE(I)-NODELINE(I)-1
       INNUMNP(I+1)=INNUMNP(I)+NUMELENP(I)
       ELEMENTTYPE(I)=LINE(16:)
       READ(IIN,"(A200)",iostat=stat) LINE
       NUMELE(I)=0
       DO WHILE (LINE(1:1)/="*")
           NUMELE(I)=NUMELE(I)+1
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       ENDPARTLINE(I)=ELEMENTLINE(I)+NUMELE(I)+1
       IF (ISMOP(I)) THEN
           DO WHILE (TRIM(LINE)/="*End Part")
               ENDPARTLINE(I)=ENDPARTLINE(I)+1
               READ(IIN,"(A200)",iostat=stat) LINE
           END DO
       ELSE
           DO WHILE (TRIM(LINE)/="*End Instance")
               ENDPARTLINE(I)=ENDPARTLINE(I)+1
               READ(IIN,"(A200)",iostat=stat) LINE
           END DO
       END IF
   END DO
   
   !
   NUMNP=INNUMNP(NUMEG+1)
   REDUCED_NUMNP=NUMNP
   
   ALLOCATE(ID(6,NUMNP),X(NUMNP),Y(NUMNP),Z(NUMNP),N(NUMNP))
   ID=0
   
   DO I=1,NUMNP
       N(I)=I
   END DO
   
   DO I=1,NUMEG
       NUMNPE=NUMELENP(I)
       ALLOCATE(EX(NUMNPE),EY(NUMNPE),EZ(NUMNPE))
       REWIND IIN
       DO J=1,NODELINE(I)
           READ(IIN,*)
       END DO
       DO J=1,NUMNPE
           READ(IIN,"(I,3(',',F))",iostat=stat)K,EX(J),EY(J),EZ(J)
       END DO
       !
       IF (ISMOP(I)) THEN
           DO J=1,ASSEMBLYLINE-NODELINE(I)-NUMNPE+2
               READ(IIN,*)
           END DO
           DO J=1,I-1
               READ(IIN,"(A200)",iostat=stat) LINE
               DO WHILE (LINE(1:13)/="*End Instance")
                   READ(IIN,"(A200)",iostat=stat) LINE
               END DO
               READ(IIN,*)
               READ(IIN,*)
           END DO
           !
           READ(IIN,"(A200)",iostat=stat) LINE
           DO WHILE (LINE(1:13)/="*End Instance")
               CALL CHAR_REPEAT(',',LINE,J) 
               BACKSPACE(IIN,iostat=stat)
               IF (J==2) THEN
                   READ(IIN,"(F,2(',',F))",iostat=stat)DX(I), DY(I), DZ(I)
                   EX=EX+DX(I)
                   EY=EY+DY(I)
                   EZ=EZ+DZ(I)
               ELSEIF (J==6) THEN
                   READ(IIN,"(F,6(',',F))",iostat=stat)X1(I), Y1(I), Z1(I), X2(I), Y2(I), Z2(I), THETA(I)
                   DO K=1,NUMNPE
                       CALL REVOLVE(EX(K),EY(K),EZ(K))
                   END DO
               END IF
               READ(IIN,"(A200)",iostat=stat) LINE
           END DO
       END IF
       !
       DO J=1,NUMNPE
           X(INNUMNP(I)+J)=EX(J)
           Y(INNUMNP(I)+J)=EY(J)
           Z(INNUMNP(I)+J)=EZ(J)
       END DO
       DEALLOCATE(EX,EY,EZ)
       !
       IF ((ELEMENTTYPE(I)/="S4R").AND.(ELEMENTTYPE(I)/="B31").AND.(ELEMENTTYPE(I)/="CPS4R")) THEN
           DO J=1,NUMNPE
               ID(4:6,INNUMNP(I)+J)=1
           END DO
       ELSE IF(ELEMENTTYPE(I)=="S4R") THEN
           DO J=1,NUMNPE
               ID(6,INNUMNP(I)+J)=1
           END DO
       ELSE IF(ELEMENTTYPE(I)=="CPS4R") THEN
           DO J=1,NUMNPE
               ID(3:5,INNUMNP(I)+J)=1
           END DO
       END IF
   END DO
   
   NUMTIE=0
   DO WHILE (LINE(1:13)/="*End Assembly")
       IF (LINE(1:4)=="*Tie") THEN
           NUMTIE=NUMTIE+1
       END IF
       READ(IIN,"(A200)",iostat=stat) LINE
   END DO
   
   ALLOCATE(MASTER(NUMTIE), SLAVE(NUMTIE))
   
   DO I=1,NUMTIE*3-1
       BACKSPACE(IIN,iostat=stat)
   END DO
   
   DO I=1,NUMTIE
       READ(IIN,"(A200)",iostat=stat) LINE
       J=INDEX(LINE,",")
       SLAVE(I)=LINE(1:J-1)
       MASTER(I)=LINE(J+2:)
       READ(IIN,*)
       READ(IIN,*)
   END DO
   
   REWIND IIN
   DO J=1,ASSEMBLYLINE
       READ(IIN,*)
   END DO
   
   READ(IIN,"(A200)",iostat=stat) LINE
   SETLINE=ASSEMBLYLINE+1
   DO WHILE ((LINE(1:5)/="*Nset").AND.(LINE(1:6)/="*Elset"))
       SETLINE=SETLINE+1
       READ(IIN,"(A200)",iostat=stat) LINE
   END DO
   
   SURFACELINE=SETLINE
   DO WHILE (LINE(1:8)/="*Surface")
       SURFACELINE=SURFACELINE+1
       READ(IIN,"(A200)",iostat=stat) LINE
   END DO
   
   DO I=1,NUMTIE
       REWIND IIN
       DO J=1,SURFACELINE-1
           READ(IIN,*)
       END DO
       J=0
       DO WHILE (J==0)
           READ(IIN,"(A200)",iostat=stat) LINE
           J=INDEX(LINE,TRIM(MASTER(I)))
       END DO
       READ(IIN,"(A200)",iostat=stat) LINE
       J=INDEX(LINE,",")
       MASTER(I)=LINE(1:J-1)
       J=0
       DO WHILE (J==0)
           READ(IIN,"(A200)",iostat=stat) LINE
           J=INDEX(LINE,TRIM(SLAVE(I)))
       END DO
       READ(IIN,"(A200)",iostat=stat) LINE
       J=INDEX(LINE,",")
       SLAVE(I)=LINE(1:J-1)
   END DO
   
   DO I=1,NUMTIE
       REWIND IIN
       DO J=1,SETLINE-1
           READ(IIN,*)
       END DO
       READ(IIN,"(A200)",iostat=stat) LINE
       NUMSET=0
       DO WHILE (LINE(1:8)/="*Surface")
           IF ((LINE(1:(13+LEN_TRIM(MASTER(I))))=="*Nset, nset="//TRIM(MASTER(I))//",").OR.(LINE(1:(13+LEN_TRIM(SLAVE(I))))=="*Nset, nset="//TRIM(SLAVE(I))//",")) THEN
               NUMSET=NUMSET+1
           END IF
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       !
       ALLOCATE(TIESETNAME(NUMSET), TIESETLIST(NUMSET))
       !
       DO J=1,SURFACELINE-SETLINE+1
           BACKSPACE(IIN,iostat=stat)
       END DO
       READ(IIN,"(A200)",iostat=stat) LINE
       DO WHILE ((LINE(1:(13+LEN_TRIM(MASTER(I))))/="*Nset, nset="//TRIM(MASTER(I))//",").AND.(LINE(1:(13+LEN_TRIM(SLAVE(I))))/="*Nset, nset="//TRIM(SLAVE(I))//","))
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       CALL SETIN(NUMSETNODE,J)
       ALLOCATE(COMBINENODE(NUMSETNODE,NUMSET), KILLNODE(NUMSETNODE*(NUMSET-1)))
       DEALLOCATE(SET)
       
       DO K=1,J+1
           BACKSPACE(IIN,iostat=stat)
       END DO
                 
       READ(IIN,"(A200)",iostat=stat) LINE
       J=0
       DO WHILE (LINE(1:8)/="*Surface")
           IF ((LINE(1:(13+LEN_TRIM(MASTER(I))))=="*Nset, nset="//TRIM(MASTER(I))//",").OR.&
               (LINE(1:(13+LEN_TRIM(SLAVE(I))))=="*Nset, nset="//TRIM(SLAVE(I))//",")) THEN
               J=J+1
               K=INDEX(LINE,"instance=")
               TIESETNAME(J)=LINE(K+9:)
               K=1
               DO WHILE (INSTANCENAME(K)/=TIESETNAME(J))
                   K=K+1
               END DO
               TIESETLIST(J)=K
               CALL SETIN(NUMSETNODE,K)
               IF (J==1) THEN
                   COMBINENODE(:,1)=SET(:)
               ELSE
                   DO K=1,NUMSETNODE
                       DO L=1,NUMSETNODE
                           IF ((ABS(X(N(INNUMNP(TIESETLIST(J))+SET(L)))-X(N(INNUMNP(TIESETLIST(1))+COMBINENODE(K,1))))<1.D-10).AND.(ABS(Y(N(INNUMNP(TIESETLIST(J))+SET(L)))-Y(N(INNUMNP(TIESETLIST(1))+COMBINENODE(K,1))))<1.D-10).AND.(ABS(Z(N(INNUMNP(TIESETLIST(J))+SET(L)))-Z(N(INNUMNP(TIESETLIST(1))+COMBINENODE(K,1))))<1.D-10)) THEN
                               COMBINENODE(K,J)=SET(L)
                           END IF
                       END DO
                   END DO
               END IF
               DEALLOCATE(SET)
           END IF
           READ(IIN,"(A200)",iostat=stat) LINE
       END DO
       !
       MINLIST=MINLOC(TIESETLIST,1)
       DO K=1,NUMSETNODE
           ID(4:6,N(INNUMNP(TIESETLIST(MINLIST))+COMBINENODE(K,MINLIST)))=1
       END DO
       DO J=1,NUMSET
           IF (J<MINLIST) THEN
               DO K=1,NUMSETNODE
                   KILLNODE((J-1)*NUMSETNODE+K)=INNUMNP(TIESETLIST(J))+COMBINENODE(K,J)
               END DO
           ELSE IF (J>MINLIST) THEN
               DO K=1,NUMSETNODE
                   KILLNODE((J-2)*NUMSETNODE+K)=INNUMNP(TIESETLIST(J))+COMBINENODE(K,J)
               END DO
           END IF
       END DO  
       CALL DELETE_NODE(KILLNODE, NUMSETNODE*(NUMSET-1))
       DO J=1,NUMSET
           IF (J/=MINLIST) THEN
               DO K=1,NUMSETNODE
                   N(INNUMNP(TIESETLIST(J))+COMBINENODE(K,J))=N(INNUMNP(TIESETLIST(MINLIST))+COMBINENODE(K,MINLIST))
               END DO
           END IF
       END DO
       
       DEALLOCATE(TIESETNAME, TIESETLIST, COMBINENODE, KILLNODE)
   END DO
   
   MATERIALLINE=SURFACELINE
   DO WHILE (TRIM(LINE)/="** MATERIALS")
       MATERIALLINE=MATERIALLINE+1
       READ(IIN,"(A200)",iostat=stat) LINE
   END DO 
   
   BOUNDARYLINE=MATERIALLINE
   DO WHILE (TRIM(LINE)/="** BOUNDARY CONDITIONS")
       BOUNDARYLINE=BOUNDARYLINE+1
       READ(IIN,"(A200)",iostat=stat) LINE
   END DO
      
   NUMLINE=BOUNDARYLINE
   DO WHILE (TRIM(LINE)/="** LOADS")
       NUMLINE=NUMLINE+1
       READ(IIN,"(A200)",iostat=stat) LINE
       IF (TRIM(LINE)=="*Boundary") THEN
           READ(IIN,"(A200)",iostat=stat) LINE
           NUMDIRACTION=0
           DO WHILE (LINE(1:2)/="**")
               NUMDIRACTION=NUMDIRACTION+1
               READ(IIN,"(A200)",iostat=stat) LINE
           END DO
           NUMLINE=NUMLINE+NUMDIRACTION+1
           ALLOCATE(FIXEDDIRACTION(NUMDIRACTION))
           DO I=1,NUMDIRACTION+1
               BACKSPACE(IIN,iostat=stat)
           END DO
           DO I=1,NUMDIRACTION
               READ(IIN,"(A200)",iostat=stat) LINE
               J=INDEX(LINE,",")
               BOUNDARYSETNAME=LINE(1:J-1)
               READ(LINE(J+2:),"(I)")FIXEDDIRACTION(I)
           END DO
           NUMLINE=NUMLINE-1
           !
           DO I=1,NUMLINE-SETLINE+1
               BACKSPACE(IIN,iostat=stat)
           END DO
           READ(IIN,"(A200)",iostat=stat) LINE
           DO WHILE (LINE(1:8)/="*Surface")
               IF (LINE(1:12+LEN_TRIM(BOUNDARYSETNAME))=="*Nset, nset="//TRIM(BOUNDARYSETNAME)) THEN
                   J=INDEX(LINE,"instance=")
                   BSETINSTANCENAME=LINE(J+9:)
                   J=1
                   DO WHILE (INSTANCENAME(J)/=BSETINSTANCENAME)
                       J=J+1
                   END DO
                   BSETLIST=J
                   CALL SETIN(NUMBSETNODE,K)
                   DO J=1,NUMBSETNODE
                       DO K=1,NUMDIRACTION
                           ID(FIXEDDIRACTION(K),N(INNUMNP(BSETLIST)+SET(J)))=1
                       END DO
                   END DO
                   DEALLOCATE(SET)
               END IF
               READ(IIN,"(A200)",iostat=stat) LINE
           END DO
           DO I=1,NUMLINE-SURFACELINE-1
               READ(IIN,*)
           END DO
           READ(IIN,"(A200)",iostat=stat) LINE
           DEALLOCATE(FIXEDDIRACTION)
       END IF
   END DO
   LOADLINE=NUMLINE
   
   !
   NLOAD=0
   ISGRAVITY=.FALSE.
   DO WHILE (TRIM(LINE)/="** OUTPUT REQUESTS")
       NUMLINE=NUMLINE+1
       READ(IIN,"(A200)",iostat=stat) LINE
       IF (TRIM(LINE)=="*Dload") THEN
           ISGRAVITY=.TRUE.
           NUMLINE=NUMLINE+1
           READ(IIN,"(A200)",iostat=stat) LINE
           I=INDEX(LINE,",")
           I=INDEX(LINE(I+1:),",")
           READ(LINE(I+2:),"(F,3(',',F))") GRAVITY, GX, GY, GZ
       ELSE IF (TRIM(LINE)=="*Cload") THEN
           NLOAD=NLOAD+1
       END IF    
   END DO
   
   IF (NLOAD>0) THEN
       NLOAD=0
       DO I=1,NUMLINE-LOADLINE+1
           BACKSPACE(IIN,iostat=stat)
       END DO
       READ(IIN,"(A200)",iostat=stat) LINE
       NUMLINE=LOADLINE
       DO WHILE (TRIM(LINE)/="** OUTPUT REQUESTS")
           NUMLINE=NUMLINE+1
           READ(IIN,"(A200)",iostat=stat) LINE
           IF (TRIM(LINE)=="*Cload") THEN
               READ(IIN,"(A200)",iostat=stat) LINE
               NUMLDIR=0
               DO WHILE (LINE(1:2)/="**")
                   NUMLDIR=NUMLDIR+1
                   READ(IIN,"(A200)",iostat=stat) LINE
               END DO
               NUMLINE=NUMLINE+NUMLDIR+1
               ALLOCATE(LOAD(NUMLDIR), LDIRACTION(NUMLDIR))
               DO I=1,NUMLDIR+1
                   BACKSPACE(IIN,iostat=stat)
               END DO
               DO I=1,NUMLDIR
                   READ(IIN,"(A200)",iostat=stat) LINE
                   J=INDEX(LINE,",")
                   LOADSETNAME=LINE(1:J-1)
                   READ(LINE(J+2:),"(I,',',F)")LDIRACTION(I),LOAD(I)
               END DO
               NUMLINE=NUMLINE-1
               !
               DO I=1,NUMLINE-SETLINE+1
                   BACKSPACE(IIN,iostat=stat)
               END DO
               READ(IIN,"(A200)",iostat=stat) LINE
               DO WHILE (LINE(1:8)/="*Surface")
                   IF (LINE(1:12+LEN_TRIM(LOADSETNAME))=="*Nset, nset="//TRIM(LOADSETNAME)) THEN
                       J=INDEX(LINE,"instance=")
                       LSETINSTANCENAME=LINE(J+9:)
                       J=1
                       DO WHILE (INSTANCENAME(J)/=BSETINSTANCENAME)
                           J=J+1
                       END DO
                       BSETLIST=J
                       CALL SETIN(NUMLSETNODE,K)
                       NLOAD=NLOAD+NUMLSETNODE*NUMLDIR
                       DO J=1,NUMLSETNODE
                           DO K=1,NUMLDIR
                               WRITE(ILOAD,"(I10,I10,F20.12)")SET(J),LDIRACTION(K),LOAD(K)
                           END DO
                       END DO
                       DEALLOCATE(SET)
                   END IF
               END DO
               DO I=1,NUMLINE-SURFACELINE-1
                   READ(IIN,*)
               END DO
               READ(IIN,"(A200)",iostat=stat) LINE
               DEALLOCATE(LOAD, LDIRACTION)
           END IF
       END DO
   END IF
   
   IF (ISGRAVITY) THEN
       WRITE (IOUT,"(4(I10),3F20.12)",iostat=stat) REDUCED_NUMNP, NUMEG, NLCASE, 1, GRAVITY*GX, GRAVITY*GY, GRAVITY*GZ
   ELSE
       WRITE (IOUT,"(4(I10))",iostat=stat) REDUCED_NUMNP, NUMEG, NLCASE, 1
   END IF
   DO I=1,REDUCED_NUMNP
       WRITE (IOUT,"(7I10,3F20.12)",iostat=stat) I, (ID(J,I), J=1,6), X(I), Y(I), Z(I)
   END DO
   
   WRITE (IOUT,"(2I10)",iostat=stat)1,NLOAD
   IF (NLOAD>0) THEN
       REWIND ILOAD
       DO I=1,NLOAD
           READ (ILOAD,"(A200)",iostat=stat) LINE
           WRITE (IOUT,"(A200)",iostat=stat) LINE
       END DO
   ELSE
       WRITE (IOUT,*)
   END IF
   
   DO I=1,NUMEG
       REWIND IIN
       DO J=1,ELEMENTLINE(I)+NUMELE(I)
           READ(IIN,*)
       END DO
       READ(IIN,"(A200)",iostat=stat) LINE
       NUMLINE=ELEMENTLINE(I)+NUMELE(I)+1
       NUMMAT=0
       IF (ISMOP(I)) THEN
           DO WHILE(TRIM(LINE)/="*End Part")
               IF (LINE(1:10)=="** Section") NUMMAT=NUMMAT+1
               READ(IIN,"(A200)",iostat=stat) LINE
               NUMLINE=NUMLINE+1
           END DO
       ELSE
           DO WHILE(TRIM(LINE)/="*End Instance")
               IF (LINE(1:10)=="** Section") NUMMAT=NUMMAT+1
               READ(IIN,"(A200)",iostat=stat) LINE
               NUMLINE=NUMLINE+1
           END DO
       END IF
       DO J=1,NUMLINE-ELEMENTLINE(I)-NUMELE(I)
           BACKSPACE(IIN,iostat=stat)
       END DO
       READ(IIN,"(A200)",iostat=stat) LINE
       NUMLINE=ELEMENTLINE(I)+NUMELE(I)+1
       NPAR(2)=NUMELE(I)
       NPAR(3)=NUMMAT
       SELECT CASE (TRIM(ELEMENTTYPE(I)))
       CASE ("T3D2")
           NPAR(1)=1
           CALL TRUSS
       CASE ("CPS4R")
           NPAR(1)=2
           CALL QUADS
       !CASE ("")
       !    NPAR(1)=3
       !    CALL TRIANGLE
       CASE ("C3D8R")
           NPAR(1)=4
           CALL HEXA8
       CASE ("B31")
           NPAR(1)=5
           CALL BEAMS
       !CASE ("")
       !    NPAR(1)=6
       !    CALL PLATE
       CASE ("S4R")
           NPAR(1)=7
           CALL SHELL
       !CASE ("")
       !    NPAR(1)=8
       !    CALL EIGHT_QUAD
       !CASE ("")
       !    NPAR(1)=9
       !    CALL NINE_QUAD
       !CASE ("")
       !    NPAR(1)=10
       !    CALL SIXTRIANGLE
       !CASE ("")
       !    NPAR(1)=11
       !    CALL TETRAHEDRON
       !CASE ("")
       !    NPAR(1)=12
       !    CALL TBEAMS
       !CASE ("")
       !    NPAR(1)=13
       !    CALL INFINITE
       CASE DEFAULT
           PRINT *, "*** STOP *** ABAQUS INP FILE DOES NOT EXIST !"
           STOP
       END SELECT
       
   END DO
   
   DEALLOCATE(NODELINE, ELEMENTLINE, ENDPARTLINE, NUMELENP, NUMELE, PARTNAME, INSTANCENAME, ELEMENTTYPE, INNUMNP, ISMOP, &
              DX, DY, DZ, X1, Y1, Z1, X2, Y2, Z2, THETA)
   DEALLOCATE(ID,X,Y,Z,N)
   DEALLOCATE(MASTER, SLAVE)
   CALL CLOSEFILES()

CONTAINS
    
SUBROUTINE OPENFILES()
! . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                   .
! .   Open input data file and results output file    .
! . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  LOGICAL :: EX
  CHARACTER*80 FileInp

  if(COMMAND_ARGUMENT_COUNT().ne.1) then
     stop 'Usage: InputFileName'
  else
     call GET_COMMAND_ARGUMENT(1,FileInp)
  end if

  INQUIRE(FILE = FileInp, EXIST = EX)
  IF (.NOT. EX) THEN
     PRINT *, "*** STOP *** ABAQUS INP FILE DOES NOT EXIST !"
     STOP
  END IF

  OPEN(IIN   , FILE = FileInp,  &
               FORM = "FORMATTED", &
               STATUS = "OLD", &
               ACCESS = "SEQUENTIAL", &
               ACTION = "READ")
  OPEN(IOUT  , FILE = "STAP90.in", &
               FORM = "FORMATTED", &
               STATUS = "REPLACE", &
               ACCESS = "SEQUENTIAL", &
               ACTION = "WRITE")
  OPEN(ILOAD , FILE = "LOAD.TMP", &
               FORM = "FORMATTED", &
               STATUS = "REPLACE", &
               ACCESS = "SEQUENTIAL", &
               ACTION = "READWRITE")
  
END SUBROUTINE OPENFILES
    
SUBROUTINE CLOSEFILES()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   Close all data files                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  CLOSE(IIN)
  CLOSE(IOUT)
  
END SUBROUTINE CLOSEFILES

SUBROUTINE CHAR_REPEAT(CHAR,A,NUM)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   计算一行中某个字符出现的次数                                    .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  CHARACTER*1 :: CHAR
  CHARACTER*200 :: A
  INTEGER :: I,NUM
  
  NUM=0
  DO I=1,200
      IF (A(I:I)==CHAR) THEN
          NUM=NUM+1
      END IF
  END DO
  
END SUBROUTINE CHAR_REPEAT

SUBROUTINE REVOLVE(X,Y,Z)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   (X,Y,Z)绕某个轴旋转THETA角                                      .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  real(8) :: X, Y, Z
  real(8) :: A, B, C, NX, NY, NZ, T, L, NL
  
  A=X2(I)-X1(I)
  B=Y2(I)-Y1(I)
  C=Z2(I)-Z1(I)
  X=X-X1(I)
  Y=Y-Y1(I)
  Z=Z-Z1(I)
  T=(A*X+B*Y+C*Z)/(A**2+B**2+C**2)
  X=X-A*T
  Y=Y-B*T
  Z=Z-C*T
  IF ((X**2+Y**2+Z**2)>=EPSILON(1.D0)) THEN
      NX=B*Z-C*Y
      NY=C*X-A*Z
      NZ=A*Y-B*X
      L=SQRT(X**2+Y**2+Z**2)
      NL=SQRT(NX**2+NY**2+NZ**2)
      NX=NX*L/NL
      NY=NY*L/NL
      NZ=NZ*L/NL
      X=X*COSD(THETA(I))+NX*SIND(THETA(I))
      Y=Y*COSD(THETA(I))+NY*SIND(THETA(I))
      Z=Z*COSD(THETA(I))+NZ*SIND(THETA(I))
  END IF
  X=X+X1(I)+A*T
  Y=Y+Y1(I)+B*T
  Z=Z+Z1(I)+C*T
  
END SUBROUTINE REVOLVE

SUBROUTINE SETIN(SIZE,HANG)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取点集数据                                                    .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  INTEGER :: SIZE, I, J, HANG, A, B, DX
  INTEGER,ALLOCATABLE :: NUM(:), SUMNUM(:)
  CHARACTER*3 :: CHAR_NUM_1
  
  BACKSPACE(IIN,iostat=stat)
  READ(IIN,"(A200)",iostat=stat) LINE
  I=INDEX(LINE, "generate")
  IF (I>0) THEN
      READ(IIN,"(I,',',I,',',I)",iostat=stat)A, B, DX
      SIZE=(B-A)/DX+1
      ALLOCATE(SET(SIZE))
      DO I=1,SIZE
          SET(I)=A+(I-1)*DX
      END DO
  ELSE
      READ(IIN,"(A200)",iostat=stat) LINE
      HANG=0
      DO WHILE (LINE(1:1)/="*")
          HANG=HANG+1
          READ(IIN,"(A200)",iostat=stat) LINE
      END DO
  
      ALLOCATE(NUM(HANG),SUMNUM(HANG+1))
      DO I=1,HANG+1
          BACKSPACE(IIN,iostat=stat)
      END DO
      SIZE=0
      SUMNUM(1)=0
      DO I=1,HANG
          READ(IIN,"(A200)",iostat=stat) LINE
          CALL CHAR_REPEAT(",",LINE,NUM(I))
          NUM(I)=NUM(I)+1
          IF (NUM(I)==2) THEN
              J=LEN_TRIM(LINE)
              IF (LINE(J:J)==",") THEN
                  NUM(I)=NUM(I)-1
              END IF
          END IF
          SUMNUM(I+1)=SUMNUM(I)+NUM(I)
      END DO
  
      SIZE=SUMNUM(HANG+1)
      ALLOCATE(SET(SIZE))
      DO I=1,HANG
          BACKSPACE(IIN,iostat=stat)
      END DO
      DO I=1,HANG
          IF (NUM(I)==1) THEN
              READ(IIN,"(I,',')",iostat=stat) (SET(SUMNUM(I)+J),J=1,NUM(I))
          ELSE
              WRITE(CHAR_NUM_1,"(I3)") NUM(I)-1
              CHAR_NUM_1=ADJUSTL(CHAR_NUM_1)
              READ(IIN,"(I,"//CHAR_NUM_1//"(',',I))",iostat=stat) (SET(SUMNUM(I)+J),J=1,NUM(I))
          END IF
      END DO
      DEALLOCATE(NUM,SUMNUM)
  END IF
  
END SUBROUTINE SETIN

SUBROUTINE DELETE_NODE(NODE, NUM)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   删除多余节点                                                    .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  INTEGER :: NODE(NUM), NUM, I, J, A, TOU, WEI  
  
  !
  DO I=1,NUM
      NODE(I)=N(NODE(I))
  END DO
  
  ! 排序
  DO I=1,NUM-1
      DO J=I+1,NUM
          IF (NODE(I)>NODE(J)) THEN
              A=NODE(I)
              NODE(I)=NODE(J)
              NODE(J)=A
          END IF
      END DO
  END DO
  
  DO I=1,NUM
      TOU=NODE(I)-I+1
      IF (I==NUM) THEN
          WEI=REDUCED_NUMNP-NUM
      ELSE
          WEI=NODE(I+1)-I-1
      END IF
      DO J=TOU,WEI
          X(J)=X(J+I)
          Y(J)=Y(J+I)
          Z(J)=Z(J+I)
          ID(1:6,J)=ID(1:6,J+I)
      END DO 
      DO J=1,NUMNP
          IF (I==NUM) THEN
              IF (N(J)>NODE(I)) N(J)=N(J)-I
          ELSE
              IF ((N(J)>NODE(I)).AND.(N(J)<NODE(I+1))) N(J)=N(J)-I
          END IF
      END DO
  END DO
  
  REDUCED_NUMNP=REDUCED_NUMNP-NUM
  
END SUBROUTINE DELETE_NODE

SUBROUTINE MATERIALIN(DENSITY, E, POISSON)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取材料属性                                                    .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  
  real(8) :: DENSITY, E, POISSON
  integer :: J, K, IE
  
  READ(IIN,"(A200)",iostat=stat) LINE
  DO WHILE(TRIM(LINE)/="** BOUNDARY CONDITIONS")
      IF (TRIM(LINE)=="*Material, name="//TRIM(MATERIALNAME)) THEN
          DO WHILE (LINE(1:2)/="**")
              IF (TRIM(LINE)=="*Density") THEN
                  READ(IIN,"(F)",iostat=stat) DENSITY
              ELSE IF (TRIM(LINE)=="*Elastic") THEN
                  READ(IIN,"(A200)",iostat=stat) LINE
                  J=INDEX(LINE,",")
                  K=INDEX(LINE,".")
                  BACKSPACE(IIN,iostat=stat)
                  IF (K<J) THEN
                      READ(IIN,"(F,',',F)",iostat=stat) E,POISSON
                  ELSE 
                      READ(IIN,"(E,',',F)",iostat=stat) IE,POISSON
                      E=REAL(IE)
                  END IF
              END IF
              READ(IIN,"(A200)",iostat=stat) LINE
          END DO
      END IF
      READ(IIN,"(A200)",iostat=stat) LINE
  END DO
  
END SUBROUTINE MATERIALIN

SUBROUTINE ELEMENTOUT(ELEN,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取并输出单元数据                                              .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  integer :: ELEN, MATP(NUMELE(I)), IEN(ELEN)
  integer :: J, K, L, NUMELEI
  CHARACTER*3 :: CHAR_ELEN
  
  NUMELEI=NUMELE(I)
  WRITE (CHAR_ELEN,"(I3)") ELEN
  CHAR_ELEN=ADJUSTL(CHAR_ELEN)
  DO J=1,NUMELEI
      READ (IIN,"(I,"//CHAR_ELEN//"(',',I))") L, (IEN(K),K=1,ELEN)
      WRITE (IOUT,"(I10,I10,"//CHAR_ELEN//"I10))") J, (N(INNUMNP(I)+IEN(K)),K=1,ELEN), MATP(J)
  END DO
  
END SUBROUTINE ELEMENTOUT

SUBROUTINE TRUSS()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取并输出TRUSS单元数据                                         .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  integer :: MATP(NUMELE(I))
  real(8) :: E(NUMMAT), AREA(NUMMAT), DENSITY(NUMMAT), POISSON   ! 
  integer :: J, K, IMAT, NUMSETNODE
  character*20 :: ENDMARK
  
  WRITE (IOUT,"(3I10)") (NPAR(J),J=1,3)
  IMAT=0
  IF (ISMOP(I)) THEN
      ENDMARK="*End Part"
  ELSE
      ENDMARK="*End Instance"
  END IF
  DO WHILE(TRIM(LINE)/=ENDMARK)
      IF (LINE(1:10)=="** Section") THEN
          IMAT=IMAT+1
          READ(IIN,"(A200)",iostat=stat) LINE
          J=INDEX(LINE,"elset=")
          K=INDEX(LINE,"material=")
          SECTIONSETNAME=LINE(J+6:K-3)
          MATERIALNAME=LINE(K+9:)
          READ(IIN,"(F)",iostat=stat) AREA(IMAT)
          NUMLINE=NUMLINE+2
          DO J=1,NUMLINE-ELEMENTLINE(I)-NUMELE(I)
              BACKSPACE(IIN,iostat=stat)
          END DO
          READ(IIN,"(A200)",iostat=stat) LINE
          DO WHILE(TRIM(LINE)/=ENDMARK)
              IF (LINE(1:14+LEN_TRIM(SECTIONSETNAME))=="*Elset, elset="//TRIM(SECTIONSETNAME)) THEN
                  CALL SETIN(NUMSETNODE,J)
                  DO J=1,NUMSETNODE
                      MATP(SET(J))=IMAT
                  END DO
                  DEALLOCATE(SET)
              END IF
              READ(IIN,"(A200)",iostat=stat) LINE
          END DO
          DO J=1,MATERIALLINE-ENDPARTLINE(I)
              READ(IIN,*)
          END DO
          CALL MATERIALIN(DENSITY(IMAT), E(IMAT), POISSON)
          DO J=1,BOUNDARYLINE-NUMLINE
              BACKSPACE(IIN,iostat=stat)
          END DO
      END IF
      READ(IIN,"(A200)",iostat=stat) LINE
      NUMLINE=NUMLINE+1
  END DO
  
  DO J=1,NUMMAT
      WRITE (IOUT,"(I10,E20.12,5F20.12)") J, E(J), AREA(J), DENSITY(J)
  END DO
  
  DO J=1,ENDPARTLINE(I)-ELEMENTLINE(I)
      BACKSPACE(IIN,iostat=stat)
  END DO
  CALL ELEMENTOUT(2,MATP)
  
END SUBROUTINE TRUSS

SUBROUTINE HEXA8()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取并输出8H单元数据                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  integer :: MATP(NUMELE(I))
  real(8) :: E(NUMMAT), POISSON(NUMMAT), DENSITY(NUMMAT) 
  integer :: J, K, IMAT, NUMSETNODE
  character*20 :: ENDMARK
  
  WRITE (IOUT,"(3I10)") (NPAR(J),J=1,3)
  IMAT=0
  IF (ISMOP(I)) THEN
      ENDMARK="*End Part"
  ELSE
      ENDMARK="*End Instance"
  END IF
  DO WHILE(TRIM(LINE)/=ENDMARK)
      IF (LINE(1:10)=="** Section") THEN
          IMAT=IMAT+1
          READ(IIN,"(A200)",iostat=stat) LINE
          J=INDEX(LINE,"elset=")
          K=INDEX(LINE,"material=")
          SECTIONSETNAME=LINE(J+6:K-3)
          MATERIALNAME=LINE(K+9:)
          NUMLINE=NUMLINE+1
          DO J=1,NUMLINE-ELEMENTLINE(I)-NUMELE(I)
              BACKSPACE(IIN,iostat=stat)
          END DO
          READ(IIN,"(A200)",iostat=stat) LINE
          DO WHILE(TRIM(LINE)/=ENDMARK)
              IF (LINE(1:14+LEN_TRIM(SECTIONSETNAME))=="*Elset, elset="//TRIM(SECTIONSETNAME)) THEN
                  CALL SETIN(NUMSETNODE,J)
                  DO J=1,NUMSETNODE
                      MATP(SET(J))=IMAT
                  END DO
                  DEALLOCATE(SET)
              END IF
              READ(IIN,"(A200)",iostat=stat) LINE
          END DO
          DO J=1,MATERIALLINE-ENDPARTLINE(I)
              READ(IIN,*)
          END DO
          CALL MATERIALIN(DENSITY(IMAT), E(IMAT), POISSON(IMAT))
          DO J=1,BOUNDARYLINE-NUMLINE
              BACKSPACE(IIN,iostat=stat)
          END DO
      END IF
      READ(IIN,"(A200)",iostat=stat) LINE
      NUMLINE=NUMLINE+1
  END DO
  
  DO J=1,NUMMAT
      WRITE (IOUT,"(I10,E20.12,5F20.12)") J, E(J), POISSON(J), DENSITY(J)
  END DO
  
  DO J=1,ENDPARTLINE(I)-ELEMENTLINE(I)
      BACKSPACE(IIN,iostat=stat)
  END DO
  CALL ELEMENTOUT(8,MATP)
  
END SUBROUTINE HEXA8

SUBROUTINE BEAMS()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取并输出欧拉梁单元数据                                        .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  integer :: MATP(NUMELE(I)), INDX(NUMMAT)
  real(8) :: E(NUMMAT), AREA(NUMMAT), G(NUMMAT), EIY(NUMMAT), EIZ(NUMMAT), DENSITY(NUMMAT), YMAX(NUMMAT), ZMAX(NUMMAT)
  real(8) :: POISSON(NUMMAT), A(NUMMAT), B(NUMMAT), T1(NUMMAT), T2(NUMMAT), T3(NUMMAT), T4(NUMMAT), APPN1(3,NUMMAT), VN1(6)
  real(8) :: IY(NUMMAT), IZ(NUMMAT)
  integer :: J, K, L, IMAT, NUMSETNODE
  character*20 :: ENDMARK
  
  WRITE (IOUT,"(3I10)") (NPAR(J),J=1,3)
  IMAT=0
  IF (ISMOP(I)) THEN
      ENDMARK="*End Part"
  ELSE
      ENDMARK="*End Instance"
  END IF
  DO WHILE(TRIM(LINE)/=ENDMARK)
      IF (LINE(1:10)=="** Section") THEN
          IMAT=IMAT+1
          READ(IIN,"(A200)",iostat=stat) LINE
          J=INDEX(LINE,"elset=")
          K=INDEX(LINE,"material=")
          L=INDEX(LINE(K:),",")
          SECTIONSETNAME=LINE(J+6:K-3)
          MATERIALNAME=LINE(K+9:K+L-2)
          INDX(IMAT)=1
          READ(IIN,"(F,5(',',F))",iostat=stat) A(IMAT), B(IMAT), T1(IMAT), T2(IMAT), T3(IMAT), T4(IMAT)
          CALL BEAM_SECTION_BOX(AREA(IMAT), IY(IMAT), IZ(IMAT), YMAX(IMAT), ZMAX(IMAT), A(IMAT), B(IMAT), T1(IMAT), T2(IMAT), T3(IMAT), T4(IMAT))
          READ(IIN,"(F2.0,2(',',F2.0))",iostat=stat) (APPN1(J,IMAT),J=1,3)
          VN1=0
          VN1(4:6)=APPN1(:,IMAT)+VN1(1:3)
          CALL REVOLVE(VN1(1),VN1(2),VN1(3))
          CALL REVOLVE(VN1(4),VN1(5),VN1(6))
          APPN1(:,IMAT)=VN1(4:6)-VN1(1:3)
          NUMLINE=NUMLINE+3
          DO J=1,NUMLINE-ELEMENTLINE(I)-NUMELE(I)
              BACKSPACE(IIN,iostat=stat)
          END DO
          READ(IIN,"(A200)",iostat=stat) LINE
          DO WHILE(TRIM(LINE)/=ENDMARK)
              IF (LINE(1:14+LEN_TRIM(SECTIONSETNAME))=="*Elset, elset="//TRIM(SECTIONSETNAME)) THEN
                  CALL SETIN(NUMSETNODE,J)
                  DO J=1,NUMSETNODE
                      MATP(SET(J))=IMAT
                  END DO
                  DEALLOCATE(SET)
              END IF
              READ(IIN,"(A200)",iostat=stat) LINE
          END DO
          DO J=1,MATERIALLINE-ENDPARTLINE(I)
              READ(IIN,*)
          END DO
          CALL MATERIALIN(DENSITY(IMAT), E(IMAT), POISSON(IMAT))
          G(IMAT)=E(IMAT)/2/(1+POISSON(IMAT))
          DO J=1,BOUNDARYLINE-NUMLINE
              BACKSPACE(IIN,iostat=stat)
          END DO
      END IF
      READ(IIN,"(A200)",iostat=stat) LINE
      NUMLINE=NUMLINE+1
  END DO
  
  DO J=1,NUMMAT
      WRITE (IOUT,"(I10,E20.12,F20.12,I10,3E20.12,4F20.12)") J, E(J), AREA(J), INDX(J), G(J), E(J)*IY(J), E(J)*IZ(J), DENSITY(J), YMAX(J), ZMAX(J)
  END DO
  
  DO J=1,ENDPARTLINE(I)-ELEMENTLINE(I)
      BACKSPACE(IIN,iostat=stat)
  END DO
  CALL ELEMENTOUT_BEAM(MATP,APPN1)
  
  
END SUBROUTINE BEAMS

SUBROUTINE BEAM_SECTION_BOX(AREA, IY, IZ, YMAX, ZMAX, A, B, T1, T2, T3, T4)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   由梁的截面参数求梁的截面面积及惯性矩                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  real(8) :: AREA, IY, IZ, YMAX, ZMAX, A, B, T1, T2, T3, T4
  real(8) :: SY, SZ, YC, ZC
  
  AREA=A*(T2+T4)+B*(T1+T3)
  SY=B*A*T2+0.5D0*(B**2)*(T1+T3)
  SZ=A*B*T1+0.5D0*(A**2)*(T2+T4)
  YC=SY/AREA
  ZC=SZ/AREA
  IF (2*YC>A) THEN
      YMAX=YC
  ELSE
      YMAX=A-YC
  END IF
  IF (2*ZC>A) THEN
      ZMAX=ZC
  ELSE
      ZMAX=A-ZC
  END IF
  IY=(B**2)*A*T2+0.3333333333333D0*(B**3)*(T1+T3)-(YC**2)*AREA
  IZ=(A**2)*B*T2+0.3333333333333D0*(A**3)*(T1+T3)-(ZC**2)*AREA
  
END SUBROUTINE BEAM_SECTION_BOX

SUBROUTINE ELEMENTOUT_BEAM(MATP,APPN1)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取并输出梁单元数据                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  integer :: MATP(NUMELE(I)), IEN(2), II, JJ, HGI, HGJ
  real(8) :: APPN1(3,NUMMAT), T(3), T_T, COR(6), LY, LZ
  integer :: J, K, L, NUMELEI
  
  NUMELEI=NUMELE(I)
  DO J=1,NUMELEI
      READ (IIN,"(I,2(',',I))") L, (IEN(K),K=1,2)
      II=N(INNUMNP(I)+IEN(1))
      JJ=N(INNUMNP(I)+IEN(2))
      HGI=0
      HGJ=0
      !IF (ID(4,II)==1) HGI=1
      !IF (ID(4,JJ)==1) HGJ=1
      T(1)=X(JJ)-X(II)
      T(2)=Y(JJ)-Y(II)
      T(3)=Z(JJ)-Z(II)
      T_T=SQRT(T(1)**2+T(2)**2+T(3)**2)
      T=T/T_T
      COR(4)=T(2)*APPN1(3,MATP(J))-T(3)*APPN1(2,MATP(J))
      COR(5)=T(3)*APPN1(1,MATP(J))-T(1)*APPN1(3,MATP(J))
      COR(6)=T(1)*APPN1(2,MATP(J))-T(2)*APPN1(1,MATP(J))
      LZ=SQRT(COR(4)**2+COR(5)**2+COR(6)**2)
      COR(4)=COR(4)/LZ
      COR(5)=COR(5)/LZ
      COR(6)=COR(6)/LZ
      COR(1)=COR(5)*T(3)-COR(6)*T(2)
      COR(2)=COR(6)*T(1)-COR(4)*T(3)
      COR(3)=COR(4)*T(2)-COR(5)*T(1)
      WRITE (IOUT,"(6I10,6F20.12))") J, II, JJ, MATP(J), HGI, HGJ, (COR(K),K=1,6)
  END DO
  
END SUBROUTINE ELEMENTOUT_BEAM

SUBROUTINE SHELL()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取并输出壳单元数据                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  integer :: MATP(NUMELE(I))
  real(8) :: E(NUMMAT), POISSON(NUMMAT), THICKNESS(NUMMAT), DENSITY(NUMMAT)    ! 
  integer :: J, K, IMAT, NUMSETNODE
  character*20 :: ENDMARK
  
  WRITE (IOUT,"(3I10)") (NPAR(J),J=1,3)
  IMAT=0
  IF (ISMOP(I)) THEN
      ENDMARK="*End Part"
  ELSE
      ENDMARK="*End Instance"
  END IF
  DO WHILE(TRIM(LINE)/=ENDMARK)
      IF (LINE(1:10)=="** Section") THEN
          IMAT=IMAT+1
          READ(IIN,"(A200)",iostat=stat) LINE
          J=INDEX(LINE,"elset=")
          K=INDEX(LINE,"material=")
          SECTIONSETNAME=LINE(J+6:K-3)
          MATERIALNAME=LINE(K+9:)
          READ(IIN,"(F)",iostat=stat) THICKNESS(IMAT)
          NUMLINE=NUMLINE+2
          DO J=1,NUMLINE-ELEMENTLINE(I)-NUMELE(I)
              BACKSPACE(IIN,iostat=stat)
          END DO
          READ(IIN,"(A200)",iostat=stat) LINE
          DO WHILE(TRIM(LINE)/=ENDMARK)
              IF (LINE(1:14+LEN_TRIM(SECTIONSETNAME))=="*Elset, elset="//TRIM(SECTIONSETNAME)) THEN
                  CALL SETIN(NUMSETNODE,J)
                  DO J=1,NUMSETNODE
                      MATP(SET(J))=IMAT
                  END DO
                  DEALLOCATE(SET)
              END IF
              READ(IIN,"(A200)",iostat=stat) LINE
          END DO
          DO J=1,MATERIALLINE-ENDPARTLINE(I)
              READ(IIN,*)
          END DO
          CALL MATERIALIN(DENSITY(IMAT), E(IMAT), POISSON(IMAT))
          DO J=1,BOUNDARYLINE-NUMLINE
              BACKSPACE(IIN,iostat=stat)
          END DO
      END IF
      READ(IIN,"(A200)",iostat=stat) LINE
      NUMLINE=NUMLINE+1
  END DO
  
  DO J=1,NUMMAT
      WRITE (IOUT,"(I10,E20.12,6F20.12)") J, E(J), POISSON(J), THICKNESS(J), DENSITY(J)
  END DO
  
  DO J=1,ENDPARTLINE(I)-ELEMENTLINE(I)
      BACKSPACE(IIN,iostat=stat)
  END DO
  CALL ELEMENTOUT(4,MATP)  
  
END SUBROUTINE SHELL

SUBROUTINE QUADS()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   读取并输出4Q单元数据                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  IMPLICIT NONE
  integer :: MATP(NUMELE(I))
  real(8) :: E(NUMMAT), POISSON(NUMMAT), DENSITY(NUMMAT)   ! 
  integer :: J, K, IMAT, NUMSETNODE
  character*20 :: ENDMARK
  
  WRITE (IOUT,"(3I10)") (NPAR(J),J=1,3)
  IMAT=0
  IF (ISMOP(I)) THEN
      ENDMARK="*End Part"
  ELSE
      ENDMARK="*End Instance"
  END IF
  DO WHILE(TRIM(LINE)/=ENDMARK)
      IF (LINE(1:10)=="** Section") THEN
          IMAT=IMAT+1
          READ(IIN,"(A200)",iostat=stat) LINE
          J=INDEX(LINE,"elset=")
          K=INDEX(LINE,"material=")
          SECTIONSETNAME=LINE(J+6:K-3)
          MATERIALNAME=LINE(K+9:)
          NUMLINE=NUMLINE+1
          DO J=1,NUMLINE-ELEMENTLINE(I)-NUMELE(I)
              BACKSPACE(IIN,iostat=stat)
          END DO
          READ(IIN,"(A200)",iostat=stat) LINE
          DO WHILE(TRIM(LINE)/=ENDMARK)
              IF (LINE(1:14+LEN_TRIM(SECTIONSETNAME))=="*Elset, elset="//TRIM(SECTIONSETNAME)) THEN
                  CALL SETIN(NUMSETNODE,J)
                  DO J=1,NUMSETNODE
                      MATP(SET(J))=IMAT
                  END DO
                  DEALLOCATE(SET)
              END IF
              READ(IIN,"(A200)",iostat=stat) LINE
          END DO
          DO J=1,MATERIALLINE-ENDPARTLINE(I)
              READ(IIN,*)
          END DO
          CALL MATERIALIN(DENSITY(IMAT), E(IMAT), POISSON(IMAT))
          DO J=1,BOUNDARYLINE-NUMLINE
              BACKSPACE(IIN,iostat=stat)
          END DO
      END IF
      READ(IIN,"(A200)",iostat=stat) LINE
      NUMLINE=NUMLINE+1
  END DO
  
  DO J=1,NUMMAT
      WRITE (IOUT,"(I10,E20.12,5F20.12)") J, E(J), POISSON(J), DENSITY(J)
  END DO
  
  DO J=1,ENDPARTLINE(I)-ELEMENTLINE(I)
      BACKSPACE(IIN,iostat=stat)
  END DO
  CALL ELEMENTOUT(4,MATP)
  
END SUBROUTINE QUADS


END PROGRAM PREPROCESS