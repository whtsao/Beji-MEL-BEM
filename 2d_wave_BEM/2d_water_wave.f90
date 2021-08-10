!===============================================================================
!	SOLVE 2D WATER WAVE PROBLEM WITH ONE SIDE WAVEMAKER AND ONESIDE FREE BY BEM
!	WEN-HAUI TSAO 2021 LSU
!===============================================================================
      PROGRAM WATER_WAVE
      IMPLICIT NONE
	  INTEGER I,J,NT,NGA,NTIM,NFIELD,NNODE,NELM,ICON,NITER
	  INTEGER NELEM(4),NS(4),BUNTYP(4)
	  INTEGER,ALLOCATABLE::LN(:,:)
	  REAL DEP,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,DELTTIME,SOR
	  REAL TIME,DIS,VEL,ACC,P_ATM,FOR,DDIS,DTEMP,WC,E1,E2,ETOL
	  REAL COOR(4,2),SIDE_L(4)
	  REAL,ALLOCATABLE::NODE(:,:),NORM(:,:),JCB(:),LENG(:),PHI(:),PPHI(:),PHIT(:),PPHIT(:)
	  REAL,ALLOCATABLE::KER1(:,:),KER2(:,:),DP(:,:),DPDS(:),PR(:),DPDT(:),ACCMO(:)
	  REAL,ALLOCATABLE::PHIT_TEMP(:)
	  REAL,ALLOCATABLE::WT(:),RT(:),SHA1(:),SHA2(:),SH(:,:)

      OPEN(UNIT=1,FILE='1.IPT',STATUS='OLD')
      OPEN(UNIT=5,FILE='IO.DAT')
      OPEN(UNIT=6,FILE='S.DAT')
      OPEN(UNIT=7,FILE='W.DAT')
      OPEN(UNIT=8,FILE='P.DAT')
      OPEN(UNIT=9,FILE='F.DAT')
      OPEN(UNIT=10,FILE='E.DAT')
      OPEN(UNIT=11,FILE='DOMAIN.DAT')
	  OPEN(UNIT=21,FILE='ERR.DAT')
      OPEN(UNIT=99,FILE='TEST.TXT')

	CALL INPUT(COOR,NFIELD,NNODE,NELM,NELEM,NS,BUNTYP,DEP,NGA,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,NTIM,DELTTIME,SOR,ICON,NITER,ETOL,WC)
	ALLOCATE(LN(NELM,2),NODE(NNODE,2),NORM(NELM,2),JCB(NELM),LENG(NELM),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE))
	ALLOCATE(KER1(NNODE,NNODE),KER2(NNODE,NNODE),DP(NNODE,2),DPDS(NNODE),DPDT(NNODE),PR(NNODE),ACCMO(NNODE))
	ALLOCATE(PHIT_TEMP(NNODE))
	ALLOCATE(WT(NGA),RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA))
	LN=0
	NODE=0.0
	NORM=0.0
	JCB=0.0
	LENG=0.0
	PHI=0.0
	PPHI=0.0
	PHIT=0.0
	PPHIT=0.0
	KER1=0.0
	KER2=0.0
	DP=0.0
	DPDS=0.0
	DPDT=0.0
	PR=0.0
	ACCMO=0.0
	PHIT_TEMP=0.0
	WT=0.0
	RT=0.0
	SHA1=0.0
	SHA2=0.0
	SH=0.0
    CALL GAUSS(WT,RT,NGA)
    CALL SHAP(SHA1,SHA2,SH,NGA,RT)
	CALL LENGTH(COOR,SIDE_L)
	CALL MESH(NNODE,NELM,NELEM,LN,COOR,SIDE_L,NODE) 
    WRITE(*,*) 'PASS MESH'

DO NT=1,NTIM
	WRITE(*,*) NT,'TH'
     TIME=(NT-1)*DELTTIME

!---WAVEMAKER DIS, VEL, ACC, PHASE LAG (IN RADIUS)
	 DIS=AMP*SIN(OMEGA*TIME+PSI) !-AMP_B*COS(OMEGA_B*TIME+PSI) !
     VEL=AMP*OMEGA*COS(OMEGA*TIME+PSI) !AMP_B*OMEGA_B*SIN(OMEGA_B*TIME+PSI) !
     ACC=-AMP*OMEGA**2*SIN(OMEGA*TIME+PSI) !AMP_B*OMEGA_B**2*COS(OMEGA_B*TIME+PSI) !
	 DDIS=DIS-DTEMP
	 DTEMP=DIS
	 WRITE(5,"(7(E15.8,1X))") TIME,DIS,VEL,ACC

!---REMESH LOCATION AND POTENTIAL OF A FREE-SURFACE NODE AT CURRENT TIME, I.E. (NT-1)TH TIME-STEP (BUT THE FOLLOWING CALCULATION IS FOR THE NEXT TIME-STEP)
      CALL REMESH(NNODE,NELEM,NODE,NS,TIME,DEP,AMP)
      WRITE(*,*) 'PASS REMESH'

!---BUILD KERNEL OF ZONE1, ZONE2, ZONE3 AND ASSEMBLE THEM INTO BIG KERNEL
	   CALL KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT)
       WRITE(*,*) 'PASS KERNEL'

!---CALCULATE PRESSURE, FORCE, ENERGY AND POWER ON THE BOUNDARY AND IN THE DOMAIN
!---AT CURRENT TIME, I.E. (NT-1)TH TIME-STEP (BUT THE FOLLOWING CALCULATION IS FOR THE NEXT TIME-STEP)
	  CALL PRESSURE(ICON,TIME,THO,GRAV,DEP,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)

!      CALL FORCE(TIME,WIDTH,NNODE,NS,NODE,PR,FOR)
!	  CALL POWER(NNODE1,NELM1,LN1,NS1,NODE1,PHI1,PPHI1,PHIT1,PPHIT1,PR1,&
!				&NNODE2,NELM2,LN2,NS2,NODE2,PHI2,PPHI2,PHIT2,PPHIT2,PR2,&
!				&NNODE3,NELM3,LN3,NS3,NODE3,PHI3,PPHI3,PHIT3,PPHIT3,PR3,&
!				&TIME,GRAV,THO,WIDTH,DEP,BROOT)
	  CALL DOMAIN(NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)
!       WRITE(*,*) 'PASS PFE'

!**************IMPLICIT SOLVER FOR PHIT AND PPHIT ON THE ABSORPTION SIDE**************
DO J=1,NITER

	  PHIT_TEMP=PHIT
!---APPLY BC FOR SOLVING FREE-SURFACE PPHI
       CALL BOUND(NNODE,NS,PHI,PPHI,PHIT,VEL,WC)
	   CALL SOLVE_BACK(PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)

!---FIRST-ORDER TAYLOR SERIES EXPANSION (ALSO GET TANGENTIAL VELOCITY ON ZONE BOUNDARY)
	  CALL TAYES1(NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,DEP,GRAV,MU,VEL,WC)

!---APPLY BC FOR SOLVING FREE-SURFACE PPHIT
	  CALL ACCBC(NNODE,NELM,NELEM,LN,NORM,JCB,PPHI,DP,ACCMO)
      CALL BOUNDT(NNODE,NELM,NS,NELEM,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
	  CALL SOLVE_BACK(PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)

CALL PHIT_CONVERGE(NELEM(2)+1,PHIT_TEMP(NS(1)+1:NS(2)),PHIT(NS(1)+1:NS(2)),E1,E2)
	IF(E1<=ETOL.AND.E2<ETOL)THEN
		WRITE(*,*) 'PASS CONVERGED',J
		WRITE(21,*) TIME,J,E1,E2
		GOTO 205
	ELSE IF(J>=NITER)THEN
		WRITE(*,*) 'CONVERGE FAIL'
		STOP
	END IF

END DO
!**************************************************************************************

205 CONTINUE

!---SECOND-ORDER TAYLOR SERIES EXPANSION
      CALL TAYES2(PHI,PPHI,PHIT,PPHIT,DPDS,DPDT,DP,NNODE,NELEM,NODE,NELM,NORM,JCB,LN,DELTTIME,GRAV,ACC)
      WRITE(*,*) 'PASS TAYES'

END DO

STOP
END
!**********************************************************************
      SUBROUTINE HEADLINE(ID,IREAD)
!**********************************************************************
      CHARACTER*2 ID
      ID =  '*'
      DO WHILE (ID .EQ. '*')
      READ(IREAD,'(A1)') ID
      END DO
      RETURN
      END
!**********************************************************************
      SUBROUTINE INPUT(COOR,NFIELD,NNODE,NELM,NELEM,NS,BUNTYP,DEP,NGA,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,NTIM,DELTTIME,SOR,ICON,NITER,ETOL,WC)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER I,J,NFIELD,NNODE,NELM,NGA,NTIM,ICON,NITER
	  INTEGER NELEM(4),NS(4),BUNTYP(4)
	  REAL    COOR(4,2),DEP,ENDTIME,DELTTIME,SOR,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,ETOL,WC,L0,T
      CHARACTER*2 ID
         ID = '*'

!=========ZONE 1==========
!---NODAL COORDINATES
         CALL HEADLINE(ID,1)
         READ(1,*)  ((COOR(I,J),J=1,2),I=1,4)
!---ELEMENT MESH NUMBER
         CALL HEADLINE(ID,1)
        READ(1,*) (NELEM(I),I=1,4)
        NNODE = 0
        NELM  = 0
        DO I=1,4
         NELM = NELM+NELEM(I)
         NNODE = NNODE+(NELEM(I)+1)
        END DO
      NS(1)=NELEM(1)+1
      NS(2)=NS(1)+NELEM(2)+1
      NS(3)=NS(2)+NELEM(3)+1
      NS(4)=NS(3)+NELEM(4)+1
	  NFIELD=(NELEM(1)-1)*(NELEM(2)-1)

!---BOUNDARY TYPE
         CALL HEADLINE(ID,1)
        READ(1,*) (BUNTYP(I),I=1,4)
!---READ THE GAUSSIAN INTEGRATION POINT NO.
         CALL HEADLINE(ID,1)
         READ(1,*)  NGA
!---READ TIME,TIME STEP,GRAV ACC.,FREQUENCY,AMPLITUDE
         CALL HEADLINE(ID,1)
         READ(1,*)  GRAV,MU,WIDTH,THO
!---READ TIME,TIME STEP,GRAV ACC.,FREQUENCY,AMPLITUDE
         CALL HEADLINE(ID,1)
         READ(1,*)  AMP,OMEGA,PSI
!---READ TIME,TIME STEP,GRAV ACC.,FREQUENCY,AMPLITUDE
         CALL HEADLINE(ID,1)
         READ(1,*)  ENDTIME,DELTTIME,SOR,ICON,NITER,ETOL
		NTIM=ENDTIME/DELTTIME+1
		DEP=COOR(4,2)

!---CALCULATE WAVE SPEED
T=2.0*ACOS(-1.0)/OMEGA
L0=GRAV*T**2/(2.0*ACOS(-1.0))
WC=L0/T

      RETURN
      END
!**********************************************************************
      SUBROUTINE LENGTH(COOR1,SIDE_L1)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      REAL  COOR1(4,2),SIDE_L1(4)
      DO I=1,3
        SIDE_L1(I)=SQRT((COOR1(I+1,1)-COOR1(I,1))**2+(COOR1(I+1,2)-COOR1(I,2))**2)
      END DO
        SIDE_L1(4)=SQRT((COOR1(4,1)-COOR1(1,1))**2+(COOR1(4,2)-COOR1(1,2))**2)

      RETURN
      END
!**********************************************************************
      SUBROUTINE MESH(NNODE1,NELM1,NELEM1,LN1,COOR1,LENG1,NODE1)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      INTEGER NEL(4),NNODE1,NELM1,NELEM1(4),LN1(NELM1,2)
      REAL    DELT,LENG1(4),COOR1(4,2),NODE1(NNODE1,2)

!** THE 1ST SURFACE (X,NORMAL : Y)********************************
      NEL(1) = NELEM1(1)+1
      DELT=LENG1(3)/NELEM1(1)
      DO I=1,NEL(1)
         NODE1(I,1)=COOR1(4,1)+(I-1)*DELT
         NODE1(I,2)=COOR1(4,2)
      END DO
      K=NEL(1)
!** THE 2ND SURFACE (Y,NORMAL : -X)********************************
      NEL(2) = NELEM1(2)+1
      DELT=LENG1(2)/NELEM1(2)
      DO I=1,NELEM1(2)+1
        NODE1(I+K,1)=COOR1(3,1)
        NODE1(I+K,2)=COOR1(3,2)-(I-1)*DELT
      END DO
      K=K+NEL(2)
!** THE 3RD SURFACE (X,NORMAL : -Y)********************************
      NEL(3) = NELEM1(3)+1
      DELT=LENG1(1)/NELEM1(3)
      DO I=1,NELEM1(3)+1
       NODE1(I+K,1)=COOR1(2,1)-(I-1)*DELT
       NODE1(I+K,2)=COOR1(2,2)
      END DO
      K=K+NEL(3)

!** THE 4TH SURFACE (Y,NORMAL : X)********************************
      NEL(4) = NELEM1(4)+1
      DELT=LENG1(4)/NELEM1(4)
      DO I=1,NELEM1(4)+1
        NODE1(I+K,1)=COOR1(1,1)
        NODE1(I+K,2)=COOR1(1,2)+(I-1)*DELT
      END DO

!----TO GIVE THE LOCAL ELEMENT NODE NUMBER
      L=1
	  N=1
      DO I=1,4
       DO J=1,NELEM1(I)
        LN1(N,1)=L
        LN1(N,2)=L+1
        L=L+1
		N=N+1
       END DO
       L=L+1
      END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE SHAP(SHA1,SHA2,SH,NGA,RT)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL (A-H,O-Z)
      INTEGER NGA
      REAL RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
      DO M=1,NGA
        SHA1(M)=0.5*(1-RT(M))
        SHA2(M)=0.5*(1+RT(M))

        SH(1,M)=SHA1(M)
        SH(2,M)=SHA2(M)
      END DO
	RETURN
	END 
!**********************************************************************
      SUBROUTINE REMESH(NNODE,NELEM,NODE,NS,TIME,DEP,AMP)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL(A-H,O-Z)
      INTEGER NNODE,NELEM(4),NS(4)
      REAL TIME,DEP,AMP,LENG,NODE(NNODE,2)

!------DUPLICATE POINT ON FS-----
    NODE(NNODE,:)=NODE(1,:)
	NODE(NS(1)+1,:)=NODE(NS(1),:)

!**** THE 1ST SURFACE (X,NORMAL : Y)********************************

!**** THE 2ND SURFACE (Y,NORMAL : -X)********************************
      LENG=NODE(NS(1)+1,2)-NODE(NS(2),2)
      DELT=LENG/NELEM(2)
	  K=1
      DO I=NS(1)+1,NS(2)
        NODE(I,1)=NODE(NS(1),1)
        NODE(I,2)=NODE(NS(1),2)-DELT*(K-1)
		K=K+1
      END DO

!**** THE 3RD SURFACE (X,NORMAL : -Y)********************************
      LENG=NODE(NS(1),1)-NODE(1,1)
      DELT=LENG/NELEM(3)
	  K=1
      DO I=NS(2)+1,NS(3)
       NODE(I,1)=NODE(NS(1),1)-DELT*(K-1)
       NODE(I,2)=0.0
	   K=K+1
      END DO

!**** THE 4TH SURFACE (Y,NORMAL : X)********************************
      LENG=NODE(1,2)-NODE(NS(3),2)
      DELT=LENG/NELEM(4)
	  K=1
      DO I=NS(3)+1,NNODE
        NODE(I,1)=NODE(1,1)
        NODE(I,2)=DELT*(K-1)
		K=K+1
      END DO

!==========OUTPUT FREE-SURFACE AND WALL AND BAFFLE BOUNDARY==========
WAV_L=NODE(1,2)-DEP
WAV_R=NODE(NELEM(1)+1,2)-DEP
WRITE(6,"(500(E15.8,1X))") NODE(:,1) 
WRITE(6,"(500(E15.8,1X))") NODE(:,2)
WRITE(7,"(5(E15.8,1X))") TIME,WAV_L,WAV_R,WAV_L/AMP,WAV_R/AMP

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUND(NNODE,NS,PHI,PPHI,PHIT,VEL,WC)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,NS(4)
	   REAL R,VEL,WC
       REAL PHI(NNODE),PPHI(NNODE),PHIT(NNODE)

      DO I=1,NS(1)
         PHI(I)=PHI(I)
         PPHI(I)=0.0
      END DO

      DO I=NS(1)+1,NS(2)
       PPHI(I)=-PHIT(I)/WC !VEL !

       PHI(I)=0.0
      END DO

      DO I=NS(2)+1,NS(3)
	   PPHI(I)=0.0
       PHI(I)=0.0
      END DO

      DO I=NS(3)+1,NNODE
       PPHI(I)=-VEL
       PHI(I)=0.0
      END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      INTEGER  NNODE,NELM,NGA,LN(NELM,2)
      REAL     KER1(NNODE,NNODE),KER2(NNODE,NNODE)
      REAL     NX,NY,H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
      REAL     LENG(NELM),NORM(NELM,2),JCB(NELM),NODE(NNODE,2)
      REAL     WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)

        KER1=0.0
        KER2=0.0
!**** CALCULATE THE JACOBIAN
      DO J=1,NELM
      LENG(J)=((NODE(LN(J,1),1)-NODE(LN(J,2),1))**2+(NODE(LN(J,1),2)-NODE(LN(J,2),2))**2)**0.5
		DO L=1,2
          PXI1(L)=(-0.5)*NODE(LN(J,1),L)+ 0.5*NODE(LN(J,2),L)
        END DO
        NX=PXI1(2)
        NY=PXI1(1)
        JCB(J)=(NX**2+NY**2)**0.5
        NORM(J,1)=-NX/JCB(J)
        NORM(J,2)=NY/JCB(J)
      END DO

!***THE SURFACE KERNELS
      DO I = 1,NNODE
       DO J=1,NELM
       DO M=1,NGA
          XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
          YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
       END DO
        DO K=1,2
         G(K)=0.0
         H(K)=0.0
         RD=((NODE(I,1)-NODE(LN(J,K),1))**2+(NODE(I,2)-NODE(LN(J,K),2))**2)**0.5
        IF (RD .LE. 0.000001) THEN
        H(K)=0.0
	    G(K)=LENG(J)/2*(1.5-LOG(LENG(J)))
		ELSE !---NON SINGULER TERM
         G(K)=0.0
         DO M=1,NGA
            H(K)=H(K)+(-1)/((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)*((XFUNC(M)-NODE(I,1))*NORM(J,1)+(YFUNC(M)-NODE(I,2))*NORM(J,2))*JCB(J)*SH(K,M)*WT(M)
            G(K)=G(K)+LOG(1./((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)**0.5)*JCB(J)*SH(K,M)*WT(M)
         END DO
      END IF
         KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
      END DO
      END DO
      END DO
!***SINGULAR OF KER1
       DO I=1,NNODE
          SIGMA=0.0
          DO J=1,NNODE
             SIGMA=KER1(I,J)+SIGMA
          END DO
          KER1(I,I)=-SIGMA
       END DO

      RETURN
      END
!**********************************************************************
       SUBROUTINE SOLVE_BACK(PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!      USING GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION
!======================================================================
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,NELEM(4),BUNTYP(4)
       REAL    PHI(NNODE),PPHI(NNODE)
       REAL    KER1(NNODE,NNODE),KER2(NNODE,NNODE),GG(NNODE,NNODE)
       REAL    H1(500,500),Q1(500),TEMP(500)
       REAL*8  SUM,A,SIG,G1(500,500),P1(500)

!********** TO MOVE THE KER1 AND KER2**********************************
!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,4
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO
!-------------------
         DO I=1,NNODE
           N=0
           DO L=1,4
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       DO I=1,NNODE
          TEMP(I)=0.0
          DO J=1,NNODE
          TEMP(I)=TEMP(I)+H1(I,J)*Q1(J)
          END DO
       END DO
!*************GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION*********
      DO I=1,NNODE
        G1(I,NNODE+1)=TEMP(I)
      END DO
      DO I=1,NNODE
         SUM=0.0
         DO K=I,NNODE
            IF (G1(K,I) .NE. 0) THEN
               IF (K .NE. I) THEN
               IF (G1(I,I) .EQ. 0) THEN
                 WRITE(*,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 STOP
               END IF
               A=G1(K,I)/G1(I,I)
               DO J=I,NNODE+1
                  G1(K,J)=G1(K,J)-A*G1(I,J)
               END DO
               END IF
            END IF
            SUM=SUM+G1(K,I)
          END DO
          IF (SUM .EQ. 0.) THEN
          WRITE(*,*) 'NO UNIQUE SOLUTION EXISTS  STOP AT GAUSSELI'
          STOP
          END IF
      END DO

      P1(NNODE)=G1(NNODE,NNODE+1)/G1(NNODE,NNODE)
      DO I=NNODE-1,1,-1
         SIG=0.0
         DO J=I+1,NNODE
            SIG=G1(I,J)*P1(J)+SIG
          END DO
         P1(I)=(G1(I,NNODE+1)-SIG)/G1(I,I)
      END DO
!=================================
       K=1
       N=0
       DO I=1,4
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
      END
!**********************************************************************
       SUBROUTINE SOLVE(SPHI,SPPHI,PHI,PPHI,KER1,KER2,BUNTYP,NNODE,SOR,NELEM)
!**********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,BUNTYP(4),NELEM(4)
       REAL    SPHI(NNODE),SPPHI(NNODE),PHI(NNODE),PPHI(NNODE)
       REAL    KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL    SOR
       REAL    TEMP(1000),H1(1000,1000),Q1(1000),P1(1000)
       REAL    G1(1000,1000),ERR(1000),G2(1000,1000),X1(1000)
       REAL    XP(1000)
       REAL*8  EMAX 

!********** TO MOVE THE KER1 AND KER2**********************************
!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,4
          DO J=K+N,NELEM(I)+N+1
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=SPHI(J)
               !P1(J)=(PPHIB(J)-PPHIA(J))+PPHIB(J)
             ELSE
               Q1(J)=SPPHI(J)
               !P1(J)=(PHIB(J)-PHIA(J))+PHIB(J)   
            END IF
          END DO
        N=N+NELEM(I)+1
       END DO
!-------------------
       DO I=1,NNODE
          N=0
          DO L=1,4
            DO J=K+N,NELEM(L)+N+1
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
            END DO
            N=N+NELEM(L)+1
          END DO
       END DO

       DO I=1,NNODE
          TEMP(I)=0.0
          DO J=1,NNODE
          TEMP(I)=TEMP(I)+H1(I,J)*Q1(J)
          END DO
       END DO

      DO I=1,NNODE
       XP(I)=P1(I)
       X1(I)=0.0
       DO J=1,NNODE
        G2(I,J)=-G1(I,J)/G1(I,I)
       END DO
        G2(I,NNODE+1)=TEMP(I)/G1(I,I)
      END DO

!*************GSI SOR INTERATIVE

      DO  L=1,300
       DO I=1,NNODE

        DO J=1,I-1
         X1(I)=X1(I)+XP(J)*G2(I,J)*SOR
        END DO

        DO J=I+1,NNODE
         X1(I)=X1(I)+XP(J)*G2(I,J)*SOR
        END DO

         X1(I)=X1(I)+G2(I,NNODE+1)*SOR+(1-SOR)*XP(I)
         XP(I)=X1(I)
       END DO
        PMAX=0.0
        EMAX=0.0
       DO I=1,NNODE
        ERR(I)=ABS(P1(I)-X1(I))
        IF (ERR(I) .GT. EMAX) THEN
         EMAX=ERR(I)
        END IF
        IF (ABS(P1(I)) .GT. PMAX) THEN
         PMAX=ABS(P1(I))
        END IF
        P1(I)=X1(I)
        X1(I)=0.0
       END DO

      IF(PMAX .EQ. 0) THEN
        IF (EMAX .LE. 0.0000001) GOTO 100
      ELSE
        ETOL=EMAX/PMAX
        IF (ETOL .LE. 0.0001) GOTO 100
      END IF
      END DO
  100 WRITE(*,*) 'L',L

       K=1
       N=0
       DO I=1,4
          DO J=K+N,NELEM(I)+N+1
            IF (BUNTYP(I) .EQ. 1) THEN
               SPPHI(J)=P1(J)
              ELSE
               SPHI(J)=P1(J)
            END IF
          END DO
          N=NELEM(I)+N+1
       END DO

      DO I=1,NNODE
       PHI(I)=SPHI(I)
       PPHI(I)=SPPHI(I)
      END DO

      RETURN
      END
!**********************************************************************
       SUBROUTINE SOLVE_OLD(SPHI,SPPHI,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP,NS,SOR)
!**********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,NS(4),NELEM(4),BUNTYP(4)
       REAL    SOR
	   REAL    SPHI(NNODE),SPPHI(NNODE),PHI(NNODE),PPHI(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL    G1(NNODE,NNODE+1),P1(NNODE),H1(NNODE,NNODE),Q1(NNODE),TEMP(NNODE)
       REAL    G2(1000,1000),X1(1000),XP(1000),ERR(1000)
       REAL*8  EMAX

	G1=0.0
	P1=0.0
	H1=0.0
	Q1=0.0
!********** TO MOVE KER1 AND KER2 INTO G1*P1(UNKNOWN)=H1*Q1(KNOWN) ****************

       K=1
!==========PUT PHI, PPHI, AND ZERO IN THE KNOWN VECTOR==========
       N=0
	   L=1
       DO I=1,4
          DO J=K+N,NELEM(I)+1+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(L)=SPHI(J)
			   L=L+1
             ELSE IF (BUNTYP(I) .EQ. 0) THEN
               Q1(L)=SPPHI(J)
			   L=L+1
            END IF 
		  END DO
        N=N+NELEM(I)+1
       END DO

!==========MOVE KER1 AND KER2 WHEN BC TYPE = 0 OR 1 ==========
	!---FOR ZONE 1---
		IROW=0
         DO I=1,NNODE
		  M=1
          N=0
          DO L=1,4
            DO J=K+N,NELEM(L)+1+N
             IF (BUNTYP(L) .EQ. 1) THEN 
               H1(IROW+I,M)=KER1(I,J)	! PHI IS GIVEN PUT KER1 IN H1
               G1(IROW+I,M)=KER2(I,J)	! PPHI IS UNKNOWN PUT KER2 IN G1
			   M=M+1
              ELSE IF (BUNTYP(L) .EQ. 0) THEN
               H1(IROW+I,M)=-KER2(I,J)	! PPHI IS GIVEN MOVE KER2 IN H1
               G1(IROW+I,M)=-KER1(I,J)	! PHI IS UNKNOWN MOVE KER1 IN G1
			   M=M+1
             END IF
            END DO
           N=NELEM(L)+1+N
          END DO
         END DO

!==========GSI SOR INTERATIVE==========
       DO I=1,NNODE
          TEMP(I)=0.0
          DO J=1,NNODE
          TEMP(I)=TEMP(I)+H1(I,J)*Q1(J)
          END DO
       END DO

      DO I=1,NNODE
       XP(I)=P1(I)
       X1(I)=0.0
       DO J=1,NNODE
        G2(I,J)=-G1(I,J)/G1(I,I)
       END DO
        G2(I,NNODE+1)=TEMP(I)/G1(I,I)
      END DO

      DO  L=1,300
       DO I=1,NNODE
        DO J=1,I-1
         X1(I)=X1(I)+XP(J)*G2(I,J)*SOR
        END DO
        DO J=I+1,NNODE
         X1(I)=X1(I)+XP(J)*G2(I,J)*SOR
        END DO
         X1(I)=X1(I)+G2(I,NNODE+1)*SOR+(1-SOR)*XP(I)
         XP(I)=X1(I)
       END DO
        PMAX=0.0
        EMAX=0.0
       DO I=1,NNODE
        ERR(I)=ABS(P1(I)-X1(I))
        IF (ERR(I) .GT. EMAX) THEN
         EMAX=ERR(I)
        END IF
		IF (ABS(P1(I)) .GT. PMAX) THEN
         PMAX=ABS(P1(I))
        END IF
		P1(I)=X1(I)
        X1(I)=0.0
       END DO
       IF(PMAX .EQ. 0) THEN
        IF (EMAX .LE. 0.0000001) GOTO 100
       ELSE
        ETOL=EMAX/PMAX
        IF (ETOL .LE. 0.0001) GOTO 100
       END IF
      END DO
  100 WRITE(*,*) 'ITER',L !

!==========MOVE BACK SOLUTION==========
	   K=1
       N=0
	   L=1
       DO I=1,4
          DO J=K+N,NELEM(I)+1+N
            IF (BUNTYP(I) .EQ. 1) THEN
               SPPHI(J)=P1(L)
			   L=L+1
             ELSE IF (BUNTYP(I) .EQ. 0) THEN
               SPHI(J)=P1(L)
			   L=L+1
            END IF 
		  END DO
        N=N+NELEM(I)+1
       END DO

!---OUTPUT NEW PHI AND PPHI
      DO I=1,NNODE
       PHI(I)=SPHI(I)
       PPHI(I)=SPPHI(I)
      END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE TAYES1(NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,DEP,GRAV,MU,VEL,WC)
!********************************************************************
    IMPLICIT INTEGER NONE
    INTEGER I,J,K,NNODE,NELM,NELEM(4),NS(4),LN(NELM,2)
    REAL WC,DEP,GRAV,MU,VEL
	REAL NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2),DPDS(NNODE),PHIT(NNODE),DPDT(NNODE)

	DPDS=0.0
	K=0
!*********************ON FREE SURFACE*********************
    DO I=K+1,K+NELEM(1)
    DO J=1,2
		IF(LN(I,J).EQ.1) THEN
            DP(LN(I,J),1)=-PPHI(NNODE)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(1)) THEN
		    DP(LN(I,J),1)=PPHI(NS(1)+1)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5*PHI(LN(I,1))+0.5*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

    DO I=1,NS(1)
	  DPDT(I)=0.5*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*PHI(I) !-GRAV*NODE(I,2)-MU*PHI(I) !
      PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
    END DO

K=K+NELEM(1)
!*********************ON RIGHT WALL*********************
    DO I=K+1,K+NELEM(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(1)+1) THEN
			DPDS(LN(I,J))=-DP(NS(1),2)
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(2)) THEN
			DPDS(LN(I,J))=0.0
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5*PHI(LN(I,1))+0.5*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=PPHI(LN(I,J))
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++HOW TO GET PHIT ON THE RIGHT WALL FOR RADIATION CONDITION
!++++++++++++++++++++++++++++++++++++++++++THE BERNOULLI EQUATION IS USELESS BECAUSE I DON'T HAVE ATM PRESSURE
!	DO I=NS(1)+1,NS(2)
!	  DPDT(I)=0.5*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*PHI(I) !-GRAV*NODE(I,2)-MU*PHI(I) !
!      PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
!    END DO

K=K+NELEM(2)
!*********************ON BOTTOM*********************
    DO I=K+1,K+NELEM(3)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)+1) THEN
            DPDS(LN(I,J))=-PPHI(NS(2))
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(3)) THEN
            DPDS(LN(I,J))=PPHI(NS(3)+1)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5*PHI(LN(I,1))+0.5*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

K=K+NELEM(3)
!*********************LEFT WALL*********************
    DO I=K+1,K+NELEM(4)
    DO J=1,2
		IF(LN(I,J).EQ.NS(3)+1) THEN
			DPDS(LN(I,J))=0.0
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NNODE) THEN
			DPDS(LN(I,J))=DP(1,2)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5*PHI(LN(I,1))+0.5*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE ACCBC(NNODE,NELM,NELEM,LN,NORM,JCB,PPHI,DP,ACCMO)
!********************************************************************
    IMPLICIT INTEGER NONE
    INTEGER I,J,K,NNODE,NELM,NELEM(4),LN(NELM,2)
    REAL DPDNX,DPDNY,DPNDS
	REAL NORM(NELM,2),JCB(NELM),PPHI(NNODE),DP(NNODE,2),ACCMO(NNODE)

	K=0
!    DO I=K+1,K+NELEM(1)
!    DO J=1,2
!	DPNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
!	DPNDX=DPNDS*NORM(I,2)	! PHINN=PHISS=0 FOR LINEAR ELEMENT
!	DPNDY=-DPNDS*NORM(I,1)
!	ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
!    END DO
!    END DO

	K=K+NELEM(1)
    DO I=K+1,K+NELEM(2)
    DO J=1,2
		DPNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
		DPNDX=DPNDS*NORM(I,2)	! PHINN=PHISS=0 FOR LINEAR ELEMENT
		DPNDY=-DPNDS*NORM(I,1)
		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
    END DO

	K=K+NELEM(2)
    DO I=K+1,K+NELEM(3)
    DO J=1,2
		DPNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
		DPNDX=DPNDS*NORM(I,2)	! PHINN=PHISS=0 FOR LINEAR ELEMENT
		DPNDY=-DPNDS*NORM(I,1)
		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
    END DO

	K=K+NELEM(3)
	DO I=K+1,K+NELEM(4)
    DO J=1,2
		DPNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
		DPNDX=DPNDS*NORM(I,2)	! PHINN=PHISS=0 FOR LINEAR ELEMENT
		DPNDY=-DPNDS*NORM(I,1)
		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
    END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUNDT(NNODE,NELM,NS,NELEM,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,NELM,NS(4),NELEM(4),LN(NELM,2)
	   REAL R,ACC,WC,DPDSS
       REAL JCB(NELM),PHIT(NNODE),PPHIT(NNODE),DPDS(NNODE),ACCMO(NNODE)

      DO I=1,NS(1)
         PHIT(I)=PHIT(I)
         PPHIT(I)=0.0
      END DO

	PPHIT=0.0
    DO I=NELEM(1)+1,NELEM(1)+NELEM(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)) THEN
		  DPDSS=0.0 ! ASSUME A FLAT GROUND SO VN=0
		ELSE
		  DPDSS=0.5*(-DPDS(LN(I,1))+DPDS(LN(I,2)))/JCB(I)
		END IF
      PPHIT(LN(I,J))=WC*DPDSS
      PHIT(LN(I,J))=0.0
    END DO
    END DO

      DO I=NS(2)+1,NS(3)
	   PPHIT(I)=0.0-ACCMO(I)
       PHIT(I)=0.0
      END DO

      DO I=NS(3)+1,NNODE
       PPHIT(I)=-ACC-ACCMO(I)
       PHIT(I)=0.0
      END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE TAYES2(PHI,PPHI,PHIT,PPHIT,DPDS,DPDT,DP,NNODE,NELEM,NODE,NELM,NORM,JCB,LN,DELTTIME,GRAV,ACC)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      INTEGER  NNODE,NELM,NELEM(4),LN(NELM,2)
      REAL DELTTIME,GRAV,ACC,D2PDT,DPTDS,DPDNDS
      REAL NORM(NELM,2),JCB(NELM)
      REAL NODE(NNODE,2),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),DP(NNODE,2),DPDS(NNODE),DPDT(NNODE)
      REAL NEW_NODE(NELEM(1)+1,2),NEW_PHI(NELEM(1)+1),D2P(2)

	DO I=1,NELEM(1)
	  DO J=1,2
		IF(LN(I,J) .EQ. 1) THEN
				D2P(1)=-PPHIT(NNODE)
				DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
				DPTDS=(D2P(1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
				D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE IF (LN(I,J) .EQ. NELEM(1)+1) THEN
				D2P(1)=PPHIT(NELEM(1)+2)
				DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
				DPTDS=(D2P(1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
				D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE
			DPTDS=(-0.5*PHIT(LN(I,1))+0.5*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
			D2P(1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		END IF

		D2PDT=DP(LN(I,J),1)*D2P(1)+DP(LN(I,J),2)*D2P(2)-GRAV*DP(LN(I,J),2)
		NEW_PHI(LN(I,J))=PHI(LN(I,J))+DELTTIME*DPDT(LN(I,J))+D2PDT*DELTTIME**2/2
		NEW_NODE(LN(I,J),1)=NODE(LN(I,J),1)+DP(LN(I,J),1)*DELTTIME+0.5*D2P(1)*DELTTIME**2
		NEW_NODE(LN(I,J),2)=NODE(LN(I,J),2)+DP(LN(I,J),2)*DELTTIME+0.5*D2P(2)*DELTTIME**2
	  END DO
	END DO

	DO I=1,NELEM(1)+1
	NODE(I,:)=NEW_NODE(I,:)
	PHI(I)=NEW_PHI(I)
	END DO

      RETURN
      END
!********************************************************************
SUBROUTINE PRESSURE(ICON,TIME,THO,GRAV,DEP,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,ICON,NNODE,NS(4)
      REAL TIME,DEP,THO,GRAV,P1,P2,P3,P_ATM
      REAL NODE(NNODE,2),PHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL CP1(NS(1)),CP2(NS(1)),CP3(NS(1)),CP(NS(1))

!----ATOM PRESSURE (BERNOULLI CONSTANT) ON THE FREE SURFACE
	DO I=ICON,ICON ! NS(1) !
	CP1(I)=THO*PHIT(I)
	CP2(I)=THO*0.5*(DP(I,1)**2+DP(I,2)**2)
	CP3(I)=THO*GRAV*(NODE(I,2)-DEP)
	CP(I)=CP1(I)+CP2(I)+CP3(I)
	ENDDO
	P_ATM=CP(ICON)

!----PRESSURE ON ZONE 1 BOUNDARY
	DO I=NS(1)+2,NNODE-1
	P1=THO*PHIT(I)
	P2=THO*0.5*(DP(I,1)**2+DP(I,2)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP)
	PR(I)=P_ATM-(P1+P2+P3)
	END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE FORCE(TIME,WIDTH,NNODE1,NS1,NODE1,PR1,FOR)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      INTEGER  NNODE1,NS1(4)
      REAL TIME,WIDTH,FB,FW,FOR
      REAL NODE1(NNODE1,2),PR1(NNODE1)

	FOR=0.0

!---HORIZONTAL FORCE ON WALL
	DO I=NS1(1)+1,NS1(2)-1
	DX=ABS(NODE1(I,2)-NODE1(I+1,2))
	FOR=FOR+0.5*DX*(PR1(I)+PR1(I+1))
	END DO
	DO I=NS1(4)+1,NNODE1-1
	DX=ABS(NODE1(I,2)-NODE1(I+1,2))
	FOR=FOR-0.5*DX*(PR1(I)+PR1(I+1))
	END DO
	FOR=FOR*WIDTH

WRITE(9,"(5(E15.8,1X))") TIME,FOR

      RETURN
      END
!**********************************************************************
      SUBROUTINE POWER(NNODE1,NELM1,LN1,NS1,NODE1,PHI1,PPHI1,PHIT1,PPHIT1,PR1,&
					  &NNODE2,NELM2,LN2,NS2,NODE2,PHI2,PPHI2,PHIT2,PPHIT2,PR2,&
					  &NNODE3,NELM3,LN3,NS3,NODE3,PHI3,PPHI3,PHIT3,PPHIT3,PR3,&
					  &TIME,GRAV,THO,WIDTH,DEP,BROOT)
!**********************************************************************
     IMPLICIT INTEGER (I-N)
     IMPLICIT REAL    (A-H,O-Z)
	 INTEGER NNODE1,NELM1,LN1(NELM1,2),NS1(4),NNODE2,NELM2,LN2(NELM2,2),NS2(3),NNODE3,NELM3,LN3(NELM3,2),NS3(3)
	 REAL NODE1(NNODE1,2),PHI1(NNODE1),PPHI1(NNODE1),PHIT1(NNODE1),PPHIT1(NNODE1),PR1(NNODE1)
	 REAL NODE2(NNODE2,2),PHI2(NNODE2),PPHI2(NNODE2),PHIT2(NNODE2),PPHIT2(NNODE2),PR2(NNODE2)
	 REAL NODE3(NNODE3,2),PHI3(NNODE3),PPHI3(NNODE3),PHIT3(NNODE3),PPHIT3(NNODE3),PR3(NNODE3)
	 REAL TIME,GRAV,THO,WIDTH,DEP,BROOT(2),POW_W,POW_B,POW,WORK_W,WORK_B,WORK,KE,DKEDT,PE,DPEDT,TE,DEDT

	POW_W=0.0
	WORK_W=0.0
!----CALCULATE THE POWER=SUM(P_WALL*VN) AND WORK=SUM(P_WALL*DDIS) DUE TO TANK
	DO I=NS1(1)+1,NS1(2)-1
	DX=ABS(NODE1(I,2)-NODE1(I+1,2))
	POW_W=POW_W-0.5*DX*(PR1(I)*PPHI1(I)+PR1(I+1)*PPHI1(I+1))
	WORK_W=WORK_W-0.5*DX*(PR1(I)+PR1(I+1))*DDIS
	END DO
	DO I=NS1(4)+1,NNODE1-1
	DX=ABS(NODE1(I,2)-NODE1(I+1,2))
	POW_W=POW_W+0.5*DX*(PR1(I)*PPHI1(I)+PR1(I+1)*PPHI1(I+1))
	WORK_W=WORK_W+0.5*DX*(PR1(I)+PR1(I+1))*DDIS
	END DO
	DO I=NS2(1)+1,NS2(2)-1
	DX=ABS(NODE2(I,2)-NODE2(I+1,2))
	POW_W=POW_W-0.5*DX*(PR2(I)*PPHI2(I)+PR2(I+1)*PPHI2(I+1))
	WORK_W=WORK_W-0.5*DX*(PR2(I)+PR2(I+1))*DDIS
	END DO
	DO I=NS3(3)+1,NNODE3-1
	DX=ABS(NODE3(I,2)-NODE3(I+1,2))
	POW_W=POW_W+0.5*DX*(PR3(I)*PPHI3(I)+PR3(I+1)*PPHI3(I+1))
	WORK_W=WORK_W+0.5*DX*(PR3(I)+PR3(I+1))*DDIS
	END DO

	POW_W=POW_W*WIDTH
	WORK_W=WORK_W*WIDTH

	POW_B=0.0
	WORK_B=0.0
!----CALCULATE THE POWER=SUM(P_WALL*VN) AND WORK=SUM(P_WALL*DDIS) DUE TO BAFFLE
	DO I=NS2(3)+1,NNODE2-1
	DX=ABS(NODE2(I,2)-NODE2(I+1,2))
	POW_B=POW_B+0.5*DX*(PR2(I)*PPHI2(I)+PR2(I+1)*PPHI2(I+1))
	H=SQRT((NODE2(I,1)-NODE2(I+1,1))**2+(NODE2(I,2)-NODE2(I+1,2))**2)
	A1=H*(2.0*PR2(I+1)+PR2(I))/(PR2(I+1)+PR2(I))/3.0
	A2=SQRT((NODE2(I,1)-BROOT(1))**2+(NODE2(I,2)-BROOT(2))**2)
	ARM=A1+A2
	WORK_B=WORK_B+0.5*DX*(PR2(I)+PR2(I+1))*DDIS_B*ARM
	END DO

	DO I=NS2(1)+1,NS3(2)-1
	DX=ABS(NODE3(I,2)-NODE3(I+1,2))
	POW_B=POW_B-0.5*DX*(PR3(I)*PPHI3(I)+PR3(I+1)*PPHI3(I+1))
	H=SQRT((NODE3(I,1)-NODE3(I+1,1))**2+(NODE3(I,2)-NODE3(I+1,2))**2)
	A1=H*(2.0*PR3(I)+PR3(I+1))/(PR3(I)+PR3(I+1))/3.0
	A2=SQRT((NODE3(I+1,1)-BROOT(1))**2+(NODE3(I+1,2)-BROOT(2))**2)
	ARM=A1+A2
	WORK_B=WORK_B-0.5*DX*(PR3(I)+PR3(I+1))*DDIS_B*ARM
	END DO

	POW_B=POW_B*WIDTH
	WORK_B=WORK_B*WIDTH

!---TOTAL POWER AND WORK DONE BY BAFFLE AND TANK
	WORK=WORK_W+WORK_B
	POW=POW_W+POW_B

	KE=0.0
	DKEDT=0.0
!---FLUID KINETIC ENERGY AND DKE/DT
	DO J=1,NELM1
	DX=SQRT((NODE1(LN1(J,1),1)-NODE1(LN1(J,2),1))**2+(NODE1(LN1(J,1),2)-NODE1(LN1(J,2),2))**2)
	KE=KE+(PHI1(LN1(J,1))*PPHI1(LN1(J,1))+PHI1(LN1(J,2))*PPHI1(LN1(J,2)))*DX/2.0
	
	T11=(PHIT1(LN1(J,1))*PPHI1(LN1(J,1))+PHIT1(LN1(J,2))*PPHI1(LN1(J,2)))*DX/2.0
	T12=(PHI1(LN1(J,1))*PPHIT1(LN1(J,1))+PHI1(LN1(J,2))*PPHIT1(LN1(J,2)))*DX/2.0
!	T21=(DP(LN(J,1),1)**2*PPHI(LN(J,1))+DP(LN(J,2),1)**2*PPHI(LN(J,2)))*TEMP/2.0
!	T22=(DP(LN(J,1),1)*PHI(LN(J,1))*PXN(LN(J,1))+DP(LN(J,2),1)*PHI(LN(J,2))*PXN(LN(J,2)))*TEMP/2.0
!	T23=(DP(LN(J,1),2)**2*PPHI(LN(J,1))+DP(LN(J,2),2)**2*PPHI(LN(J,2)))*TEMP/2.0
!	T24=(DP(LN(J,1),2)*PHI(LN(J,1))*PYN(LN(J,1))+DP(LN(J,2),2)*PHI(LN(J,2))*PYN(LN(J,2)))*TEMP/2.0
	DKEDT=DKEDT+T11+T12+T21+T22+T23+T24
	END DO
	KE=KE*0.5*THO*WIDTH
	DKEDT=DKEDT*0.5*THO*WIDTH


	PE=0.0
	DPEDT=0.0
!---FLUID POTENTIAL ENERGY AND DPE/DT (INCLUDING NONINERTIAL + TRANSLATIONAL)

!	DO I=1,NELEM(1)
!	A=NODE(I,:)-ORG
!	ITA1=A(1)*N(1)+A(2)*N(2)
!	X1=A(1)*S(1)+A(2)*S(2)
!	DITA1=DP(I,1)*N(1)+DP(I,2)*N(2)

!	A=NODE(I+1,:)-ORG
!	ITA2=A(1)*N(1)+A(2)*N(2)
!	X2=A(1)*S(1)+A(2)*S(2)
!	DITA2=DP(I+1,1)*N(1)+DP(I+1,2)*N(2)

!	TEMP=ABS(X2-X1) ! BECAUSE THIS IS AN INTEGRAL ALONG THE LOCAL X-DIRECTION

!	PE=PE+0.5*(ITA1**2+ITA2**2)*TEMP/2.0*COS(DIS)+(X1*ITA1+X2*ITA2)*TEMP/2.0*SIN(DIS)

!	V11=(ITA1*DITA1+ITA2*DITA2)*TEMP/2.0*COS(DIS)+0.5*(ITA1**2+ITA2**2)*TEMP/2.0*(-SIN(DIS))*VEL
!	V12=(X1*DITA1+X2*DITA2)*TEMP/2.0*SIN(DIS)+(X1*ITA1+X2*ITA2)*TEMP/2.0*COS(DIS)*VEL
!	V21=(ITA1*DP(LN(I,1),1)+ITA2*DP(LN(I,2),1))*TEMP/2.0*SIN(DIS)
!	V22=(ITA1*DP(LN(I,1),2)+ITA2*DP(LN(I,2),2))*TEMP/2.0*COS(DIS)
!	V23=(X1*DP(LN(I,1),2)+X2*DP(LN(I,2),2))*TEMP/2.0*SIN(DIS)
!	DPEDT=DPEDT+V11+V12+V21+V22+V23
!	END DO

!	PE=PE*THO*GRAV-THO*GRAV*WAMASS*0.5*DEP
!	DPEDT=DPEDT*THO*GRAV*WIDTH

!---TOTAL POWER AND TOTAL MECHANICAL ENERGY
	TE=KE+PE
	DEDT=DKEDT+DPEDT

!-----CALCULATE THE RIGID-BODY (A RECTANGULAR WATER) ENERGY
!	TEMP=WAMASS/DEP
!	INER=WAMASS*(DEP**2+TEMP**2)/12.0+WAMASS*(DEP*0.5)**2
!	RGBT=0.5*THO*INER*VEL**2+THO*GRAV*WAMASS*0.5*DEP*(COS(DIS)-1.0)
!	DRGBT=THO*INER*ACC*VEL-THO*GRAV*WAMASS*0.5*DEP*SIN(DIS)*VEL

	WRITE(51,"(7(E15.8,1X))") TIME,KE,PE,TE,WORK_W,WORK_B,WORK
	WRITE(52,"(7(E15.8,1X))") TIME,DKEDT,DPEDT,DEDT,POW_W,POW_B,POW

      RETURN
      END
!********************************************************************
      SUBROUTINE DOMAIN(NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,L,M,NGA,NFIELD,NNODE,NELM,NELEM(4),NS(4),LN(NELM,2)
	  REAL THO,GRAV,DEP,P_ATM
	  REAL DX,DY,PI2,TEMP,P1,P2,P3
	  REAL NODE(NNODE,2),NORM(NELM,2),JCB(NELM)
	  REAL PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL DNODE(NFIELD,2),DVX(NFIELD),DVY(NFIELD),DPHIT(NFIELD),DPR(NFIELD)
	  REAL KER1(NFIELD,NNODE),KER2(NFIELD,NNODE)
      REAL H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
	  REAL WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
	
	PI2=2.0*ACOS(-1.0)

!----CREATE DOMAIN POINT
	L=1
	DO I=2,NS(1)-1
	K=NS(3)-I+1
	DX=(NODE(K,1)-NODE(I,1))/NELEM(2)
	DY=(NODE(K,2)-NODE(I,2))/NELEM(2)
		DO J=2,NELEM(2)
		DNODE(L,1)=NODE(I,1)+DX*(J-1)
		DNODE(L,2)=NODE(I,2)+DY*(J-1)
		L=L+1
		END DO
	END DO

!----CALCULATE JACOBIAN
!      DO J=1,NELM
!      LENG(J)=((NODE(LN(J,1),1)-NODE(LN(J,2),1))**2+(NODE(LN(J,1),2)-NODE(LN(J,2),2))**2)**0.5
!	  DO L=1,2
!          PXI1(L)=(-0.5)*NODE(LN(J,1),L)+ 0.5*NODE(LN(J,2),L)
!        END DO
!        NX=PXI1(2)
!        NY=PXI1(1)
!        JCB(J)=(NX**2+NY**2)**0.5
!        NORM(J,1)=-NX/JCB(J)
!        NORM(J,2)=NY/JCB(J)
!      END DO

	KER1=0.0
	KER2=0.0
	DVX=0.D0
!----CALCULATE X VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.0
         H(K)=0.0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((YFUNC(M)-DNODE(I,2))**2-(XFUNC(M)-DNODE(I,1))**2)*NORM(J,1)-2.0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,2))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(XFUNC(M)-DNODE(I,1))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVX=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2

	KER1=0.0
	KER2=0.0
	DVY=0.D0
!----CALCULATE Y VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.0
         H(K)=0.0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((XFUNC(M)-DNODE(I,1))**2-(YFUNC(M)-DNODE(I,2))**2)*NORM(J,2)-2.0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,1))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(YFUNC(M)-DNODE(I,2))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVY=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2
	
	KER1=0.0
	KER2=0.0
	DPHIT=0.D0
!----CALCULATE PARTIAL POTENTIAL OVER TIME BY BIE
      DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.0
         H(K)=0.0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(-1)/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*((XFUNC(M)-DNODE(I,1))*NORM(J,1)+(YFUNC(M)-DNODE(I,2))*NORM(J,2))*TEMP
            G(K)=G(K)+LOG(1./((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**0.5)*TEMP	 
		 END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DPHIT=(MATMUL(KER2,PPHIT)-MATMUL(KER1,PHIT))/PI2

!----CALCULATE PRESSURE DISTRIBUTION IN DOMAIN
	DO I=1,NFIELD
	P1=THO*DPHIT(I)
	P2=THO*0.5*(DVX(I)**2+DVY(I)**2)
	P3=THO*GRAV*(DNODE(I,2)-DEP)
	DPR(I)=P_ATM-(P1+P2+P3)
	END DO

	WRITE(11,'(3000(1X,F15.7))') NODE(:,1),DNODE(:,1)
	WRITE(11,'(3000(1X,F15.7))') NODE(:,2),DNODE(:,2)
	WRITE(11,'(3000(1X,F15.7))') DP(:,1),DVX
	WRITE(11,'(3000(1X,F15.7))') DP(:,2),DVY
	WRITE(11,'(3000(1X,F15.7))') PR,DPR

      RETURN
      END
!********************************************************************
      SUBROUTINE GAUSS(WT,RT,NGA)
!********************************************************************
      INTEGER NGA
      REAL   WT(NGA),RT(NGA)

      SELECT CASE(NGA)
       CASE(3)
        WT(1)=0.55555555
        WT(2)=0.88888889
        WT(3)=0.55555555
        RT(1)=0.77459667
        RT(2)=0.0
        RT(3)=-0.77459667
       CASE(4)
        WT(1)=0.65214515
        WT(2)=0.34785484
        WT(3)=0.34785484
        WT(4)=0.65214515
        RT(1)=0.33998104
        RT(2)=0.86113631
        RT(3)=-0.86113631
        RT(4)=-0.33998104
       CASE(5)
        WT(1)=0.23692689
        WT(2)=0.47862867
        WT(3)=0.56888889
        WT(4)=0.47862867
        WT(5)=0.23692689
        RT(1)=0.90617985
        RT(2)=0.53846931
        RT(3)=0.0
        RT(4)=-0.53846931
        RT(5)=-0.90617985
	 CASE(6)
	  WT(1)=0.17132449
	  WT(2)=0.36076157
	  WT(3)=0.46791393
	  WT(4)=0.46791393
	  WT(5)=0.36076157
	  WT(6)=0.17132449
	  RT(1)=0.93246951
	  RT(2)=0.66120938
	  RT(3)=0.23861918
	  RT(4)=-0.23861918
	  RT(5)=-0.66120938
	  RT(6)=-0.9346951
       CASE(8)
        WT(1)=0.10122853
        WT(2)=0.22238103
        WT(3)=0.31370664
        WT(4)=0.36268378
        WT(8)=0.10122853
        WT(7)=0.22238103
        WT(6)=0.31370664
        WT(5)=0.36268378
        RT(1)=0.96028985
        RT(2)=0.79666647
        RT(3)=0.52553240
        RT(4)=0.18343464
        RT(8)=-0.96028985
        RT(7)=-0.79666647
        RT(6)=-0.52553240
        RT(5)=-0.18343464
       CASE(10)
        WT(1)=0.06667134
        WT(2)=0.14945134
        WT(3)=0.21908636
        WT(4)=0.26926671
        WT(5)=0.29552422
        WT(10)=0.06667134
        WT(9)=0.14945134
        WT(8)=0.21908636
        WT(7)=0.26926671
        WT(6)=0.29552422
        RT(1)=0.97390652
        RT(2)=0.86506336
        RT(3)=0.67940956
        RT(4)=0.43339539
        RT(5)=0.14887433
        RT(10)=-0.97390652
        RT(9)=-0.86506336
        RT(8)=-0.67940956
        RT(7)=-0.43339539
        RT(6)=-0.14887433
      END SELECT

      RETURN
      END
!********************************************************************
	SUBROUTINE SIMP(N,X,W,A,B)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N
	REAL A,B,X(10000),W(10000)

	DO I=1,N
	X(I)=A+(B-A)*(I-1)/(N-1)
	W(I)=(B-A)/(N-1)/3.D0
	END DO

	DO I=2,N-1,2
	W(I)=4.D0*W(I)
	END DO

	DO I=3,N-2,2
	W(I)=2.D0*W(I)
	END DO

	RETURN
	END
!********************************************************************
	SUBROUTINE PHIT_CONVERGE(N,P1,P2,E1,E2)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N
	REAL P1(N),P2(N),E1,E2

	E1=MAXVAL(ABS(P1-P2))

	E2=0.0
	DO I=1,N
	E2=(P1(I)-P2(I))**2
	END DO
	E2=SQRT(E2/N)

	RETURN
	END
