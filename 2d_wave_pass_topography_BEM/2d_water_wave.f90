!===============================================================================
!	SOLVE 2D WATER WAVE PROBLEM WITH ONE SIDE WAVEMAKER AND ONESIDE FREE BY BEM
!	WEN-HAUI TSAO 2021 LSU
!===============================================================================
      PROGRAM WATER_WAVE
      IMPLICIT NONE
	  INTEGER I,J,NT,NGA,NTIM,NFIELD,NPL,WTYP,NNODE,NELM,ICON,NITER
	  INTEGER,ALLOCATABLE::NELEM(:),ME(:),NS(:),BUNTYP(:),LN(:,:)
	  REAL DEP,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,DELTTIME,SOR
	  REAL TIME,DIS,VEL,ACC,P_ATM,FOR,DDIS,DTEMP,WC,E1,E2,ETOL
	  REAL,ALLOCATABLE::COOR(:,:),SIDE_L(:)
	  REAL,ALLOCATABLE::NODE(:,:),NORM(:,:),JCB(:),LENG(:),PHI(:),PPHI(:),PHIT(:),PPHIT(:)
	  REAL,ALLOCATABLE::KER1(:,:),KER2(:,:),DP(:,:),DPDS(:),PR(:),DPDT(:),ACCMO(:)
	  REAL,ALLOCATABLE::PHIT_TEMP(:)
	  REAL,ALLOCATABLE::WT(:),RT(:),SHA1(:),SHA2(:),SH(:,:)
	  REAL,ALLOCATABLE::DI(:),VE(:),AC(:)

      OPEN(UNIT=1,FILE='1.IPT',STATUS='OLD')
      OPEN(UNIT=2,FILE='2.IPT',STATUS='OLD')
      OPEN(UNIT=5,FILE='IO.DAT')
      OPEN(UNIT=6,FILE='S.DAT')
      OPEN(UNIT=7,FILE='W.DAT')
      OPEN(UNIT=8,FILE='P.DAT')
      OPEN(UNIT=9,FILE='F.DAT')
      OPEN(UNIT=10,FILE='E.DAT')
      OPEN(UNIT=11,FILE='DOMAIN.DAT')
	  OPEN(UNIT=21,FILE='ERR.DAT')
      OPEN(UNIT=99,FILE='TEST.TXT')

!---TOPOGRAGHY AND WAVE TYPE
	CALL WAVE_GEN(NPL,WTYP)
	ALLOCATE(NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL),COOR(NPL,2),SIDE_L(NPL))
	NELEM=0
	NS=0
	BUNTYP=0
	COOR=0.0
	SIDE_L=0.0

!---INPUT ALL KINDS OF PARAMETERS
	CALL INPUT(NPL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,DEP,NGA,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,NTIM,DELTTIME,SOR,ICON,NITER,ETOL,WC)
	ALLOCATE(LN(NELM,2),NODE(NNODE,2),NORM(NELM,2),JCB(NELM),LENG(NELM),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE))
	ALLOCATE(KER1(NNODE,NNODE),KER2(NNODE,NNODE),DP(NNODE,2),DPDS(NNODE),DPDT(NNODE),PR(NNODE),ACCMO(NNODE))
	ALLOCATE(PHIT_TEMP(NNODE))
	ALLOCATE(WT(NGA),RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA))
	ALLOCATE(DI(NTIM),VE(NTIM),AC(NTIM))
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
	DI=0.0
	VE=0.0
	AC=0.0

!---SHAPE FUNCTION AND MESH
    CALL GAUSS(WT,RT,NGA)
    CALL SHAP(SHA1,SHA2,SH,NGA,RT)
	CALL LENGTH(NPL,COOR,SIDE_L)
	CALL MESH(NPL,NNODE,NELM,NELEM,LN,COOR,SIDE_L,NODE)
    WRITE(*,*) 'PASS MESH'

!---GENERATE WAVES
	SELECT CASE (WTYP)
	CASE(1)
		DO I=1,NTIM
		 DI(I)=AMP*SIN(OMEGA*TIME+PSI)
		 VE(I)=AMP*OMEGA*COS(OMEGA*TIME+PSI)
		 AC(I)=-AMP*OMEGA**2*SIN(OMEGA*TIME+PSI)
		 END DO
	CASE(2)
		NT=OMEGA/DELTTIME
		DO I=1,NT+1
		TIME=(I-1)*DELTTIME
		DI(I)=0.5*AMP-0.5*AMP*COS(ACOS(-1.0)/OMEGA*TIME) !STR/DUR*TIME
		VE(I)=ACOS(-1.0)/OMEGA*0.5*AMP*SIN(ACOS(-1.0)/OMEGA*TIME) !STR/DUR
		AC(I)=(ACOS(-1.0)/OMEGA)**2*0.5*AMP*COS(ACOS(-1.0)/OMEGA*TIME) !0.0
		END DO
		DI(NT+2:NTIM)=DI(NT+1)
	END SELECT

DO NT=1,3 !NTIM
	WRITE(*,*) NT,'TH'
     TIME=(NT-1)*DELTTIME

!---WAVEMAKER DIS, VEL, ACC, PHASE LAG (IN RADIUS)
	 DIS=DI(NT)
     VEL=VE(NT)
     ACC=AC(NT)
	 WRITE(5,"(7(E15.8,1X))") TIME,DIS,VEL,ACC

!---REMESH LOCATION AND POTENTIAL OF A FREE-SURFACE NODE AT CURRENT TIME, I.E. (NT-1)TH TIME-STEP (BUT THE FOLLOWING CALCULATION IS FOR THE NEXT TIME-STEP)
      CALL REMESH(NPL,NNODE,NELEM,NODE,NS,TIME,DEP,AMP)
      WRITE(*,*) 'PASS REMESH'

!---BUILD KERNEL OF ZONE1, ZONE2, ZONE3 AND ASSEMBLE THEM INTO BIG KERNEL
	   CALL KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT)
       WRITE(*,*) 'PASS KERNEL'

!---CALCULATE PRESSURE, FORCE, ENERGY AND POWER ON THE BOUNDARY AND IN THE DOMAIN
!---AT CURRENT TIME, I.E. (NT-1)TH TIME-STEP (BUT THE FOLLOWING CALCULATION IS FOR THE NEXT TIME-STEP)
	  CALL PRESSURE(ICON,TIME,THO,GRAV,DEP,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
	  CALL DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)

!**************IMPLICIT SOLVER FOR PHIT AND PPHIT ON THE ABSORPTION SIDE**************
DO J=1,NITER

	  PHIT_TEMP=PHIT
!---APPLY BC FOR SOLVING FREE-SURFACE PPHI
       CALL BOUND(NPL,NNODE,NELEM,NS,BUNTYP,PHI,PPHI,PHIT,VEL,WC)
	   CALL SOLVE_BACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)

!---FIRST-ORDER TAYLOR SERIES EXPANSION (ALSO GET TANGENTIAL VELOCITY ON ZONE BOUNDARY)
	  CALL TAYES1(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,DEP,GRAV,MU,VEL,WC)

!---APPLY BC FOR SOLVING FREE-SURFACE PPHIT
	  CALL ACCBC(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PPHI,DP,ACCMO)
      CALL BOUNDT(NPL,NNODE,NELM,NS,NELEM,BUNTYP,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
	  CALL SOLVE_BACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)

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
      SUBROUTINE INPUT(NPL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,DEP,NGA,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,NTIM,DELTTIME,SOR,ICON,NITER,ETOL,WC)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER I,J,NPL,NFIELD,NNODE,NELM,NGA,NTIM,ICON,NITER
	  INTEGER NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL)
	  REAL    COOR(NPL,2),DEP,ENDTIME,DELTTIME,SOR,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,ETOL,WC,L0,T
      CHARACTER*2 ID
         ID = '*'
!---NODAL COORDINATES
         CALL HEADLINE(ID,1)
         READ(1,*)  ((COOR(I,J),J=1,2),I=1,NPL)
!---ELEMENT MESH NUMBER
         CALL HEADLINE(ID,1)
        READ(1,*) (NELEM(I),I=1,NPL)
        NNODE = 0
        NELM  = 0
		ME(1) = NELEM(1)
		NS	  = 0
        DO I=1,NPL
         NELM = NELM+NELEM(I)
         NNODE = NNODE+(NELEM(I)+1)
        END DO

		DO I=2,NPL
		 ME(I)=ME(I-1)+NELEM(I)
		END DO

      NS(1)=NELEM(1)+1
	  DO I=2,NPL
      NS(I)=NS(I-1)+NELEM(I)+1
	  END DO

	  NFIELD=(NELEM(1)-1)*(NELEM(2)-1)

!---BOUNDARY TYPE
         CALL HEADLINE(ID,1)
        READ(1,*) (BUNTYP(I),I=1,NPL)
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
		DEP=COOR(NPL,2)

!---CALCULATE WAVE SPEED
T=2.0*ACOS(-1.0)/OMEGA
L0=GRAV*T**2/(2.0*ACOS(-1.0))
WC=SQRT(GRAV*DEP) !L0/T

      RETURN
      END
!**********************************************************************
      SUBROUTINE WAVE_GEN(NPL,WTYP)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER NPL,WTYP
      CHARACTER*2 ID
         ID = '*'
!---NUMBER OF PLANES
         CALL HEADLINE(ID,2)
         READ(2,*) NPL
!---WAVE GENERATION: 1=PERIODIC; 2=SOLITARY
         CALL HEADLINE(ID,2)
         READ(2,*) WTYP

      RETURN
      END
!**********************************************************************
      SUBROUTINE LENGTH(NPL,COOR1,SIDE_L1)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
	  INTEGER NPL
      REAL  COOR1(NPL,2),SIDE_L1(NPL)
      DO I=1,NPL-1
        SIDE_L1(I)=SQRT((COOR1(I+1,1)-COOR1(I,1))**2+(COOR1(I+1,2)-COOR1(I,2))**2)
      END DO
        SIDE_L1(NPL)=SQRT((COOR1(NPL,1)-COOR1(1,1))**2+(COOR1(NPL,2)-COOR1(1,2))**2)

      RETURN
      END
!**********************************************************************
      SUBROUTINE MESH(NPL,NNODE,NELM,NELEM,LN,COOR,LENG,NODE)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      INTEGER NPL,NEL(NPL),NNODE,NELM,NELEM(NPL),LN(NELM,2)
      REAL    SX,SY,NORM,DELT,LENG(NPL),COOR(NPL,2),NODE(NNODE,2)

K=0
DO I=1,NPL-1
	J=NPL-I
    NEL(I) = NELEM(I)+1
    DELT=LENG(J)/NELEM(I)
	SX=COOR(J,1)-COOR(J+1,1)
	SY=COOR(J,2)-COOR(J+1,2)
	NORM=SQRT(SX**2+SY**2)
	SX=SX/NORM
	SY=SY/NORM
    DO L=1,NELEM(I)+1
       NODE(L+K,1)=COOR(J+1,1)+(L-1)*DELT*SX
       NODE(L+K,2)=COOR(J+1,2)+(L-1)*DELT*SY
    END DO
    K=K+NEL(I)
END DO

    NEL(NPL) = NELEM(NPL)+1
    DELT=LENG(NPL)/NELEM(NPL)
	SX=COOR(NPL,1)-COOR(1,1)
	SY=COOR(NPL,2)-COOR(1,2)
	NORM=SQRT(SX**2+SY**2)
	SX=SX/NORM
	SY=SY/NORM
    DO I=1,NELEM(NPL)+1
      NODE(I+K,1)=COOR(1,1)+(I-1)*DELT*SX
      NODE(I+K,2)=COOR(1,2)+(I-1)*DELT*SY
    END DO

!----TO GIVE THE LOCAL ELEMENT NODE NUMBER
      L=1
	  N=1
      DO I=1,NPL
       DO J=1,NELEM(I)
        LN(N,1)=L
        LN(N,2)=L+1
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
      SUBROUTINE REMESH(NPL,NNODE,NELEM,NODE,NS,TIME,DEP,AMP)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL(A-H,O-Z)
      INTEGER NPL,NNODE,NELEM(NPL),NS(NPL)
      REAL TIME,DEP,AMP,LENG,NODE(NNODE,2),SX,SY,NORM

!------ENSURE DUPLICATE POINT ON END NODE OF FREE SURFACE-----
    NODE(NNODE,:)=NODE(1,:)
	NODE(NS(1)+1,:)=NODE(NS(1),:)

!------BOTTOM END NODE GOES WITH FREE SURFACE -----
	NODE(NS(NPL-1),1)=NODE(1,1)
	NODE(NS(2),1)=NODE(NS(1),1)

!----- REMESH FOR ALL PLANES EXCEPT THE FREE SURFACE
DO I=2,NPL
	SX=NODE(NS(I),1)-NODE(NS(I-1),1)
	SY=NODE(NS(I),2)-NODE(NS(I-1),2)
	NORM=SQRT(SX**2+SY**2)
	SX=SX/NORM
	SY=SY/NORM
	DELT=NORM/NELEM(I)
	  K=1
      DO L=NS(I-1)+1,NS(I)
        NODE(L,1)=NODE(NS(I-1),1)+DELT*(K-1)*SX
        NODE(L,2)=NODE(NS(I-1),2)+DELT*(K-1)*SY
		K=K+1
      END DO
END DO

!==========OUTPUT FREE-SURFACE AND WALL AND BAFFLE BOUNDARY==========
WAV_L=NODE(1,2)-DEP
WAV_R=NODE(NELEM(1)+1,2)-DEP
WRITE(6,"(500(E15.8,1X))") NODE(:,1) 
WRITE(6,"(500(E15.8,1X))") NODE(:,2)
!WRITE(7,"(5(E15.8,1X))") TIME,WAV_L,WAV_R,WAV_L/AMP,WAV_R/AMP

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUND(NPL,NNODE,NELEM,NS,BUNTYP,PHI,PPHI,PHIT,VEL,WC)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),NS(NPL),BUNTYP(NPL)
	   REAL R,VEL,WC
       REAL PHI(NNODE),PPHI(NNODE),PHIT(NNODE)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			   PHI(J)=PHI(J)
			   PPHI(J)=0.0
            ELSE
			  IF (I==2)THEN
			   PPHI(J)=-PHIT(J)/WC
			   PHI(J)=0.0
			  ELSE IF (I==NPL)THEN
			   PPHI(J)=-VEL
			   PHI(J)=0.0
			  ELSE
			   PPHI(J)=0.0
			   PHI(J)=0.0
			  END IF
            END IF
          END DO
          N=N+(NELEM(I)+1)
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
       SUBROUTINE SOLVE_BACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!      USING GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION
!======================================================================
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
       REAL    PHI(NNODE),PPHI(NNODE)
       REAL    KER1(NNODE,NNODE),KER2(NNODE,NNODE),GG(NNODE,NNODE)
       REAL    H1(500,500),Q1(500),TEMP(500)
       REAL*8  SUM,A,SIG,G1(500,500),P1(500)

!********** TO MOVE THE KER1 AND KER2**********************************
!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,NPL
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
           DO L=1,NPL
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
       DO I=1,NPL
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
       SUBROUTINE SOLVE(NPL,PHI,PPHI,KER1,KER2,BUNTYP,NNODE,SOR,NELEM)
!**********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NPL,NNODE,BUNTYP(NPL),NELEM(NPL)
       REAL    PHI(NNODE),PPHI(NNODE)
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
       DO I=1,NPL
          DO J=K+N,NELEM(I)+N+1
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
               !P1(J)=(PPHIB(J)-PPHIA(J))+PPHIB(J)
             ELSE
               Q1(J)=PPHI(J)
               !P1(J)=(PHIB(J)-PHIA(J))+PHIB(J)   
            END IF
          END DO
        N=N+NELEM(I)+1
       END DO
!-------------------
       DO I=1,NNODE
          N=0
          DO L=1,NPL
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
       DO I=1,NPL
          DO J=K+N,NELEM(I)+N+1
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
              ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=NELEM(I)+N+1
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
      SUBROUTINE TAYES1(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,DEP,GRAV,MU,VEL,WC)
!********************************************************************
    IMPLICIT INTEGER NONE
    INTEGER I,J,K,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),NS(NPL),LN(NELM,2)
    REAL WC,DEP,GRAV,MU,VEL
	REAL NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2),DPDS(NNODE),PHIT(NNODE),DPDT(NNODE)

	DPDS=0.0
!*********************ON FREE SURFACE*********************
    DO I=1,ME(1)
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

!*********************ON RIGHT WALL*********************
    DO I=ME(1)+1,ME(2)
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

!*********************LEFT WALL*********************
    DO I=ME(NPL-1)+1,ME(NPL)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NPL-1)+1) THEN
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

!*********************BOTTOM*********************
    DO I=ME(2)+1,ME(NPL-1)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)+1) THEN
          DPDS(LN(I,J))=-PPHI(NS(2))
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
          DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(NPL-1)) THEN
          DPDS(LN(I,J))=PPHI(NS(NPL-1)+1)
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
      SUBROUTINE ACCBC(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PPHI,DP,ACCMO)
!********************************************************************
    IMPLICIT INTEGER NONE
    INTEGER I,J,K,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),LN(NELM,2)
    REAL DPDNX,DPDNY,DPNDS
	REAL NORM(NELM,2),JCB(NELM),PPHI(NNODE),DP(NNODE,2),ACCMO(NNODE)

DO K=1,NPL-1
  DO I=ME(K)+1,ME(K+1)
    DO J=1,2
		DPNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
		DPNDX=DPNDS*NORM(I,2)	! PHINN=PHISS=0 FOR LINEAR ELEMENT
		DPNDY=-DPNDS*NORM(I,1)
		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
  END DO
END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUNDT(NPL,NNODE,NELM,NS,NELEM,BUNTYP,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,NELM,NS(NPL),NELEM(NPL),BUNTYP(NPL),LN(NELM,2)
	   REAL R,ACC,WC
       REAL JCB(NELM),PHIT(NNODE),PPHIT(NNODE),DPDS(NNODE),ACCMO(NNODE),DPDSS(NNODE)

	PPHIT=0.0

    DO I=NELEM(1)+1,NELEM(1)+NELEM(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)) THEN
		  DPDSS(LN(I,J))=0.0 ! ASSUME A FLAT GROUND SO VN=0
		ELSE
		  DPDSS(LN(I,J))=0.5*(-DPDS(LN(I,1))+DPDS(LN(I,2)))/JCB(I)
		END IF
    END DO
    END DO

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			 PHIT(J)=PHIT(J)
			 PPHIT(J)=0.0
            ELSE
			  IF (I==2)THEN
				PPHIT(J)=WC*DPDSS(J)
				PHIT(J)=0.0
			  ELSE IF (I==NPL)THEN
			    PPHIT(J)=-ACC-ACCMO(J)
			    PHIT(J)=0.0
			  ELSE
			    PPHIT(J)=0.0-ACCMO(J)
			    PHIT(J)=0.0
			  END IF
            END IF
          END DO
          N=N+(NELEM(I)+1)
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
!********************************************************************
      SUBROUTINE DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,L,M,IL,IR,NPL,NGA,NFIELD,NNODE,NELM,NELEM(NPL),NS(NPL),LN(NELM,2)
	  REAL THO,GRAV,DEP,P_ATM
	  REAL HB,DX,DY,PI2,TEMP,P1,P2,P3
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
	  CALL BWLOC(NODE(I,1),NS(NPL-1)-NS(2),NODE(NS(2)+1:NS(NPL-1),1),NS(2),IL,IR)
	  HB=0.5*(NODE(IL,2)+NODE(IR,2))+0.015 ! keep it a little far away from the boundary
	  DY=-NODE(I,2)/NELEM(2)
		DO J=2,NELEM(2)
		  TEMP=NODE(I,2)+DY*(J-1)
			IF(TEMP>HB)THEN
			DNODE(L,1)=NODE(I,1)
			DNODE(L,2)=NODE(I,2)+DY*(J-1)
			L=L+1
			END IF
		END DO
	END DO

!--- SET A DUMMY NODE TO USELESS DNODE
	DO I=L,NFIELD
	DNODE(I,:)=DNODE(1,:)
	END DO

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
!********************************************************************
	SUBROUTINE BWLOC(PX,N,X,IST,IL,IR)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,IST,IL,IR
	REAL PX,X(N)

	DO I=1,N-1
		IF(X(I)>=PX.AND.X(I+1)<=PX)THEN
		IR=IST+I
		IL=IST+I+1
		GOTO 777
		END IF
	END DO
777 CONTINUE

	RETURN
	END
