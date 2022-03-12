!===============================================================================
!	SOLVE 2D WATER WAVE PROBLEM WITH ONE SIDE WAVEMAKER AND ONESIDE FREE BY BEM
!	WEN-HAUI TSAO 2021 LSU
!===============================================================================
      PROGRAM WATER_WAVE
      IMPLICIT NONE
	  INTEGER I,J,NT,NGA,NTIM,NFIELD,NPL,WTYP,NNODE,NELM,ICON,NITER,NWG,E1LOC,OUTYP
	  INTEGER,ALLOCATABLE::NELEM(:),ME(:),NS(:),BUNTYP(:),LN(:,:)
	  REAL*8 DEP,D_OUT,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,DELTTIME,WGX(10)
	  REAL*8 TIME,DIS,VEL,ACC,P_ATM,FOR,DDIS,DTEMP,WC,E1,E2,ETOL
	  REAL*8,ALLOCATABLE::COOR(:,:),SIDE_L(:)
	  REAL*8,ALLOCATABLE::NODE(:,:),NORM(:,:),JCB(:),LENG(:),PHI(:),PPHI(:),PHIT(:),PPHIT(:)
	  REAL*8,ALLOCATABLE::KER1(:,:),KER2(:,:),DP(:,:),DPDS(:),PR(:),DPDT(:),ACCMO(:)
	  REAL*8,ALLOCATABLE::PHIT_TEMP(:)
	  REAL*8,ALLOCATABLE::WT(:),RT(:),SHA1(:),SHA2(:),SH(:,:)
      REAL*8,ALLOCATABLE::G1(:,:),H1(:,:),EYE(:)
	  REAL*8,ALLOCATABLE::DI(:),VE(:),AC(:)

      OPEN(UNIT=1,FILE='1.ipt',STATUS='OLD')
      OPEN(UNIT=2,FILE='2.ipt',STATUS='OLD')
      OPEN(UNIT=5,FILE='IO.DAT')
      OPEN(UNIT=6,FILE='S.DAT')
      OPEN(UNIT=7,FILE='WG.DAT')
!      OPEN(UNIT=8,FILE='P.DAT')
!      OPEN(UNIT=9,FILE='F.DAT')
!      OPEN(UNIT=10,FILE='E.DAT')
!      OPEN(UNIT=11,FILE='DOMAIN.DAT')
	  OPEN(UNIT=21,FILE='ERR.DAT')
	  OPEN(UNIT=22,FILE='ABORT.TXT')
	  OPEN(UNIT=23,FILE='CFL.DAT')
      OPEN(UNIT=99,FILE='TEST.TXT')

!---TOPOGRAGHY AND WAVE TYPE
	CALL INPUT_2(NPL,WTYP,OUTYP,NWG,WGX)
	ALLOCATE(NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL),COOR(NPL,2),SIDE_L(NPL))
	NELEM=0
	NS=0
	BUNTYP=0
	COOR=0.D0
	SIDE_L=0.D0

!---INPUT ALL KINDS OF PARAMETERS
	CALL INPUT_1(NPL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,DEP,NGA,GRAV,MU,WIDTH,THO,&
			  &AMP,OMEGA,PSI,NTIM,DELTTIME,ICON,NITER,ETOL)
	ALLOCATE(LN(NELM,2),NODE(NNODE,2),NORM(NELM,2),JCB(NELM),LENG(NELM))
	ALLOCATE(PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),PHIT_TEMP(NNODE))
	ALLOCATE(KER1(NNODE,NNODE),KER2(NNODE,NNODE),DP(NNODE,2))
	ALLOCATE(DPDS(NNODE),DPDT(NNODE),PR(NNODE),ACCMO(NNODE))
	ALLOCATE(WT(NGA),RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA))
	ALLOCATE(DI(NTIM),VE(NTIM),AC(NTIM))
    ALLOCATE(G1(NNODE,NNODE),H1(NNODE,NNODE),EYE(NNODE))
	LN=0
	NODE=0.D0
	NORM=0.D0
	JCB=0.D0
	LENG=0.D0
	PHI=0.D0
	PPHI=0.D0
	PHIT=0.D0
	PPHIT=0.D0
	KER1=0.D0
	KER2=0.D0
	DP=0.D0
	DPDS=0.D0
	DPDT=0.D0
	PR=0.D0
	ACCMO=0.D0
	PHIT_TEMP=0.D0
	WT=0.D0
	RT=0.D0
	SHA1=0.D0
	SHA2=0.D0
	SH=0.D0
	DI=0.D0
	VE=0.D0
	AC=0.D0
    G1=0.D0
    H1=0.D0
    EYE=1.D0
    
!---PREPARE IF OUTLET IS A WALL (IF A WALL, NO ITERATION NEEDED)
    IF (OUTYP==0)THEN
        NITER=1
    END IF
        
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
		 TIME=(I-1)*DELTTIME
		 DI(I)=AMP*DCOS(OMEGA*TIME+PSI)-AMP !AMP*SIN(OMEGA*TIME+PSI)
		 VE(I)=-AMP*OMEGA*DSIN(OMEGA*TIME+PSI) !AMP*OMEGA*COS(OMEGA*TIME+PSI)
		 AC(I)=-AMP*OMEGA**2*DCOS(OMEGA*TIME+PSI) !-AMP*OMEGA**2*SIN(OMEGA*TIME+PSI)
		 END DO
	CASE(2)
		NT=OMEGA/DELTTIME
		DO I=1,NT+1
		TIME=(I-1)*DELTTIME
		DI(I)=0.5D0*AMP-0.5D0*AMP*DCOS(DACOS(-1.D0)/OMEGA*TIME) !STR/DUR*TIME
		VE(I)=DACOS(-1.D0)/OMEGA*0.5D0*AMP*DSIN(DACOS(-1.D0)/OMEGA*TIME) !STR/DUR
		AC(I)=(DACOS(-1.D0)/OMEGA)**2*0.5D0*AMP*DCOS(DACOS(-1.D0)/OMEGA*TIME) !0.D0
		END DO
		DI(NT+2:NTIM)=DI(NT+1)
    END SELECT  
    
DO NT=1,NTIM
	WRITE(*,*) NT,'TH'
     TIME=(NT-1)*DELTTIME
!---WAVEMAKER DIS, VEL, ACC, PHASE LAG (IN RADIUS)
	 DIS=DI(NT)
     VEL=VE(NT)
     ACC=AC(NT)
	 WRITE(5,"(7(E15.8,1X))") TIME,DIS,VEL,ACC

!---CALCULATE WAVE SPEED FOR RADIATION CONDITION
	D_OUT=NODE(NS(1),2)-COOR(NPL-2,2)
	CALL WAVE_SPD(GRAV,OMEGA,D_OUT,WC)

!---REMESH LOCATION AND POTENTIAL OF A FREE-SURFACE NODE AT CURRENT (NT-1)TH TIME
!---THE FOLLOWING CALCULATION IS FOR THE NEXT TIME-STEP)
      CALL REMESH(NPL,NNODE,NELEM,NODE,NS,TIME,DEP,AMP,NWG,WGX)
      WRITE(*,*) 'PASS REMESH'

!---BUILD KERNEL OF ZONE1, ZONE2, ZONE3 AND ASSEMBLE THEM INTO BIG KERNEL
	   CALL KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE)
       WRITE(*,*) 'PASS KERNEL'

!---CALCULATE PRESSURE ON THE BOUNDARY
	  CALL PRESSURE(ICON,TIME,THO,GRAV,DEP,NPL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)

!---CALCULATE PRESSURE AND VELOCITY IN THE DOMAIN
!	  CALL DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
!				 &SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)

!**************IMPLICIT SOLVER FOR PHIT AND PPHIT ON THE ABSORPTION SIDE**************
DO J=1,NITER

	  PHIT_TEMP=PHIT
!---APPLY BC FOR SOLVING FREE-SURFACE PPHI
       CALL BOUND(NPL,NNODE,NELEM,NS,BUNTYP,OUTYP,PHI,PPHI,PHIT,VEL,WC)
!	   CALL SOLVE_LAPACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!	   CALL SOLVE_BACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
        CALL SOLVE_LAPACK2_1(NPL,PHI,PPHI,KER1,KER2,G1,H1,NNODE,NELEM,BUNTYP)

        !---FIRST-ORDER TAYLOR SERIES EXPANSION (ALSO GET TANGENTIAL VELOCITY ON ZONE BOUNDARY)
	  CALL TAYES1(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,&
				 &DEP,GRAV,MU,VEL,WC)

!---APPLY BC FOR SOLVING FREE-SURFACE PPHIT
	  CALL ACCBC(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PPHI,DP,ACCMO)
      CALL BOUNDT(NPL,NNODE,NELM,NS,NELEM,BUNTYP,OUTYP,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
!	  CALL SOLVE_LAPACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)
!	  CALL SOLVE_BACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)
      CALL SOLVE_LAPACK2_2(NPL,PHIT,PPHIT,G1,H1,NNODE,NELEM,BUNTYP)
      
IF (OUTYP==1)THEN
CALL CONVERGE(NELEM(2)+1,PHIT_TEMP(NS(1)+1:NS(2)),PHIT(NS(1)+1:NS(2)),E1LOC,E1,E2)
	IF(E1<=ETOL.AND.E2<ETOL)THEN
		WRITE(*,*) 'PASS CONVERGED',J
		WRITE(21,*) TIME,J,E1,E2
		GOTO 205
	ELSE IF(J>=NITER)THEN
		WRITE(22,*) 'CG FAIL',TIME,E1LOC,E1,E2
		WRITE(*,*) 'CONVERGE FAIL'
		STOP
	END IF
END IF

END DO
!**************************************************************************************

205 CONTINUE

!---CHECK THE CFL NUMBER	
!	  CALL COURANT(TIME,DELTTIME,NNODE,NELM,LN,NODE,DP,JCB)

!---SECOND-ORDER TAYLOR SERIES EXPANSION
      CALL TAYES2(PHI,PPHI,PHIT,PPHIT,DPDS,DPDT,DP,NNODE,NELEM,NODE,NELM,NORM,JCB,LN,&
				 &DELTTIME,GRAV,ACC)
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
SUBROUTINE INPUT_1(NPL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,DEP,NGA,GRAV,MU,WIDTH,THO,&
				&AMP,OMEGA,PSI,NTIM,DELTTIME,ICON,NITER,ETOL)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER I,J,NPL,NFIELD,NNODE,NELM,NGA,NTIM,ICON,NITER
	  INTEGER NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL)
	  REAL*8 COOR(NPL,2),DEP,ENDTIME,DELTTIME,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,ETOL
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

	  NFIELD=(NELEM(1)-1)*(NELEM(NPL)-1)

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
         READ(1,*)  ENDTIME,DELTTIME,ICON,NITER,ETOL
		NTIM=ENDTIME/DELTTIME+1
		DEP=COOR(NPL,2)

      RETURN
      END
!**********************************************************************
      SUBROUTINE INPUT_2(NPL,WTYP,OUTYP,NWG,WGX)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER NPL,WTYP,OUTYP,NWG
	  REAL*8 WGX(10)
      CHARACTER*2 ID
         ID = '*'
!---NUMBER OF PLANES
         CALL HEADLINE(ID,2)
         READ(2,*) NPL
!---WAVE GENERATION: 1=PERIODIC; 2=SOLITARY
         CALL HEADLINE(ID,2)
         READ(2,*) WTYP
!---WAVE GENERATION: 0=WALL; 1=RADIATION
         CALL HEADLINE(ID,2)
         READ(2,*) OUTYP
!---WAVE GAUGE NUMBER
         CALL HEADLINE(ID,2)
         READ(2,*) NWG
		 IF (NWG>10)THEN
		 WRITE(*,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 WRITE(22,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 STOP
		 END IF
!---WAVE GAUGE LOCATION
         CALL HEADLINE(ID,2)
         READ(2,*) WGX(1:NWG)

      RETURN
      END
!**********************************************************************
      SUBROUTINE LENGTH(NPL,COOR1,SIDE_L1)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
	  INTEGER NPL
      REAL*8  COOR1(NPL,2),SIDE_L1(NPL)
      DO I=1,NPL-1
        SIDE_L1(I)=DSQRT((COOR1(I+1,1)-COOR1(I,1))**2+(COOR1(I+1,2)-COOR1(I,2))**2)
      END DO
        SIDE_L1(NPL)=DSQRT((COOR1(NPL,1)-COOR1(1,1))**2+(COOR1(NPL,2)-COOR1(1,2))**2)

      RETURN
      END
!**********************************************************************
      SUBROUTINE MESH(NPL,NNODE,NELM,NELEM,LN,COOR,LENG,NODE)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NPL,NEL(NPL),NNODE,NELM,NELEM(NPL),LN(NELM,2)
      REAL*8 SX,SY,NORM,DELT,LENG(NPL),COOR(NPL,2),NODE(NNODE,2)

K=0
DO I=1,NPL-1
	J=NPL-I
    NEL(I) = NELEM(I)+1
    DELT=LENG(J)/NELEM(I)
	SX=COOR(J,1)-COOR(J+1,1)
	SY=COOR(J,2)-COOR(J+1,2)
	NORM=DSQRT(SX**2+SY**2)
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
	NORM=DSQRT(SX**2+SY**2)
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
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NGA
      REAL*8 RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
      DO M=1,NGA
        SHA1(M)=0.5D0*(1-RT(M))
        SHA2(M)=0.5D0*(1+RT(M))

        SH(1,M)=SHA1(M)
        SH(2,M)=SHA2(M)
      END DO
	RETURN
	END 
!**********************************************************************
      SUBROUTINE REMESH(NPL,NNODE,NELEM,NODE,NS,TIME,DEP,AMP,NWG,WGX)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NPL,NNODE,NELEM(NPL),NS(NPL),NWG,IL,IR
      REAL*8 TIME,DEP,AMP,LENG,NODE(NNODE,2),SX,SY,NORM,WGX(10),WGY(10)

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
	NORM=DSQRT(SX**2+SY**2)
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
WRITE(6,"(2000(E15.8,1X))") NODE(:,1) 
WRITE(6,"(2000(E15.8,1X))") NODE(:,2)

!==========OUTPUT WAVE ELEVATION AT THE WAVE GAUGES (IN CM)==========
DO I=1,NWG
CALL BWLOC(-WGX(I),NS(1),-NODE(1:NS(1),1),0,IR,IL)
TEMP=(WGX(I)-NODE(IL,1))/(NODE(IR,1)-NODE(IL,1))
WGY(I)=NODE(IL,2)+TEMP*(NODE(IR,2)-NODE(IL,2))-DEP
END DO
WRITE(7,"(11(E15.8,1X))") TIME,WGY(1:NWG)*100.D0

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUND(NPL,NNODE,NELEM,NS,BUNTYP,OUTYP,PHI,PPHI,PHIT,VEL,WC)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8 (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),NS(NPL),BUNTYP(NPL),OUTYP
	   REAL*8 R,VEL,WC
       REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			   PHI(J)=PHI(J)
			   PPHI(J)=0.D0
            ELSE
			  IF (I==2)THEN
                  IF(OUTYP==0)THEN
                      PPHI(J)=0.D0
                      PHI(J)=0.D0
                  ELSE
                      PPHI(J)=-PHIT(J)/WC
                      PHI(J)=0.D0
                  END IF
			  ELSE IF (I==NPL)THEN
			   PPHI(J)=-VEL
			   PHI(J)=0.D0
			  ELSE
			   PPHI(J)=0.D0
			   PHI(J)=0.D0
			  END IF
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE)
!********************************************************************
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,NNODE,NELM,NGA,LN(NELM,2)
      REAL*8 RD,SIGMA(NNODE),EYE(NNODE)
      REAL*8  KER1(NNODE,NNODE),KER2(NNODE,NNODE)
      REAL*8  NX,NY,H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
      REAL*8  LENG(NELM),NORM(NELM,2),JCB(NELM),NODE(NNODE,2)
      REAL*8  WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)

        KER1=0.D0
        KER2=0.D0
!**** CALCULATE THE JACOBIAN
      DO J=1,NELM
      LENG(J)=DSQRT((NODE(LN(J,1),1)-NODE(LN(J,2),1))**2+(NODE(LN(J,1),2)-NODE(LN(J,2),2))**2)
		DO L=1,2
          PXI1(L)=(-0.5D0)*NODE(LN(J,1),L)+ 0.5D0*NODE(LN(J,2),L)
        END DO
        NX=PXI1(2)
        NY=PXI1(1)
        JCB(J)=DSQRT(NX**2+NY**2)
        NORM(J,1)=-NX/JCB(J)
        NORM(J,2)=NY/JCB(J)
      END DO
      
!$omp parallel do private(I,J,K,M,XFUNC,YFUNC,RD,G,H)
!***THE SURFACE KERNELS
      DO I = 1,NNODE
       DO J=1,NELM
       DO M=1,NGA
          XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
          YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
       END DO
        DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
         RD=DSQRT((NODE(I,1)-NODE(LN(J,K),1))**2+(NODE(I,2)-NODE(LN(J,K),2))**2)
        IF (RD .LE. 0.0000001D0) THEN
        H(K)=0.D0
	    G(K)=LENG(J)/2*(1.5D0-DLOG(LENG(J)))
		ELSE !---NON DSINGULER TERM
         G(K)=0.D0
         DO M=1,NGA
            H(K)=H(K)+(-1.D0)/((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)*&
                    &((XFUNC(M)-NODE(I,1))*NORM(J,1)+(YFUNC(M)-NODE(I,2))*NORM(J,2))&
                    &*JCB(J)*SH(K,M)*WT(M)
            G(K)=G(K)+DLOG(1.D0/((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)**0.5D0)*JCB(J)*SH(K,M)*WT(M)
         END DO
      END IF
         KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
      END DO
      END DO
      END DO
!$omp end parallel do
      
!***DSINGULAR OF KER1
    DO I=1,NNODE
    KER1(I,I)=0.D0
    END DO      
      CALL DGEMM('N','N',NNODE,1,NNODE,1.D0,KER1,NNODE,EYE,NNODE,0.D0,SIGMA,NNODE)
    DO I=1,NNODE
    KER1(I,I)=-SIGMA(I)
    END DO  
!       DO I=1,NNODE
!          SIGMA=0.D0
!          DO J=1,NNODE
!             SIGMA=KER1(I,J)+SIGMA
!          END DO
!          KER1(I,I)=-SIGMA
!       END DO

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
       IMPLICIT REAL*8    (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
       REAL*8    PHI(NNODE),PPHI(NNODE)
       REAL*8    KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8    H1(NNODE,NNODE),Q1(NNODE),TEMP(NNODE)
       REAL*8  SUM,A,SIG,G1(NNODE,NNODE),P1(NNODE)

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
          TEMP(I)=0.D0
          DO J=1,NNODE
          TEMP(I)=TEMP(I)+H1(I,J)*Q1(J)
          END DO
       END DO
!*************GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION*********
      DO I=1,NNODE
        G1(I,NNODE+1)=TEMP(I)
      END DO
      DO I=1,NNODE
         SUM=0.D0
         DO K=I,NNODE
            IF (G1(K,I) .NE. 0) THEN
               IF (K .NE. I) THEN
               IF (G1(I,I) .EQ. 0.D0) THEN
                 WRITE(*,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 WRITE(22,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
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
          IF (SUM .EQ. 0.D0) THEN
          WRITE(*,*) 'NO UNIQUE SOLUTION EXISTS  STOP AT GAUSSELI'
          WRITE(22,*) 'NO UNIQUE SOLUTION EXISTS  STOP AT GAUSSELI'
          STOP
          END IF
      END DO

      P1(NNODE)=G1(NNODE,NNODE+1)/G1(NNODE,NNODE)
      DO I=NNODE-1,1,-1
         SIG=0.D0
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
       SUBROUTINE SOLVE_LAPACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!======================================================================
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8    (A-H,O-Z)
       INTEGER NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8    PHI(NNODE),PPHI(NNODE)
       REAL*8    KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8    H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)

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
          P1(I)=0.D0
          DO J=1,NNODE
          P1(I)=P1(I)+H1(I,J)*Q1(J)
          END DO
       END DO

!*************SOLVE BY CALLING LAPACK*********
CALL DGESV(NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

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
       SUBROUTINE SOLVE_LAPACK2_1(NPL,PHI,PPHI,KER1,KER2,G1,H1,NNODE,NELEM,BUNTYP)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
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
       
!-----MOVE KER1 AND KER2---- 
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

       P1=MATMUL(H1,Q1)
       
!*************SOLVE BY CALLING LAPACK*********
CALL DGETRF(NNODE,NNODE,G1,NNODE,IPIV,INFO)
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)
       
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
       SUBROUTINE SOLVE_LAPACK2_2(NPL,PHI,PPHI,G1,H1,NNODE,NELEM,BUNTYP)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
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

       P1=MATMUL(H1,Q1)

!*************SOLVE BY CALLING LAPACK*********
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

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
!********************************************************************
SUBROUTINE TAYES1(NPL,NNODE,NELM,NELEM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,&
	  &DPDS,DP,PHIT,DPDT,DEP,GRAV,MU,VEL,WC)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),NS(NPL),LN(NELM,2)
    REAL*8 WC,DEP,GRAV,MU,VEL
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2)
	REAL*8 DPDS(NNODE),PHIT(NNODE),DPDT(NNODE)

	DPDS=0.D0
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
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

    DO I=1,NS(1)
	  DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*PHI(I) !-GRAV*NODE(I,2)-MU*PHI(I) !
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
			DPDS(LN(I,J))=0.D0
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=PPHI(LN(I,J))
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************LEFT WALL*********************
    DO I=ME(NPL-1)+1,ME(NPL)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NPL-1)+1) THEN
			DPDS(LN(I,J))=0.D0
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NNODE) THEN
			DPDS(LN(I,J))=DP(1,2)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
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
		  DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
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
    IMPLICIT NONE
    INTEGER I,J,K,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),LN(NELM,2)
    REAL*8 DPNDX,DPNDY,DPNDS
	REAL*8 NORM(NELM,2),JCB(NELM),PPHI(NNODE),DP(NNODE,2),ACCMO(NNODE)

DO K=1,NPL-1
  DO I=ME(K)+1,ME(K+1)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
		DPNDX=DPNDS*NORM(I,2)	! PHINN=PHISS=0 FOR LINEAR ELEMENT
		DPNDY=-DPNDS*NORM(I,1)
		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
  END DO
END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE BOUNDT(NPL,NNODE,NELM,NS,NELEM,BUNTYP,OUTYP,LN,PHIT,PPHIT,DPDS,JCB,WC,ACC,ACCMO)
!********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8    (A-H,O-Z)
       INTEGER NNODE,NELM,NS(NPL),NELEM(NPL),BUNTYP(NPL),LN(NELM,2),OUTYP
	   REAL*8 R,ACC,WC
       REAL*8 JCB(NELM),PHIT(NNODE),PPHIT(NNODE),DPDS(NNODE),ACCMO(NNODE),DPDSS(NNODE)

	PPHIT=0.D0

    DO I=NELEM(1)+1,NELEM(1)+NELEM(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(2)) THEN
		  DPDSS(LN(I,J))=0.D0 ! ASSUME A FLAT GROUND SO VN=0
		ELSE
		  DPDSS(LN(I,J))=0.5D0*(-DPDS(LN(I,1))+DPDS(LN(I,2)))/JCB(I)
		END IF
    END DO
    END DO

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
			 PHIT(J)=PHIT(J)
			 PPHIT(J)=0.D0
            ELSE
			  IF (I==2)THEN
                  IF (OUTYP==0)THEN
                      PPHIT(J)=0.D0
                      PHIT(J)=0.D0 
                  ELSE
                      PPHIT(J)=WC*DPDSS(J)
                      PHIT(J)=0.D0
                  END IF
			  ELSE IF (I==NPL)THEN
			    PPHIT(J)=-ACC-ACCMO(J)
			    PHIT(J)=0.D0
			  ELSE
			    PPHIT(J)=0.D0-ACCMO(J)
			    PHIT(J)=0.D0
			  END IF
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO

      RETURN
      END
!**********************************************************************
SUBROUTINE TAYES2(PHI,PPHI,PHIT,PPHIT,DPDS,DPDT,DP,NNODE,NELEM,NODE,NELM,NORM,&
				 &JCB,LN,DELTTIME,GRAV,ACC)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8    (A-H,O-Z)
      INTEGER  NNODE,NELM,NELEM(4),LN(NELM,2)
      REAL*8 DELTTIME,GRAV,ACC,D2PDT,DPTDS,DPDNDS
      REAL*8 NORM(NELM,2),JCB(NELM)
      REAL*8 NODE(NNODE,2),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE)
	  REAL*8 DP(NNODE,2),DPDS(NNODE),DPDT(NNODE)
      REAL*8 NEW_NODE(NELEM(1)+1,2),NEW_PHI(NELEM(1)+1),D2P(2)

	DO I=1,NELEM(1)
	  DO J=1,2
		IF(LN(I,J) .EQ. 1) THEN
				D2P(1)=-PPHIT(NNODE)
				DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
				DPTDS=(D2P(1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
				D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE IF (LN(I,J) .EQ. NELEM(1)+1) THEN
				D2P(1)=PPHIT(NELEM(1)+2)
				DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
				DPTDS=(D2P(1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
				D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE
			DPTDS=(-0.5D0*PHIT(LN(I,1))+0.5D0*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
			D2P(1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		END IF

		D2PDT=DP(LN(I,J),1)*D2P(1)+DP(LN(I,J),2)*D2P(2)-GRAV*DP(LN(I,J),2)
		NEW_PHI(LN(I,J))=PHI(LN(I,J))+DELTTIME*DPDT(LN(I,J))+D2PDT*DELTTIME**2/2
		NEW_NODE(LN(I,J),1)=NODE(LN(I,J),1)+DP(LN(I,J),1)*DELTTIME+0.5D0*D2P(1)*DELTTIME**2
		NEW_NODE(LN(I,J),2)=NODE(LN(I,J),2)+DP(LN(I,J),2)*DELTTIME+0.5D0*D2P(2)*DELTTIME**2
	  END DO
	END DO

	DO I=1,NELEM(1)+1
	NODE(I,:)=NEW_NODE(I,:)
	PHI(I)=NEW_PHI(I)
	END DO

      RETURN
      END
!********************************************************************
SUBROUTINE PRESSURE(ICON,TIME,THO,GRAV,DEP,NPL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,ICON,NPL,NNODE,NS(NPL)
      REAL*8 TIME,DEP,THO,GRAV,P1,P2,P3,P_ATM
      REAL*8 NODE(NNODE,2),PHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 CP1(NS(1)),CP2(NS(1)),CP3(NS(1)),CP(NS(1))

!----ATOM PRESSURE (BERNOULLI CONSTANT) ON THE FREE SURFACE
	DO I=ICON,ICON ! NS(1) !
	CP1(I)=THO*PHIT(I)
	CP2(I)=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	CP3(I)=THO*GRAV*(NODE(I,2)-DEP)
	CP(I)=CP1(I)+CP2(I)+CP3(I)
	ENDDO
	P_ATM=CP(ICON)

!----PRESSURE ON ZONE 1 BOUNDARY
	DO I=NS(1)+2,NNODE-1
	P1=THO*PHIT(I)
	P2=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP)
	PR(I)=P_ATM-(P1+P2+P3)
	END DO

      RETURN
      END
!********************************************************************
SUBROUTINE DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
				 &SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,L,M,IL,IR,NPL,NGA,NFIELD,NNODE,NELM,NELEM(NPL),NS(NPL),LN(NELM,2)
	  REAL*8 THO,GRAV,DEP,P_ATM
	  REAL*8 HB,DX,DY,PI2,TEMP,P1,P2,P3
	  REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM)
	  REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 DNODE(NFIELD,2),DVX(NFIELD),DVY(NFIELD),DPHIT(NFIELD),DPR(NFIELD)
	  REAL*8 KER1(NFIELD,NNODE),KER2(NFIELD,NNODE)
      REAL*8 H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
	  REAL*8 WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
	
	PI2=2.D0*DACOS(-1.D0)

!----CREATE DOMAIN POINT
	L=1
	DO I=2,NS(1)-1
	  CALL BWLOC(NODE(I,1),NS(NPL-1)-NS(2),NODE(NS(2)+1:NS(NPL-1),1),NS(2),IL,IR)
	  HB=NODE(IL,2)+(NODE(I,1)-NODE(IL,1))/(NODE(IR,1)-NODE(IL,1))*(NODE(IR,2)-NODE(IL,2))+0.01D0 ! keep it a little far away from the boundary
	  DY=-NODE(I,2)/NELEM(NPL)
		DO J=2,NELEM(NPL)
		  TEMP=NODE(I,2)+DY*(J-1)
			IF(TEMP>HB)THEN
			DNODE(L,1)=NODE(I,1)
			DNODE(L,2)=TEMP !NODE(I,2)+DY*(J-1)
			L=L+1
			END IF
		END DO
	END DO

!--- SET A DUMMY NODE TO USELESS DNODE
	DO I=L,NFIELD
	DNODE(I,:)=DNODE(1,:)
	END DO

	KER1=0.D0
	KER2=0.D0
	DVX=0.D0
!----CALCULATE X VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((YFUNC(M)-DNODE(I,2))**2-(XFUNC(M)-DNODE(I,1))**2)*NORM(J,1)-&
					&2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,2))/&
					&((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(XFUNC(M)-DNODE(I,1))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVX=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2

	KER1=0.D0
	KER2=0.D0
	DVY=0.D0
!----CALCULATE Y VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((XFUNC(M)-DNODE(I,1))**2-(YFUNC(M)-DNODE(I,2))**2)*NORM(J,2)-&
					&2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,1))/&
					&((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(YFUNC(M)-DNODE(I,2))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVY=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2
	
	KER1=0.D0
	KER2=0.D0
	DPHIT=0.D0
!----CALCULATE PARTIAL POTENTIAL OVER TIME BY BIE
      DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(-1.D0)/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*&
					&((XFUNC(M)-DNODE(I,1))*NORM(J,1)+(YFUNC(M)-DNODE(I,2))*NORM(J,2))*TEMP
            G(K)=G(K)+DLOG(1.D0/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**0.5D0)*TEMP	 
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
	P2=THO*0.5D0*(DVX(I)**2+DVY(I)**2)
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
      REAL*8   WT(NGA),RT(NGA)

      SELECT CASE(NGA)
       CASE(3)
        WT(1)=0.55555555
        WT(2)=0.88888889
        WT(3)=0.55555555
        RT(1)=0.77459667
        RT(2)=0.D0
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
        RT(3)=0.D0
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
        WT(1)=0.1012285362903763D0
        WT(2)=0.2223810344533745D0
        WT(3)=0.3137066458778873D0
        WT(4)=0.3626837833783620D0
        WT(8)=0.1012285362903763D0
        WT(7)=0.2223810344533745D0
        WT(6)=0.3137066458778873D0
        WT(5)=0.3626837833783620D0
        RT(1)=0.9602898564975363D0
        RT(2)=0.7966664774136267D0
        RT(3)=0.5255324099163290D0
        RT(4)=0.1834346424956498D0
        RT(8)=-0.9602898564975363D0
        RT(7)=-0.7966664774136267D0
        RT(6)=-0.5255324099163290D0
        RT(5)=-0.1834346424956498D0
       CASE(10)
        WT(1)=0.D06667134
        WT(2)=0.14945134
        WT(3)=0.21908636
        WT(4)=0.26926671
        WT(5)=0.29552422
        WT(10)=0.D06667134
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
	SUBROUTINE CONVERGE(N,P1,P2,E1LOC,E1,E2)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,K(1),E1LOC
	REAL*8 P1(N),P2(N),E1,E2

	E1=MAXVAL(DABS(P1-P2))
	E1LOC=MAXLOC(DABS(P1-P2),1)

	E2=0.D0
	DO I=1,N
	E2=E2+(P1(I)-P2(I))**2
	END DO
	E2=DSQRT(E2/N)

	RETURN
	END
!********************************************************************
	SUBROUTINE BWLOC(PX,N,X,IST,IL,IR)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,IST,IL,IR
	REAL*8 PX,X(N)

	DO I=1,N-1
		IF(X(I)>PX.AND.X(I+1)<=PX)THEN
		IR=IST+I
		IL=IST+I+1
		GOTO 777
		END IF
	END DO
777 CONTINUE

	RETURN
	END
!********************************************************************
	SUBROUTINE WAVE_SPD(GRAV,OMEGA,D,C)
!********************************************************************
	IMPLICIT NONE
	INTEGER I
	REAL*8 K,K2,GRAV,OMEGA,D,C,PI,F0,F1
	PI=DACOS(-1.D0)

	K=1.D0
	DO I=1,100
	F0=K*DTANH(K*D)-OMEGA**2/GRAV
	F1=DTANH(K*D)+K-K*(DTANH(K*D)**2)
	K2=K-(F0)/(F1)
		IF((K2-K)/K<=0.000001D0) THEN
		GOTO 717
		END IF
	K=K2
	END DO
	717 CONTINUE

	C=DSQRT(GRAV*DTANH(K*D)/K)

	RETURN
	END
!********************************************************************
SUBROUTINE COURANT(TIME,DELTTIME,NNODE,NELM,LN,NODE,DP,JCB)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,CFLOC,NNODE,NELM,LN(NELM,2)
    REAL*8 TIME,DELTTIME,U,V,VE
	REAL*8 NODE(NNODE,2),DP(NNODE,2),JCB(NELM),CN(NELM),CFL

  DO I=1,NELM
    U=DSQRT(DP(LN(I,1),1)**2+DP(LN(I,1),2)**2)
    V=DSQRT(DP(LN(I,2),1)**2+DP(LN(I,2),2)**2)
    VE=MAX(U,V)
    CN(I)=0.5D0*VE*DELTTIME/JCB(I)
  END DO
  CFL=MAXVAL(CN)
  CFLOC=MAXLOC(CN,1)
  
	WRITE(23,*) TIME,CFL

    IF (CFL>=250.D0)THEN
    WRITE(22,*) TIME,"CFL=",CFL,"@ ELEMENT",CFLOC
    STOP
    END IF    

	RETURN
	END