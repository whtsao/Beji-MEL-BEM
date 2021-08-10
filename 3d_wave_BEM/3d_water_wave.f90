!===============================================================================
!	SOLVE 3D WATER WAVE PROBLEM WITH ONE SIDE WAVEMAKER AND ONESIDE FREE BY BEM
!	WEN-HAUI TSAO 2021 LSU
!===============================================================================
      PROGRAM WATER_WAVE
      IMPLICIT NONE
      INTEGER I,J,NT,NFIELD,NNODE,ENUM,NTIM,ICON,NITER
	  INTEGER NELEM(6,2),NEL(6),NE(6),ME(6),BUNTYP(6)
	  INTEGER,ALLOCATABLE::LN(:,:),CBC(:,:),ACBC(:,:),DN(:,:)
	  REAL PI,QPI,EPI
	  REAL TIME,DIS,VEL,ACC,E1,E2,GRAV,MU,THO,AMP,OMEGA,PSI,DELTTIME,SOR,ETOL,DEP,WC,P_ATM
      REAL COOR(8,3),LENG(12),W2(7),RT1(7),RT2(7),SHA1(7),SHA2(7),SHA3(7),SHA4(7),SH(4,7)
	  REAL,ALLOCATABLE::NODE(:,:),PHI(:),PPHI(:),PHIT(:),PPHIT(:),PHIT_TEMP(:),KER1(:,:),KER2(:,:)
	  REAL,ALLOCATABLE::VELX(:),VELY(:),VELZ(:),DT(:),ACCMO(:),DPDNN(:),PR(:)
	  REAL,ALLOCATABLE::X2ND(:),Y2ND(:),Z2ND(:),PHI2ND(:)
	  REAL,ALLOCATABLE::NORMX(:,:),NORMY(:,:),NORMZ(:,:)
	  REAL,ALLOCATABLE::RM1(:,:),RM2(:,:),THI(:,:)
      REAL,ALLOCATABLE::S1X(:,:),S1Y(:,:),S1Z(:,:)
      REAL,ALLOCATABLE::S2X(:,:),S2Y(:,:),S2Z(:,:)
      REAL,ALLOCATABLE::DPDS1(:,:),DPDS2(:,:)
	  REAL,ALLOCATABLE::NX(:,:),NY(:,:),NZ(:,:),JCB(:,:)
	  REAL XI1(4),XI2(4)
	  DATA XI1/-1,1,1,-1/
      DATA XI2/-1,-1,1,1/
	  PI = ACOS(-1.0)
      QPI= 1./4/PI
      EPI= PI/8

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

!**********************************************************************
      CALL INPUT(NFIELD,NNODE,ENUM,NELEM,NEL,NE,ME,BUNTYP,NTIM,ICON,NITER,COOR,GRAV,MU,THO,AMP,OMEGA,PSI,DELTTIME,SOR,ETOL,DEP,WC)
      WRITE(*,*) 'PASS INPUT'
	  ALLOCATE(LN(ENUM,4),CBC(NEL(1),2),ACBC(NNODE,2),DN(2,(NELEM(1,1)-1)*(NELEM(1,2)-1)))
	  ALLOCATE(NODE(NNODE,3),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),PHIT_TEMP(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE))
	  ALLOCATE(VELX(NNODE),VELY(NNODE),VELZ(NNODE),DT(NNODE),ACCMO(NNODE),DPDNN(NNODE),PR(NNODE))
      ALLOCATE(X2ND(NEL(1)),Y2ND(NEL(1)),Z2ND(NEL(1)),PHI2ND(NEL(1)))
	  ALLOCATE(NORMX(ME(1),4),NORMY(ME(1),4),NORMZ(ME(1),4))
	  ALLOCATE(RM1(ME(1),4),RM2(ME(1),4),THI(ME(1),4))
      ALLOCATE(S1X(ME(1),4),S1Y(ME(1),4),S1Z(ME(1),4))
      ALLOCATE(S2X(ME(1),4),S2Y(ME(1),4),S2Z(ME(1),4))
      ALLOCATE(DPDS1(ME(1),4),DPDS2(ME(1),4))
	  ALLOCATE(NX(ENUM,7),NY(ENUM,7),NZ(ENUM,7),JCB(ENUM,7))
	  LN=0
	  CBC=0
	  ACBC=0
	  DN=0
	  NODE=0.0
	  PHI=0.0
	  PPHI=0.0
	  PHIT=0.0
	  PPHIT=0.0
	  PHIT_TEMP=0.0
	  KER1=0.0
	  KER2=0.0
	  VELX=0.0
	  VELY=0.0
	  VELZ=0.0
	  DT=0.0
	  ACCMO=0.0
	  DPDNN=0.0
	  PR=0.0
	  X2ND=0.0
	  Y2ND=0.0
	  Z2ND=0.0
	  PHI2ND=0.0
	  RM1=0.0
	  RM2=0.0
	  THI=0.0
	  S1X=0.0
	  S1Y=0.0
	  S1Z=0.0
      S2X=0.0
	  S2Y=0.0
	  S2Z=0.0
      DPDS1=0.0
	  DPDS2=0.0
	  NX=0.0
	  NY=0.0
	  NZ=0.0
	  JCB=0.0

!**********************************************************************
      CALL GAUSS(W2,RT1,RT2)
      CALL LENGTH(COOR,LENG)
      CALL SHAP(SHA1,SHA2,SHA3,SHA4,SH,RT1,RT2)

!**********************************************************************
      CALL MESH(NNODE,ENUM,NELEM,NEL,NE,LN,COOR,LENG,NODE,CBC,ACBC,DN)
      WRITE(*,*) 'PASS MESH'

DO NT=1,NTIM
	TIME=(NT-1)*DELTTIME
	WRITE(*,*) NT,'TH'

	DIS=AMP*SIN(OMEGA*TIME+PSI) !-AMP_B*COS(OMEGA_B*TIME+PSI) !
	VEL=AMP*OMEGA*COS(OMEGA*TIME+PSI) !AMP_B*OMEGA_B*SIN(OMEGA_B*TIME+PSI) !
	ACC=-AMP*OMEGA**2*SIN(OMEGA*TIME+PSI) !AMP_B*OMEGA_B**2*COS(OMEGA_B*TIME+PSI) !

!**********************************************************************
	CALL REMESH(NNODE,NELEM,NEL,NODE,COOR,PHI,PHIT,VELX,VELY,VELZ,DT,X2ND,Y2ND,Z2ND,PHI2ND,TIME,DELTTIME)
	WRITE(*,*) 'PASS REMESH'

    CALL KERNEL(KER1,KER2,NNODE,ENUM,NELEM,ME,LN,NODE,SHA1,SHA2,SHA3,SHA4,SH,W2,RT1,RT2,NX,NY,NZ,JCB,QPI,EPI)
	WRITE(*,*) 'PASS KERNEL'

!**********************************************************************
	CALL PRESSURE(ICON,THO,GRAV,DEP,NNODE,NE,NODE,PHIT,VELX,VELY,VELZ,PR,P_ATM)
	CALL DOMAIN(NFIELD,NNODE,ENUM,NELEM,LN,DN,THO,GRAV,DEP,P_ATM,QPI,&
			   &NODE,NX,NY,NZ,JCB,PHI,PPHI,PHIT,PPHIT,VELX,VELY,VELZ,&
			   &PR,SH,W2,SHA1,SHA2,SHA3,SHA4)
!**************IMPLICIT SOLVER FOR PHIT AND PPHIT ON THE ABSORPTION SIDE**************
	DO J=1,NITER

		PHIT_TEMP=PHIT
	!---APPLY BC FOR SOLVING FREE-SURFACE PPHI
		CALL BOUND(NNODE,NE,NODE,PHI,PPHI,PHIT,VEL,WC)
		CALL SOLVE(NNODE,NEL,BUNTYP,KER1,KER2,PHI,PPHI,SOR)

	!---FIRST-ORDER TAYLOR SERIES EXPANSION
		CALL TAYES1(NNODE,ENUM,LN,ME,NE,CBC,NODE,PHI,PPHI,PHIT,VELX,VELY,VELZ,&
				   &NORMX,NORMY,NORMZ,RM1,RM2,THI,S1X,S1Y,S1Z,S2X,S2Y,S2Z,&
				   &DPDS1,DPDS2,DT,GRAV,DEP,XI1,XI2)

	!---PREPARE CONVECTIVE TERM (ALSO GET VELOCITY)
		CALL ACCBC(NNODE,ENUM,NELEM,LN,NEL,ME,ACBC,NODE,PHI,PPHI,VELX,VELY,VELZ,ACCMO,DPDNN,VEL)

	!---APPLY BC FOR SOLVING FREE-SURFACE PPHIT
		CALL BOUNDT(NNODE,NE,PHIT,PPHIT,ACCMO,DPDNN,ACC,WC)
		CALL SOLVE(NNODE,NEL,BUNTYP,KER1,KER2,PHIT,PPHIT,SOR)

	!---CHECK WHEATHER PHIT ON THE RADIATION PLANE CONVERGE OR NOT
		CALL PHIT_CONVERGE(NE(5)-NE(1),PHIT_TEMP(NE(1)+1:NE(5)),PHIT(NE(1)+1:NE(5)),E1,E2)
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
      CALL TAYES2(NNODE,ENUM,LN,ME,NEL,PPHI,PHIT,PPHIT,VELX,VELY,VELZ,&
				 &NORMX,NORMY,NORMZ,S1X,S1Y,S1Z,S2X,S2Y,S2Z,RM1,RM2,THI,&
				 &DPDS1,DPDS2,X2ND,Y2ND,Z2ND,PHI2ND,GRAV)
      WRITE(*,*) 'PASS TAYES2ND'

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
      SUBROUTINE INPUT(NFIELD,NNODE,ENUM,NELEM,NEL,NE,ME,BUNTYP,NTIM,ICON,NITER,COOR,GRAV,MU,THO,AMP,OMEGA,PSI,DELTTIME,SOR,ETOL,DEP,WC)
!**********************************************************************
!     COOR    :THE EIGHT COORDINATE OF THE TANK (X,Y,Z)
!     NELEM   :THE MASH NUMBERS (X*Y OR Y*Z...) OF THE SIX SURFACE
!     BUNTYP  :THE BOUNDARY TYPE OF THE SIX SURFACE
!     BUNVAL  :THE BOUNDARY VALUE OF THE SIX SURFACES
!     NNODE   :THE TOTAL MASH NODES
!======================================================================
      IMPLICIT NONE
      INTEGER I,J,IREAD,NFIELD,NNODE,ENUM,NELEM(6,2),NEL(6),NE(6),ME(6),BUNTYP(6),NTIM,ICON,NITER
	  REAL	  GRAV,MU,THO,AMP,OMEGA,PSI,ENDTIME,DELTTIME,SOR,ETOL,DEP,WC,L0,T
      REAL    COOR(8,3)
      CHARACTER*2 ID
         ID = '*'
         IREAD = 1 
!     NODAL COORDINATES
         CALL HEADLINE(ID,IREAD)
         READ(IREAD,*)  ((COOR(I,J),J=1,3),I=1,8)
		 DEP=COOR(5,3)

!     ELEMENT MESH NUMBER
        CALL HEADLINE(ID,IREAD)
         READ(IREAD,*) ((NELEM(I,J),J=1,2),I=1,6)
      NNODE = 0
      ENUM  = 0
      DO I=1,6
         ENUM  = ENUM+NELEM(I,1)*NELEM(I,2)				! TOTAL ELEMENT NUMBER
         NNODE = NNODE+(NELEM(I,1)+1)*(NELEM(I,2)+1)	! TOTAL NODE NUMBER
      END DO

	  NFIELD=(NELEM(1,1)-1)*(NELEM(1,2)-1)*(NELEM(2,2)-1)

      NEL(1) = (NELEM(1,1)+1)*(NELEM(1,2)+1)			! NODE NUMBER OF EACH PLANE
      NEL(2) = (NELEM(2,1)+1)*(NELEM(2,2)+1)
      NEL(3) = (NELEM(3,1)+1)*(NELEM(3,2)+1)
      NEL(4) = (NELEM(4,1)+1)*(NELEM(4,2)+1)
      NEL(5) = (NELEM(5,1)+1)*(NELEM(5,2)+1)
      NEL(6) = (NELEM(6,1)+1)*(NELEM(6,2)+1)

	  NE(1)=NEL(1)
	  DO I=2,6
	  NE(I)=NE(I-1)+NEL(I)								! ACCUMULATE NODE NUMBER OF EACH PLANE
	  END DO

      ME(1)=NELEM(1,1)*NELEM(1,2)						! ACCUMULATE ELEMENT NUMBER OF EACH PLANE
      ME(2)=NELEM(2,1)*NELEM(2,2)+ME(1)
      ME(3)=NELEM(3,1)*NELEM(3,2)+ME(2)
      ME(4)=NELEM(4,1)*NELEM(4,2)+ME(3)
      ME(5)=NELEM(5,1)*NELEM(5,2)+ME(4)
      ME(6)=NELEM(6,1)*NELEM(6,2)+ME(5)

!     BOUNDARY CONDITION
        CALL HEADLINE(ID,IREAD)
        READ(IREAD,*) (BUNTYP(I),I=1,6)

!	  INPUT THE GRAV,MU,THO
         CALL HEADLINE(ID,IREAD)
         READ(IREAD,*)  GRAV,MU,THO

!	  AMP(rad),OMEGA(rad/s),PSI(rad) of wavemaker
		 CALL HEADLINE(ID,IREAD)
         READ(IREAD,*)  AMP,OMEGA,PSI

!     ENDTIME,DELTTIME,SOR,ICON,NITER,ETOL
         CALL HEADLINE(ID,IREAD)
         READ(IREAD,*)  ENDTIME,DELTTIME,SOR,ICON,NITER,ETOL
		NTIM=ENDTIME/DELTTIME+1

!---CALCULATE WAVE SPEED
	T=2.0*ACOS(-1.0)/OMEGA
	L0=GRAV*T**2/(2.0*ACOS(-1.0))
	WC=L0/T

      RETURN
      END
!**********************************************************************
      SUBROUTINE LENGTH(COOR,LENG)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      REAL  COOR(8,3),LENG(12)
      DO I=1,3
        LENG(I)=SQRT((COOR(I+1,1)-COOR(I,1))**2+(COOR(I+1,2)-COOR(I,2))**2+(COOR(I+1,3)-COOR(I,3))**2)
      END DO
        LENG(4)=SQRT((COOR(4,1)-COOR(1,1))**2+(COOR(4,2)-COOR(1,2))**2+(COOR(4,3)-COOR(1,3))**2)
      DO I=1,4
        LENG(I+4)=SQRT((COOR(I+4,1)-COOR(I,1))**2+(COOR(I+4,2)-COOR(I,2))**2+(COOR(I+4,3)-COOR(I,3))**2)
      END DO
      DO I=1,3
        LENG(I+8)=SQRT((COOR(I+5,1)-COOR(I+4,1))**2+(COOR(I+5,2)-COOR(I+4,2))**2+(COOR(I+5,3)-COOR(I+4,3))**2)
      END DO
      LENG(12)=SQRT((COOR(8,1)-COOR(5,1))**2+(COOR(8,2)-COOR(5,2))**2+(COOR(8,3)-COOR(5,3))**2)
      RETURN
      END
!**********************************************************************
      SUBROUTINE MESH(NNODE,ENUM,NELEM,NEL,NE,LN,COOR,LENG,NODE,CBC,ACBC,DN)
!**********************************************************************
!     NODE   :THE COORDINATE OF EVERY MASH NODE
!     NORM   :THE NORMAL VECTOR OF THE ELEMENT
      IMPLICIT NONE
	  INTEGER I,J,K,L,M,MS
      INTEGER NNODE,ENUM,NELEM(6,2),NEL(6),NE(6),LN(ENUM,4),CBC(NEL(1),2),ACBC(NNODE,2),DN(2,(NELEM(1,1)-1)*(NELEM(1,2)-1))
      REAL    R,DELT1,DELT2,COOR(8,3),LENG(12),NODE(NNODE,3)

!**** THE 1ST SURFACE (X-Y,NORMAL : Z)********************************
      DELT1=LENG(1)/NELEM(1,1)
      DELT2=LENG(2)/NELEM(1,2)
      K=1
      L=1
      DO I=1,NELEM(1,2)+1
         NODE(K,2)=COOR(5,2)+DELT2*(I-1) 
         DO J=1,NELEM(1,1)+1
            NODE(L,1)=COOR(5,1)+DELT1*(J-1)
            NODE(L,2)=NODE(K,2)
            NODE(L,3)=COOR(5,3)
            L=L+1
         END DO
         K=K+NELEM(1,1)+1
      END DO
!**** THE 2ND SURFACE (X-Z,NORMAL : -Y)********************************
      DELT1=LENG(1)/NELEM(2,1)
      DELT2=LENG(5)/NELEM(2,2)
      K=1+NEL(1)
      L=1+NEL(1)
      DO I=1,NELEM(2,2)+1
         NODE(K,3)=COOR(1,3)+DELT2*(I-1) 
         DO J=1,NELEM(2,1)+1
            NODE(L,1)=COOR(1,1)+DELT1*(J-1)
            NODE(L,2)=COOR(1,2)
            NODE(L,3)=NODE(K,3)
            L=L+1
         END DO
         K=K+NELEM(2,1)+1
      END DO
!**** THE 3RD SURFACE (Y-Z,NORMAL : X)********************************
      DELT1=LENG(2)/NELEM(3,1)
      DELT2=LENG(6)/NELEM(3,2)
      K=1+NEL(1)+NEL(2)
      L=1+NEL(1)+NEL(2)
      DO I=1,NELEM(3,2)+1
         NODE(K,3)=COOR(2,3)+DELT2*(I-1) 
         DO J=1,NELEM(3,1)+1
            NODE(L,1)=COOR(2,1)
            NODE(L,2)=COOR(2,2)+DELT1*(J-1)
            NODE(L,3)=NODE(K,3)
            L=L+1
         END DO
         K=K+NELEM(3,1)+1
      END DO
!**** THE 4TH SURFACE (X-Z,NORMAL : Y)********************************
      DELT1=LENG(3)/NELEM(4,1)
      DELT2=LENG(8)/NELEM(4,2)
      K=1+NEL(1)+NEL(2)+NEL(3)
      L=1+NEL(1)+NEL(2)+NEL(3)
      DO I=1,NELEM(4,2)+1
         NODE(K,3)=COOR(4,3)+DELT2*(I-1) 
         DO J=1,NELEM(4,1)+1
            NODE(L,1)=COOR(4,1)+DELT1*(J-1)
            NODE(L,2)=COOR(4,2)
            NODE(L,3)=NODE(K,3)
            L=L+1
         END DO
         K=K+NELEM(4,1)+1
      END DO
!**** THE 5TH SURFACE (Y-Z,NORMAL : -X)********************************
      DELT1=LENG(4)/NELEM(5,1)
      DELT2=LENG(8)/NELEM(5,2)
      K=1+NEL(1)+NEL(2)+NEL(3)+NEL(4)
      L=1+NEL(1)+NEL(2)+NEL(3)+NEL(4)
      DO I=1,NELEM(5,2)+1
         NODE(K,3)=COOR(1,3)+DELT2*(I-1) 
         DO J=1,NELEM(5,1)+1
            NODE(L,1)=COOR(1,1)
            NODE(L,2)=COOR(1,1)+DELT1*(J-1)
            NODE(L,3)=NODE(K,3)
            L=L+1
         END DO
         K=K+NELEM(5,1)+1
      END DO
!**** THE 6TH SURFACE (Y-Z,NORMAL : -Z)********************************
      DELT1=LENG(9)/NELEM(6,1)
      DELT2=LENG(12)/NELEM(6,2)
      K=1+NEL(1)+NEL(2)+NEL(3)+NEL(4)+NEL(5)
      L=1+NEL(1)+NEL(2)+NEL(3)+NEL(4)+NEL(5)
      DO I=1,NELEM(6,2)+1
         NODE(K,2)=COOR(1,2)+DELT2*(I-1) 
         DO J=1,NELEM(6,1)+1
            NODE(L,1)=COOR(1,1)+DELT1*(J-1)
            NODE(L,2)=NODE(K,2)
            NODE(L,3)=COOR(1,3)
            L=L+1
         END DO
         K=K+NELEM(6,1)+1
      END DO

!*****TO GIVE THE LOCAL ELEMENT NODE NUMBER
      L =1
      MS=1
      DO K=1,6
      M=MS
      DO I=1,NELEM(K,2)
         DO J=1,NELEM(K,1)
            LN(L,1)=M
            LN(L,2)=M+1
            LN(L,3)=LN(L,2)+(NELEM(K,1)+1)
            LN(L,4)=LN(L,1)+(NELEM(K,1)+1)
            L=L+1
            M=M+1
         END DO
            M=M+1
      END DO
            MS=MS+NEL(K)
      END DO

!**********DEFINE CONTINUOUS NODE
	DO I=1,NEL(1)
		DO J=NEL(1)+1,NNODE
		R=SQRT((NODE(I,1)-NODE(J,1))**2+(NODE(I,2)-NODE(J,2))**2+(NODE(I,3)-NODE(J,3))**2)
			IF (R<=0.000001.AND.NODE(I,1)==COOR(1,1))THEN
			CBC(I,1)=1
			CBC(I,2)=J
			ACBC(J,1)=1
			ACBC(J,2)=I
			ELSE IF (R<=0.000001.AND.NODE(I,2)==COOR(1,2))THEN
			CBC(I,1)=2
			CBC(I,2)=J
			ACBC(J,1)=2
			ACBC(J,2)=I
			ELSE IF (R<=0.000001.AND.NODE(I,1)==COOR(3,1))THEN
			CBC(I,1)=3
			CBC(I,2)=J
			ACBC(J,1)=3
			ACBC(J,2)=I
			ELSE IF (R<=0.000001.AND.NODE(I,2)==COOR(3,2))THEN
			CBC(I,1)=4
			CBC(I,2)=J
			ACBC(J,1)=4
			ACBC(J,2)=I
			END IF
		END DO
	END DO

!**********CREATE NODE NUMBER FOR DOMAIN MESH
	K=1
	DO I=1,NE(1)
		IF ((NODE(I,1).NE.COOR(1,1)) .AND. (NODE(I,1).NE.COOR(3,1)) .AND. (NODE(I,2).NE.COOR(1,2)) .AND. (NODE(I,2).NE.COOR(3,2)))THEN
		DN(1,K)=I
		K=K+1
		END IF
	END DO

	K=1
	DO I=NE(5)+1,NNODE
		IF ((NODE(I,1).NE.COOR(1,1)) .AND. (NODE(I,1).NE.COOR(3,1)) .AND. (NODE(I,2).NE.COOR(1,2)) .AND. (NODE(I,2).NE.COOR(3,2)))THEN
		DN(2,K)=I
		K=K+1
		END IF
	END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE SHAP(SHA1,SHA2,SHA3,SHA4,SH,RT1,RT2)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL (A-H,O-Z)
      REAL RT1(7),RT2(7)
      REAL SHA1(7),SHA2(7),SHA3(7),SHA4(7),SH(4,7)
      DO M=1,7
        SHA1(M)=0.25*(1-RT1(M))*(1-RT2(M))
        SHA2(M)=0.25*(1+RT1(M))*(1-RT2(M))
        SHA3(M)=0.25*(1+RT1(M))*(1+RT2(M))
        SHA4(M)=0.25*(1-RT1(M))*(1+RT2(M))

        SH(1,M)=SHA1(M)
        SH(2,M)=SHA2(M)
        SH(3,M)=SHA3(M)
        SH(4,M)=SHA4(M)
      END DO
	RETURN
	END
!**********************************************************************
      SUBROUTINE REMESH(NNODE,NELEM,NEL,NODE,COOR,PHI,PHIT,VELX,VELY,VELZ,DT,X2ND,Y2ND,Z2ND,PHI2ND,TIME,DELTTIME)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL(A-H,O-Z)
      INTEGER NNODE,NELEM(6,2),NEL(6)
	  REAL TIME,DELTTIME,DELT(100)
      REAL NODE(NNODE,3),COOR(8,3),PHI(NNODE),PHIT(NNODE)
      REAL VELX(NNODE),VELY(NNODE),VELZ(NNODE),DT(NNODE)
      REAL X2ND(NEL(1)),Y2ND(NEL(1)),Z2ND(NEL(1)),PHI2ND(NEL(1))

!*****THE NEXT TIME STEP
      DO I=1,NEL(1)
       PHI(I)=PHI(I)+DELTTIME*DT(I)+0.5*PHI2ND(I)*DELTTIME**2
       NODE(I,1)=NODE(I,1)+VELX(I)*DELTTIME+0.5*X2ND(I)*DELTTIME**2
       NODE(I,2)=NODE(I,2)+VELY(I)*DELTTIME+0.5*Y2ND(I)*DELTTIME**2
       NODE(I,3)=NODE(I,3)+VELZ(I)*DELTTIME+0.5*Z2ND(I)*DELTTIME**2
      END DO

! LET THE DOUBLE POINT HAVE THE SAME COOINATE
!   THE FREE S. AND 2ND SURFACE
      NI=NEL(1)+(NELEM(2,1)+1)*NELEM(2,2)
      DO I=1,NELEM(1,1)+1
       PHI(I+NI)=PHI(I)
       PHIT(I+NI)=PHIT(I)
       DO J=1,3
         NODE(I+NI,J)=NODE(I,J)
       END DO
      END DO
!   THE FREE S. AND 3RD SURFACE
      NI=NEL(1)+NEL(2)+(NELEM(3,1)+1)*NELEM(3,2)
       DO I=1,NELEM(1,2)+1
        PHI(I+NI)=PHI((NELEM(1,1)+1)*I)
        PHIT(I+NI)=PHIT((NELEM(1,1)+1)*I)
        DO J=1,3
         NODE(I+NI,J)=NODE((NELEM(1,1)+1)*I,J)
        END DO
       END DO
!  THE FREE S. AND 4TH SURFACE
      NI=NEL(1)+NEL(2)+NEL(3)+(NELEM(4,1)+1)*NELEM(4,2)
      NII=(NELEM(1,1)+1)*NELEM(1,2)
      DO I=1,NELEM(1,1)+1
      PHI(I+NI)=PHI(I+NII)
      PHIT(I+NI)=PHIT(I+NII)
       DO J=1,3
        NODE(I+NI,J)=NODE(I+NII,J)
       END DO
      END DO
!  THE FREE S. AND 5TH SURFACE
      NI=NEL(1)+NEL(2)+NEL(3)+NEL(4)+(NELEM(5,1)+1)*NELEM(5,2)
       DO I=1,NELEM(1,2)+1
         NII=(NELEM(1,1)+1)*(I-1)+1
         PHI(I+NI)=PHI(NII)
         PHIT(I+NI)=PHIT(NII)
        DO J=1,3
         NODE(I+NI,J)=NODE(NII,J)
        END DO
       END DO

!   REMESH THE BOUNDARY ELEMENTS

!**** THE 2ND SURFACE (X-Z,NORMAL : -Y)********************************
      L=1+NEL(1)

      NI=(NELEM(2,1)+1)*NELEM(2,2)
      DO I=1,NELEM(2,1)+1
       DELT(I)=(NODE(NI+L+I-1,3)-NODE(L+I-1,3))/NELEM(2,2)
      END DO

      DO I=1,NELEM(2,2)
        NI=(NELEM(2,1)+1)*(NELEM(2,2)+1-I)
        DO J=1,NELEM(2,1)+1
            NODE(L,1)=NODE(L+NI,1)
            NODE(L,2)=COOR(1,2)
            NODE(L,3)=COOR(1,3)+DELT(J)*(I-1)
            L=L+1
         END DO
      END DO
!**** THE 3RD SURFACE (Y-Z,NORMAL : -X)********************************
      L=1+NEL(1)+NEL(2)

      NI=(NELEM(3,1)+1)*NELEM(3,2)
      DO I=1,NELEM(3,1)+1
       DELT(I)=(NODE(NI+L+I-1,3)-NODE(L+I-1,3))/NELEM(3,2)
      END DO

      DO I=1,NELEM(3,2)+1
        NI=(NELEM(3,1)+1)*(NELEM(3,2)+1-I)
         DO J=1,NELEM(3,1)+1
            NODE(L,1)=NODE(L+NI,1)
            NODE(L,2)=NODE(L+NI,2)
            NODE(L,3)=COOR(1,3)+DELT(J)*(I-1)
            L=L+1
         END DO
      END DO
!**** THE 4TH SURFACE (X-Z,NORMAL : X)********************************
      L=1+NEL(1)+NEL(2)+NEL(3)

      NI=(NELEM(4,1)+1)*NELEM(4,2)
      DO I=1,NELEM(4,1)+1
       DELT(I)=(NODE(NI+L+I-1,3)-NODE(L+I-1,3))/NELEM(4,2)
      END DO

      DO I=1,NELEM(4,2)
        NI=(NELEM(4,1)+1)*(NELEM(4,2)+1-I)
        DO J=1,NELEM(4,1)+1
            NODE(L,1)=NODE(L+NI,1)
            NODE(L,2)=COOR(4,2)
            NODE(L,3)=COOR(1,3)+DELT(J)*(I-1)
            L=L+1
         END DO
      END DO
!**** THE 5TH SURFACE (Y-Z,NORMAL : -X)********************************
      L=1+NEL(1)+NEL(2)+NEL(3)+NEL(4)

      NI=(NELEM(5,1)+1)*NELEM(5,2)
      DO I=1,NELEM(3,1)+1
       DELT(I)=(NODE(NI+L+I-1,3)-NODE(L+I-1,3))/NELEM(5,2)
      END DO

      DO I=1,NELEM(5,2)+1
        NI=(NELEM(5,1)+1)*(NELEM(5,2)+1-I)
         DO J=1,NELEM(5,1)+1
            NODE(L,1)=NODE(L+NI,1)
            NODE(L,2)=NODE(L+NI,2)
            NODE(L,3)=COOR(1,3)+DELT(J)*(I-1)
            L=L+1
         END DO
      END DO
!**** THE 6TH SURFACE (Y-Z,NORMAL : -Z)********************************
      K=1
      L=1+NEL(1)+NEL(2)+NEL(3)+NEL(4)+NEL(5)
      DO I=1,NELEM(6,2)+1
         DO J=1,NELEM(6,1)+1
            NODE(L,1)=NODE(K,1)
            NODE(L,2)=NODE(K,2)
            L=L+1
            K=K+1
         END DO
      END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE BOUND(NNODE,NE,NODE,PHI,PPHI,PHIT,VEL,WC)
!**********************************************************************
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,NE(6)
       REAL    R,VEL,WC
       REAL    NODE(NNODE,3),PHI(NNODE),PPHI(NNODE),PHIT(NNODE)

!   FREE SURFACE
       DO I=1,NE(1)
        PHI(I)=PHI(I)
       END DO

!   2ND SURFACE
       DO I=NE(1)+1,NE(2)
        PPHI(I)=-PHIT(I)/WC
       END DO

!   3RD SURFACE
       DO I=NE(2)+1,NE(3)
        PPHI(I)=-PHIT(I)/WC
       END DO

!   4TH SURFACE
       DO I=NE(3)+1,NE(4)
        PPHI(I)=-PHIT(I)/WC
       END DO

!   5TH SURFACE WAVE MAKER
       DO I=NE(4)+1,NE(5)
		R=SQRT(NODE(I,1)**2+NODE(I,3)**2)
		PPHI(I)=-VEL*R
       END DO

!   6TH SURFACE
       DO I=NE(5)+1,NNODE
        PPHI(I)=0.0
      END DO

       RETURN
       END 
!**********************************************************************
      SUBROUTINE KERNEL(KER1,KER2,NNODE,ENUM,NELEM,ME,LN,NODE,SHA1,SHA2,SHA3,SHA4,SH,W2,RT1,RT2,&
					   &NORMX,NORMY,NORMZ,JCB,QPI,EPI)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      INTEGER  NNODE,ENUM,NELEM(6,2),ME(6),LN(ENUM,4)
	  REAL QPI,EPI
      REAL NODE(NNODE,3),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
      REAL NORMX(ENUM,7),NORMY(ENUM,7),NORMZ(ENUM,7),JCB(ENUM,7)
      REAL SH(4,7),W2(7),RT1(7),RT2(7),SHA1(7),SHA2(7),SHA3(7),SHA4(7)
	  REAL NX,NY,NZ,THI1,THI2,LO1,LO2,C1,S1,C2,S2,SHAPE1,SHAPE2
      REAL PXI1(3),PXI2(3),H(4),G(4),X1(2),X2(2),XFUNC(7),YFUNC(7),ZFUNC(7)
      REAL XF(2),YF(2),ZF(2),RX(2),JCBS(2)

      KER1=0.0
      KER2=0.0
!****** CALCULATE THE JACOBIAN
      DO J=1,ENUM
        DO M=1,7
        DO L=1,3
          PXI1(L)=0.25*(-(1-RT2(M))*NODE(LN(J,1),L)+(1-RT2(M))*NODE(LN(J,2),L)+(1+RT2(M))*NODE(LN(J,3),L)+(-(1+RT2(M)))*NODE(LN(J,4),L))
          PXI2(L)=0.25*(-(1-RT1(M))*NODE(LN(J,1),L)+(-(1+RT1(M)))*NODE(LN(J,2),L)+(1+RT1(M))*NODE(LN(J,3),L)+(1-RT1(M))*NODE(LN(J,4),L))
		END DO
        NX=PXI1(2)*PXI2(3)-PXI1(3)*PXI2(2)
        NY=PXI1(1)*PXI2(3)-PXI1(3)*PXI2(1)
        NZ=PXI1(1)*PXI2(2)-PXI1(2)*PXI2(1)
        JCB(J,M)=(NX**2+NY**2+NZ**2)**0.5
        NORMX(J,M)=NX/JCB(J,M)
        NORMY(J,M)=NY/JCB(J,M)
        NORMZ(J,M)=NZ/JCB(J,M)
        END DO

        DO M=1,7
			IF (J .GT. ME(5)) THEN
			  NORMZ(J,M)=-NORMZ(J,M)
			ELSE IF ((J .GT. ME(1)) .AND. (J .LE. ME(2))) THEN
			  NORMY(J,M)=-NORMY(J,M)
			ELSE IF ((J .GT. ME(4)) .AND. (J .LE. ME(5))) THEN
			  NORMX(J,M)=-NORMX(J,M)
			END IF
        END DO
      END DO

!*****THE SURFACE KERNELS******************************************
      DO I = 1,NNODE
       DO J=1,ENUM
       DO M=1,7
          XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)+SHA3(M)*NODE(LN(J,3),1)+SHA4(M)*NODE(LN(J,4),1)
          YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)+SHA3(M)*NODE(LN(J,3),2)+SHA4(M)*NODE(LN(J,4),2)
          ZFUNC(M)=SHA1(M)*NODE(LN(J,1),3)+SHA2(M)*NODE(LN(J,2),3)+SHA3(M)*NODE(LN(J,3),3)+SHA4(M)*NODE(LN(J,4),3)
        END DO

        IN=0
        DO K=1,4
         RD=((NODE(I,1)-NODE(LN(J,K),1))**2+(NODE(I,2)-NODE(LN(J,K),2))**2+(NODE(I,3)-NODE(LN(J,K),3))**2)**0.5
			IF (RD .LE. 0.000001) THEN
			 IN=IN+1
			 KN=K
			END IF
        END DO

        DO K=1,4
         G(K)=0.0
         H(K)=0.0
         IF (IN .NE. 0) THEN
!      ****CALCULATE THE SINGULAR POINT****
!      ====KER1====
          IF (LN(J,K) .EQ. I) THEN
           IDELT=1
          ELSE
           IDELT=0
          END IF

          DO M=1,7
            H(K)=H(K)+((-QPI)/(((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2+&
				&(ZFUNC(M)-NODE(I,3))**2)**0.5)**3*((XFUNC(M)-NODE(I,1))*NORMX(J,M)+&
				&(YFUNC(M)-NODE(I,2))*NORMY(J,M)+(ZFUNC(M)-NODE(I,3))*NORMZ(J,M))*&
				&JCB(J,M)*(SH(K,M)-IDELT)*W2(M))
          END DO
!      ====KER2====
!      ---- WHEN ON A CORNER POINT ----
         DO M=1,7
          THI1=EPI*(RT1(M)+1)
          C1=COS(THI1)
          S1=SIN(THI1)
          LO1 =1.0/C1*(RT2(M)+1)
          THI2=EPI*(RT1(M)+3)
          C2=COS(THI2)
          S2=SIN(THI2)
          LO2 =1.0/S2*(RT2(M)+1)
         SELECT CASE(KN)
          CASE(1)
            X1(1)=LO1*C1-1
            X2(1)=LO1*S1-1
            X1(2)=LO2*C2-1
            X2(2)=LO2*S2-1
          CASE(2)
            X1(1)=1-LO1*S1
            X2(1)=LO1*C1-1
            X1(2)=1-LO2*S2
            X2(2)=LO2*C2-1
          CASE(3)
            X1(1)=1-LO1*C1
            X2(1)=1-LO1*S1
            X1(2)=1-LO2*C2
            X2(2)=1-LO2*S2
          CASE(4)
            X1(1)=LO1*S1-1
            X2(1)=1-LO1*C1
            X1(2)=LO2*S2-1
            X2(2)=1-LO2*C2
          END SELECT

      CALL SHAFUN(SHAPE1,X1(1),X2(1),K)
      CALL SHAFUN(SHAPE2,X1(2),X2(2),K)
        DO L=1,2
			DO LS=1,3
			  PXI1(LS)=0.25*(-(1-X2(L))*NODE(LN(J,1),LS)+(1-X2(L))*NODE(LN(J,2),LS)+(1+X2(L))*NODE(LN(J,3),LS)+(-(1+X2(L)))*NODE(LN(J,4),LS))
			  PXI2(LS)=0.25*(-(1-X1(L))*NODE(LN(J,1),LS)+(-(1+X1(L)))*NODE(LN(J,2),LS)+(1+X1(L))*NODE(LN(J,3),LS)+(1-X1(L))*NODE(LN(J,4),LS))
			END DO
        NX=PXI1(2)*PXI2(3)-PXI1(3)*PXI2(2)
        NY=PXI1(1)*PXI2(3)-PXI1(3)*PXI2(1)
        NZ=PXI1(1)*PXI2(2)-PXI1(2)*PXI2(1)
        JCBS(L)=(NX**2+NY**2+NZ**2)**0.5

        XF(L)=0.25*(1-X1(L))*(1-X2(L))*NODE(LN(J,1),1)+0.25*(1+X1(L))*(1-X2(L))*NODE(LN(J,2),1)+0.25*(1+X1(L))*(1+X2(L))*NODE(LN(J,3),1)+0.25*(1-X1(L))*(1+X2(L))*NODE(LN(J,4),1)
        YF(L)=0.25*(1-X1(L))*(1-X2(L))*NODE(LN(J,1),2)+0.25*(1+X1(L))*(1-X2(L))*NODE(LN(J,2),2)+0.25*(1+X1(L))*(1+X2(L))*NODE(LN(J,3),2)+0.25*(1-X1(L))*(1+X2(L))*NODE(LN(J,4),2)
        ZF(L)=0.25*(1-X1(L))*(1-X2(L))*NODE(LN(J,1),3)+0.25*(1+X1(L))*(1-X2(L))*NODE(LN(J,2),3)+0.25*(1+X1(L))*(1+X2(L))*NODE(LN(J,3),3)+0.25*(1-X1(L))*(1+X2(L))*NODE(LN(J,4),3)
        RX(L)=((XF(L)-NODE(I,1))**2+(YF(L)-NODE(I,2))**2+(ZF(L)-NODE(I,3))**2)**0.5
        END DO

       G(K)=G(K)+(1./(RX(1)*C1)*LO1*SHAPE1*JCBS(1)*W2(M))/32+(1./(RX(2)*S2)*LO2*SHAPE2*JCBS(2)*W2(M))/32
      END DO

      ELSE
!     ==== NON SINGULER TERM ====
          DO M=1,7
            H(K)=H(K)+(-QPI/(((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2+&
					&(ZFUNC(M)-NODE(I,3))**2)**0.5)**3*((XFUNC(M)-NODE(I,1))*NORMX(J,M)+&
					&(YFUNC(M)-NODE(I,2))*NORMY(J,M)+(ZFUNC(M)-NODE(I,3))*NORMZ(J,M))*&
					&JCB(J,M)*SH(K,M)*W2(M))
            G(K)=G(K)+QPI/(((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2+&
					&(ZFUNC(M)-NODE(I,3))**2)**0.5)*JCB(J,M)*SH(K,M)*W2(M)
          END DO
      END IF

         KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
      END DO
      END DO
      END DO
!********************************************************************
       DO I=1,NNODE
          SIGMA=0.0
          DO J=1,NNODE
             SIGMA=KER1(I,J)+SIGMA
          END DO
          KER1(I,I)=-SIGMA+KER1(I,I)
       END DO

       RETURN
       END
!**********************************************************************
       SUBROUTINE SOLVE(NNODE,NEL,BUNTYP,KER1,KER2,PHI,PPHI,SOR)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI USING GAUSSIAN SEIDEL INTERATIVE
       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL    (A-H,O-Z)
       INTEGER NNODE,NEL(6),BUNTYP(6)
       REAL    SOR
       REAL    KER1(NNODE,NNODE),KER2(NNODE,NNODE),PHI(NNODE),PPHI(NNODE)
       REAL    H1(NNODE,NNODE),G1(NNODE,NNODE),Q1(NNODE),P1(NNODE),TEMP(NNODE)
       REAL    G2(NNODE,NNODE+1),X1(NNODE),XP(NNODE),ERR(NNODE)
       REAL*8  EMAX 

!********** TO MOVE THE KER1 AND KER2**********************************
!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,6
          DO J=K+N,NEL(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
               P1(J)=PPHI(J)
             ELSE
               Q1(J)=PPHI(J)
               P1(J)=PHI(J)
            END IF
          END DO
        N=N+NEL(I)
       END DO
!------KER1 KER2-------
       DO I=1,NNODE
          N=0
          DO L=1,6
            DO J=K+N,NEL(L)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
            END DO
            N=NEL(L)+N
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

!*************GSI SOR INTERATIVE*************
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
  100 CONTINUE
!*************GSI SOR INTERATIVE*************

       K=1
       N=0
       DO I=1,6
          DO J=K+N,NEL(I)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
              ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=NEL(I)+N
       END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE TAYES1(NNODE,ENUM,LN,ME,NE,CBC,NODE,PHI,PPHI,PHIT,VELX,VELY,VELZ,&
	  &NORMX,NORMY,NORMZ,RM1,RM2,THI,S1X,S1Y,S1Z,S2X,S2Y,S2Z,DPDS1,DPDS2,DT,GRAV,DEP,XI1,XI2)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL(A-H,O-Z)
      INTEGER NNODE,ENUM,LN(ENUM,4),ME(6),NE(6),CBC(NE(1),2)
      INTEGER IC(NNODE)
      REAL   NODE(NNODE,3),PHI(NNODE),PPHI(NNODE),PHIT(NNODE)
      REAL   VELX(NNODE),VELY(NNODE),VELZ(NNODE),DT(NNODE)
      REAL   DPDX(2),PXI1(3),PXI2(3),XI1(4),XI2(4)
      REAL   GRAV,DEP,NX,NY,NZ,JCB,DPX,DPY,DPZ
      REAL   NORMX(ME(1),4),NORMY(ME(1),4),NORMZ(ME(1),4)
	  REAL	 RM1(ME(1),4),RM2(ME(1),4),THI(ME(1),4)
      REAL   S1X(ME(1),4),S1Y(ME(1),4),S1Z(ME(1),4)
      REAL   S2X(ME(1),4),S2Y(ME(1),4),S2Z(ME(1),4)
      REAL   DPDS1(ME(1),4),DPDS2(ME(1),4)

       IC=0
       VELX=0.0
       VELY=0.0
       VELZ=0.0
!*********** CALCULATE PARTICLE VELOCITY ON FREE SURFACE ***********

      DO I=1,ME(1)
      DO J=1,4

!----- CALCULATE THE NORMAL VECTOR
       DO L=1,3
         PXI1(L)=0.25*(-(1-XI2(J))*NODE(LN(I,1),L)+(1-XI2(J))*NODE(LN(I,2),L)+(1+XI2(J))*NODE(LN(I,3),L)+(-(1+XI2(J)))*NODE(LN(I,4),L))
         PXI2(L)=0.25*(-(1-XI1(J))*NODE(LN(I,1),L)+(-(1+XI1(J)))*NODE(LN(I,2),L)+(1+XI1(J))*NODE(LN(I,3),L)+(1-XI1(J))*NODE(LN(I,4),L))
       END DO

        NX=PXI1(2)*PXI2(3)-PXI1(3)*PXI2(2)
        NY=PXI1(1)*PXI2(3)-PXI1(3)*PXI2(1)
        NZ=PXI1(1)*PXI2(2)-PXI1(2)*PXI2(1)
        JCB=(NX**2+NY**2+NZ**2)**0.5
        NORMX(I,J)=NX/JCB
        NORMY(I,J)=NY/JCB
        NORMZ(I,J)=NZ/JCB

!----- CALCULATE THE TANGENTIAL 1,2 VECTOR
        RM1(I,J)=(PXI1(1)**2+PXI1(2)**2+PXI1(3)**2)**0.5
        RM2(I,J)=(PXI2(1)**2+PXI2(2)**2+PXI2(3)**2)**0.5
        THI(I,J)=ACOS(1./(RM1(I,J)*RM2(I,J))*(PXI1(1)*PXI2(1)+PXI1(2)*PXI2(2)+PXI1(3)*PXI2(3)))

        S1X(I,J)=PXI1(1)/RM1(I,J)
        S1Y(I,J)=PXI1(2)/RM1(I,J)
        S1Z(I,J)=PXI1(3)/RM1(I,J)

        S2X(I,J)=NORMY(I,J)*S1Z(I,J)-NORMZ(I,J)*S1Y(I,J)
        S2Y(I,J)=NORMZ(I,J)*S1X(I,J)-NORMX(I,J)*S1Z(I,J)
        S2Z(I,J)=NORMX(I,J)*S1Y(I,J)-NORMY(I,J)*S1X(I,J)

!----- CALCULATE THE TANGENTIAL DIFFERENTIAL OF PHI
         DPDX(1)=0.25*(-(1-XI2(J))*PHI(LN(I,1))+(1-XI2(J))*PHI(LN(I,2))+(1+XI2(J))*PHI(LN(I,3))+(-(1+XI2(J)))*PHI(LN(I,4)))
         DPDX(2)=0.25*(-(1-XI1(J))*PHI(LN(I,1))+(-(1+XI1(J)))*PHI(LN(I,2))+(1+XI1(J))*PHI(LN(I,3))+(1-XI1(J))*PHI(LN(I,4)))

         DPDS1(I,J)=DPDX(1)/RM1(I,J)
         DPDS2(I,J)=-DPDX(1)*COS(THI(I,J))/RM1(I,J)/SIN(THI(I,J))+DPDX(2)/RM2(I,J)/SIN(THI(I,J))

!        RX1=ABS(NODE(LN(I,J),1)-NODE(1,1))		! DISTANT TO X = 0
!        RX2=ABS(NODE(LN(I,J),1)-NODE(NE(1),1))	! DISTANT TO X = L
!        RY1=ABS(NODE(LN(I,J),2)-NODE(1,2))		! DISTANT TO Y = 0
!        RY2=ABS(NODE(LN(I,J),2)-NODE(NE(1),2))	! DISTANT TO Y = B
!----- X DIRCTION----------
        DPX=DPDS1(I,J)*S1X(I,J)+DPDS2(I,J)*S2X(I,J)+PPHI(LN(I,J))*NORMX(I,J)

        IF (CBC(LN(I,J),1)==1) THEN
			DPX=-PPHI(CBC(LN(I,J),2))
			DPDS1(I,J)=(DPX-PPHI(LN(I,J))*NORMX(I,J)-DPDS2(I,J)*S2X(I,J))/S1X(I,J)
        ELSE IF (CBC(LN(I,J),1)==3) THEN
			DPX=PPHI(CBC(LN(I,J),2))
			DPDS1(I,J)=(DPX-PPHI(LN(I,J))*NORMX(I,J)-DPDS2(I,J)*S2X(I,J))/S1X(I,J)
		END IF

!IF (RX1 .LE. 0.00001) THEN
!DPX=-PPHI(NE(5))
!DPDS1(I,J)=(DPX-PPHI(LN(I,J))*NORMX(I,J)-DPDS2(I,J)*S2X(I,J))/S1X(I,J)
!ELSE IF (RX2 .LE. 0.00001) THEN
!DPX=PPHI(NE(3))
!DPDS1(I,J)=(DPX-PPHI(LN(I,J))*NORMX(I,J)-DPDS2(I,J)*S2X(I,J))/S1X(I,J)
!END IF

!--- Y DIRECTION----------
        DPY=DPDS1(I,J)*S1Y(I,J)+DPDS2(I,J)*S2Y(I,J)+PPHI(LN(I,J))*NORMY(I,J)

        IF (CBC(LN(I,J),1)==2) THEN
          DPY=-PPHI(CBC(LN(I,J),2))
          DPDS2(I,J)=(DPY-PPHI(LN(I,J))*NORMY(I,J)-DPDS1(I,J)*S1Y(I,J))/S2Y(I,J)
		ELSE IF (CBC(LN(I,J),1)==4) THEN
          DPY=PPHI(CBC(LN(I,J),2))
          DPDS2(I,J)=(DPY-PPHI(LN(I,J))*NORMY(I,J)-DPDS1(I,J)*S1Y(I,J))/S2Y(I,J)
		END IF

!IF (RY1 .LE. 0.00001) THEN
!DPY=0.0
!DPDS2(I,J)=(DPY-PPHI(LN(I,J))*NORMY(I,J)-DPDS1(I,J)*S1Y(I,J))/S2Y(I,J)
!ELSE IF (RY2 .LE. 0.00001) THEN
!DPY=0.0
!DPDS2(I,J)=(DPY-PPHI(LN(I,J))*NORMY(I,J)-DPDS1(I,J)*S1Y(I,J))/S2Y(I,J)
!END IF

!-- Z DIRECTION --------
        DPZ=DPDS1(I,J)*S1Z(I,J)+DPDS2(I,J)*S2Z(I,J)+PPHI(LN(I,J))*NORMZ(I,J)

      IC(LN(I,J))=IC(LN(I,J))+1
      VELX(LN(I,J))=VELX(LN(I,J))+DPX
      VELY(LN(I,J))=VELY(LN(I,J))+DPY
      VELZ(LN(I,J))=VELZ(LN(I,J))+DPZ
      END DO
      END DO

!----- FREE SURFACE VEL XYZ AND PHIT AND DT -----------
      DO I=1,NE(1)
       VELX(I)=VELX(I)/IC(I)
       VELY(I)=VELY(I)/IC(I)
       VELZ(I)=VELZ(I)/IC(I)
       DT(I)=0.5*(VELX(I)**2+VELY(I)**2+VELZ(I)**2)-GRAV*(NODE(I,3)-DEP)
       PHIT(I)=DT(I)-(VELX(I)**2+VELY(I)**2+VELZ(I)**2)
      END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE ACCBC(NNODE,ENUM,NELEM,LN,NEL,ME,ACBC,NODE,PHI,PPHI,VELX,VELY,VELZ,ACCMO,PPPHI,VEL)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL    (A-H,O-Z)
      INTEGER  NNODE,ENUM,NELEM(6,2),LN(ENUM,4),NEL(6),ME(6),ACBC(NNODE,2)
      REAL NODE(NNODE,3),PHI(NNODE),PPHI(NNODE),PPPHI(NNODE)
      REAL VELX(NNODE),VELY(NNODE),VELZ(NNODE),ACCMO(NNODE)
      REAL NX,NY,NZ,JCB,RM1,RM2,THI,XMOD,YMOD,ZMOD
	  REAL DPDX(2),DPDS(2),S1(3),S2(3),PX1(3),PX2(3),NORM(3),XI1(4),XI2(4)

!*********** 2ND PLANE ***********
      DO I=ME(1)+1,ME(2)
       DO J=1,4
!***** CALCULATE THE JACOBIAN OF XI1 AND XI2 
       DO L=1,3
          PX1(L)=0.25*(-(1-XI2(J))*NODE(LN(I,1),L)+(1-XI2(J))*NODE(LN(I,2),L)+(1+XI2(J))*NODE(LN(I,3),L)+(-(1+XI2(J)))*NODE(LN(I,4),L))
          PX2(L)=0.25*(-(1-XI1(J))*NODE(LN(I,1),L)+(-(1+XI1(J)))*NODE(LN(I,2),L)+(1+XI1(J))*NODE(LN(I,3),L)+(1-XI1(J))*NODE(LN(I,4),L))
        END DO
        NX=PX1(2)*PX2(3)-PX1(3)*PX2(2)
        NY=PX1(1)*PX2(3)-PX1(3)*PX2(1)
        NZ=PX1(1)*PX2(2)-PX1(2)*PX2(1)
        JCB=(NX**2+NY**2+NZ**2)**0.5
        NORM(1)=NX/JCB
        NORM(2)=NY/JCB
        NORM(3)=NZ/JCB

!***** CALCULATE THE  DIFFERENTIAL OF PHI ,XI1,XI2
         RM1=(PX1(1)**2+PX1(2)**2+PX1(3)**2)**0.5
         RM2=(PX2(1)**2+PX2(2)**2+PX2(3)**2)**0.5
         THI=ACOS(1./(RM1*RM2)*(PX1(1)*PX2(1)+PX1(2)*PX2(2)+PX1(3)*PX2(3)))

         DPDX(1)=0.25*(-(1-XI2(J))*PHI(LN(I,1))+(1-XI2(J))*PHI(LN(I,2))+(1+XI2(J))*PHI(LN(I,3))+(-(1+XI2(J)))*PHI(LN(I,4)))
         DPDX(2)=0.25*(-(1-XI1(J))*PHI(LN(I,1))+(-(1+XI1(J)))*PHI(LN(I,2))+(1+XI1(J))*PHI(LN(I,3))+(1-XI1(J))*PHI(LN(I,4)))

         DPDS(1)=DPDX(1)/RM1
         DPDS(2)=-DPDX(1)*COS(THI)/RM1/SIN(THI)+DPDX(2)/RM2/SIN(THI)

		 S1(1)=PX1(1)/RM1
         S1(2)=PX1(2)/RM1
         S1(3)=PX1(3)/RM1

         S2(1)=NORM(2)*S1(3)-NORM(3)*S1(2)
         S2(2)=NORM(3)*S1(1)-NORM(1)*S1(3)
         S2(3)=NORM(1)*S1(2)-NORM(2)*S1(1)

!------THE FOLLOW IS FOR THE BOUNDARYT CONDITION
         DPPDX1=0.25*(-(1-XI2(J))*PPHI(LN(I,1))+(1-XI2(J))*PPHI(LN(I,2))+(1+XI2(J))*PPHI(LN(I,3))+(-(1+XI2(J)))*PPHI(LN(I,4)))
         DPPDX2=0.25*(-(1-XI1(J))*PPHI(LN(I,1))+(-(1+XI1(J)))*PPHI(LN(I,2))+(1+XI1(J))*PPHI(LN(I,3))+(1-XI1(J))*PPHI(LN(I,4)))

         DPPDS1=DPPDX1/RM1
         DPPDS2=-DPPDX1*COS(THI)/RM1/SIN(THI)+DPPDX2/RM2/SIN(THI)

         DPDS1X1=0.25*(-(1-XI2(J))*DPDS(1)+(1-XI2(J))*DPDS(1)+(1+XI2(J))*DPDS(1)+(-(1+XI2(J)))*DPDS(1))
         DPDS2X1=0.25*(-(1-XI2(J))*DPDS(2)+(1-XI2(J))*DPDS(2)+(1+XI2(J))*DPDS(2)+(-(1+XI2(J)))*DPDS(2))
         DPDS2X2=0.25*(-(1-XI1(J))*DPDS(2)+(-(1+XI1(J)))*DPDS(2)+(1+XI1(J))*DPDS(2)+(1-XI1(J))*DPDS(2))

         DPDS1S1=DPDS1X1/RM1
         DPDS2S2=-DPDS2X1*COS(THI)/RM1/SIN(THI)-DPDS2X2/RM2/SIN(THI)

         DPDNN=-(DPDS1S1+DPDS2S2)
		 PPPHI(LN(I,J))=DPDNN

         XMOD=DPDNN*NORM(1)+DPPDS1*S1(1)+DPPDS2*S2(1)
         YMOD=DPDNN*NORM(2)+DPPDS1*S1(2)+DPPDS2*S2(2)
         ZMOD=DPDNN*NORM(3)+DPPDS1*S1(3)+DPPDS2*S2(3)

!***** CALCULATE THE VELOCITY OF LATERAL PARTICLE AND PHI OVER T
       VELX(LN(I,J))=NORM(1)*PPHI(LN(I,J))+S1(1)*DPDS(1)+S2(1)*DPDS(2)
       VELY(LN(I,J))=NORM(2)*PPHI(LN(I,J))+S1(2)*DPDS(1)+S2(2)*DPDS(2)
       VELZ(LN(I,J))=NORM(3)*PPHI(LN(I,J))+S1(3)*DPDS(1)+S2(3)*DPDS(2)

		   IF(NODE(LN(I,J),3) .EQ. 0.) THEN
			VELZ(LN(I,J))=0.0
		   END IF

		   IF (ACBC(LN(I,J),1)==2) THEN
			VELZ(LN(I,J))=VELZ(ACBC(LN(I,J),2))
		   END IF

       ACCMO(LN(I,J))=VELX(LN(I,J))*XMOD+VELY(LN(I,J))*YMOD+VELZ(LN(I,J))*ZMOD
      END DO
      END DO

!*********** 3RD PLANE ***********
      DO I=ME(2)+1,ME(3)
       DO J=1,4
!***** CALCULATE THE JACOBIAN OF XI1 AND XI2 
       DO L=1,3
          PX1(L)=0.25*(-(1-XI2(J))*NODE(LN(I,1),L)+(1-XI2(J))*NODE(LN(I,2),L)+(1+XI2(J))*NODE(LN(I,3),L)+(-(1+XI2(J)))*NODE(LN(I,4),L))
          PX2(L)=0.25*(-(1-XI1(J))*NODE(LN(I,1),L)+(-(1+XI1(J)))*NODE(LN(I,2),L)+(1+XI1(J))*NODE(LN(I,3),L)+(1-XI1(J))*NODE(LN(I,4),L))
        END DO
        NX=PX1(2)*PX2(3)-PX1(3)*PX2(2)
        NY=PX1(1)*PX2(3)-PX1(3)*PX2(1)
        NZ=PX1(1)*PX2(2)-PX1(2)*PX2(1)
        JCB=(NX**2+NY**2+NZ**2)**0.5
        NORM(1)=NX/JCB
        NORM(2)=NY/JCB
        NORM(3)=NZ/JCB

!***** CALCULATE THE  DIFFERENTIAL OF PHI ,XI1,XI2
         RM1=(PX1(1)**2+PX1(2)**2+PX1(3)**2)**0.5
         RM2=(PX2(1)**2+PX2(2)**2+PX2(3)**2)**0.5
         THI=ACOS(1./(RM1*RM2)*(PX1(1)*PX2(1)+PX1(2)*PX2(2)+PX1(3)*PX2(3)))

         DPDX(1)=0.25*(-(1-XI2(J))*PHI(LN(I,1))+(1-XI2(J))*PHI(LN(I,2))+(1+XI2(J))*PHI(LN(I,3))+(-(1+XI2(J)))*PHI(LN(I,4)))
         DPDX(2)=0.25*(-(1-XI1(J))*PHI(LN(I,1))+(-(1+XI1(J)))*PHI(LN(I,2))+(1+XI1(J))*PHI(LN(I,3))+(1-XI1(J))*PHI(LN(I,4)))

         DPDS(1)=DPDX(1)/RM1
         DPDS(2)=-DPDX(1)*COS(THI)/RM1/SIN(THI)+DPDX(2)/RM2/SIN(THI)

		 S1(1)=PX1(1)/RM1
         S1(2)=PX1(2)/RM1
         S1(3)=PX1(3)/RM1

         S2(1)=NORM(2)*S1(3)-NORM(3)*S1(2)
         S2(2)=NORM(3)*S1(1)-NORM(1)*S1(3)
         S2(3)=NORM(1)*S1(2)-NORM(2)*S1(1)

!------THE FOLLOW IS FOR THE BOUNDARYT CONDITION
         DPPDX1=0.25*(-(1-XI2(J))*PPHI(LN(I,1))+(1-XI2(J))*PPHI(LN(I,2))+(1+XI2(J))*PPHI(LN(I,3))+(-(1+XI2(J)))*PPHI(LN(I,4)))
         DPPDX2=0.25*(-(1-XI1(J))*PPHI(LN(I,1))+(-(1+XI1(J)))*PPHI(LN(I,2))+(1+XI1(J))*PPHI(LN(I,3))+(1-XI1(J))*PPHI(LN(I,4)))

         DPPDS1=DPPDX1/RM1
         DPPDS2=-DPPDX1*COS(THI)/RM1/SIN(THI)+DPPDX2/RM2/SIN(THI)

         DPDS1X1=0.25*(-(1-XI2(J))*DPDS(1)+(1-XI2(J))*DPDS(1)+(1+XI2(J))*DPDS(1)+(-(1+XI2(J)))*DPDS(1))
         DPDS2X1=0.25*(-(1-XI2(J))*DPDS(2)+(1-XI2(J))*DPDS(2)+(1+XI2(J))*DPDS(2)+(-(1+XI2(J)))*DPDS(2))
         DPDS2X2=0.25*(-(1-XI1(J))*DPDS(2)+(-(1+XI1(J)))*DPDS(2)+(1+XI1(J))*DPDS(2)+(1-XI1(J))*DPDS(2))

         DPDS1S1=DPDS1X1/RM1
         DPDS2S2=-DPDS2X1*COS(THI)/RM1/SIN(THI)-DPDS2X2/RM2/SIN(THI)

         DPDNN=-(DPDS1S1+DPDS2S2)
		 PPPHI(LN(I,J))=DPDNN

         XMOD=DPDNN*NORM(1)+DPPDS1*S1(1)+DPPDS2*S2(1)
         YMOD=DPDNN*NORM(2)+DPPDS1*S1(2)+DPPDS2*S2(2)
         ZMOD=DPDNN*NORM(3)+DPPDS1*S1(3)+DPPDS2*S2(3)

!***** CALCULATE THE VELOCITY OF LATERAL PARTICLE AND PHI OVER T
       VELX(LN(I,J))=NORM(1)*PPHI(LN(I,J))+S1(1)*DPDS(1)+S2(1)*DPDS(2)
       VELY(LN(I,J))=NORM(2)*PPHI(LN(I,J))+S1(2)*DPDS(1)+S2(2)*DPDS(2)
       VELZ(LN(I,J))=NORM(3)*PPHI(LN(I,J))+S1(3)*DPDS(1)+S2(3)*DPDS(2)

		   IF(NODE(LN(I,J),3) .EQ. 0.) THEN
			VELZ(LN(I,J))=0.0
		   END IF

		   IF (ACBC(LN(I,J),1)==3) THEN
			VELZ(LN(I,J))=VELZ(ACBC(LN(I,J),2))
		   END IF

       ACCMO(LN(I,J))=VELX(LN(I,J))*XMOD+VELY(LN(I,J))*YMOD+VELZ(LN(I,J))*ZMOD
      END DO
      END DO

!*********** 4TH PLANE ***********
      DO I=ME(3)+1,ME(4)
       DO J=1,4
!***** CALCULATE THE JACOBIAN OF XI1 AND XI2 
       DO L=1,3
          PX1(L)=0.25*(-(1-XI2(J))*NODE(LN(I,1),L)+(1-XI2(J))*NODE(LN(I,2),L)+(1+XI2(J))*NODE(LN(I,3),L)+(-(1+XI2(J)))*NODE(LN(I,4),L))
          PX2(L)=0.25*(-(1-XI1(J))*NODE(LN(I,1),L)+(-(1+XI1(J)))*NODE(LN(I,2),L)+(1+XI1(J))*NODE(LN(I,3),L)+(1-XI1(J))*NODE(LN(I,4),L))
        END DO
        NX=PX1(2)*PX2(3)-PX1(3)*PX2(2)
        NY=PX1(1)*PX2(3)-PX1(3)*PX2(1)
        NZ=PX1(1)*PX2(2)-PX1(2)*PX2(1)
        JCB=(NX**2+NY**2+NZ**2)**0.5
        NORM(1)=NX/JCB
        NORM(2)=NY/JCB
        NORM(3)=NZ/JCB

!***** CALCULATE THE  DIFFERENTIAL OF PHI ,XI1,XI2
         RM1=(PX1(1)**2+PX1(2)**2+PX1(3)**2)**0.5
         RM2=(PX2(1)**2+PX2(2)**2+PX2(3)**2)**0.5
         THI=ACOS(1./(RM1*RM2)*(PX1(1)*PX2(1)+PX1(2)*PX2(2)+PX1(3)*PX2(3)))

         DPDX(1)=0.25*(-(1-XI2(J))*PHI(LN(I,1))+(1-XI2(J))*PHI(LN(I,2))+(1+XI2(J))*PHI(LN(I,3))+(-(1+XI2(J)))*PHI(LN(I,4)))
         DPDX(2)=0.25*(-(1-XI1(J))*PHI(LN(I,1))+(-(1+XI1(J)))*PHI(LN(I,2))+(1+XI1(J))*PHI(LN(I,3))+(1-XI1(J))*PHI(LN(I,4)))

         DPDS(1)=DPDX(1)/RM1
         DPDS(2)=-DPDX(1)*COS(THI)/RM1/SIN(THI)+DPDX(2)/RM2/SIN(THI)

		 S1(1)=PX1(1)/RM1
         S1(2)=PX1(2)/RM1
         S1(3)=PX1(3)/RM1

         S2(1)=NORM(2)*S1(3)-NORM(3)*S1(2)
         S2(2)=NORM(3)*S1(1)-NORM(1)*S1(3)
         S2(3)=NORM(1)*S1(2)-NORM(2)*S1(1)

!------THE FOLLOW IS FOR THE BOUNDARYT CONDITION
         DPPDX1=0.25*(-(1-XI2(J))*PPHI(LN(I,1))+(1-XI2(J))*PPHI(LN(I,2))+(1+XI2(J))*PPHI(LN(I,3))+(-(1+XI2(J)))*PPHI(LN(I,4)))
         DPPDX2=0.25*(-(1-XI1(J))*PPHI(LN(I,1))+(-(1+XI1(J)))*PPHI(LN(I,2))+(1+XI1(J))*PPHI(LN(I,3))+(1-XI1(J))*PPHI(LN(I,4)))

         DPPDS1=DPPDX1/RM1
         DPPDS2=-DPPDX1*COS(THI)/RM1/SIN(THI)+DPPDX2/RM2/SIN(THI)

         DPDS1X1=0.25*(-(1-XI2(J))*DPDS(1)+(1-XI2(J))*DPDS(1)+(1+XI2(J))*DPDS(1)+(-(1+XI2(J)))*DPDS(1))
         DPDS2X1=0.25*(-(1-XI2(J))*DPDS(2)+(1-XI2(J))*DPDS(2)+(1+XI2(J))*DPDS(2)+(-(1+XI2(J)))*DPDS(2))
         DPDS2X2=0.25*(-(1-XI1(J))*DPDS(2)+(-(1+XI1(J)))*DPDS(2)+(1+XI1(J))*DPDS(2)+(1-XI1(J))*DPDS(2))

         DPDS1S1=DPDS1X1/RM1
         DPDS2S2=-DPDS2X1*COS(THI)/RM1/SIN(THI)-DPDS2X2/RM2/SIN(THI)

         DPDNN=-(DPDS1S1+DPDS2S2)
		 PPPHI(LN(I,J))=DPDNN

         XMOD=DPDNN*NORM(1)+DPPDS1*S1(1)+DPPDS2*S2(1)
         YMOD=DPDNN*NORM(2)+DPPDS1*S1(2)+DPPDS2*S2(2)
         ZMOD=DPDNN*NORM(3)+DPPDS1*S1(3)+DPPDS2*S2(3)

!***** CALCULATE THE VELOCITY OF LATERAL PARTICLE AND PHI OVER T
       VELX(LN(I,J))=NORM(1)*PPHI(LN(I,J))+S1(1)*DPDS(1)+S2(1)*DPDS(2)
       VELY(LN(I,J))=NORM(2)*PPHI(LN(I,J))+S1(2)*DPDS(1)+S2(2)*DPDS(2)
       VELZ(LN(I,J))=NORM(3)*PPHI(LN(I,J))+S1(3)*DPDS(1)+S2(3)*DPDS(2)

		   IF(NODE(LN(I,J),3) .EQ. 0.) THEN
			VELZ(LN(I,J))=0.0
		   END IF

		   IF (ACBC(LN(I,J),1)==4) THEN
			VELZ(LN(I,J))=VELZ(ACBC(LN(I,J),2))
		   END IF

       ACCMO(LN(I,J))=VELX(LN(I,J))*XMOD+VELY(LN(I,J))*YMOD+VELZ(LN(I,J))*ZMOD
      END DO
      END DO

!**** 5TH PLANE WAVE MAKER
      DO I=ME(4)+1,ME(5)
       DO J=1,4
!***** CALCULATE THE JACOBIAN OF XI1 AND XI2 
       DO L=1,3
          PX1(L)=0.25*(-(1-XI2(J))*NODE(LN(I,1),L)+(1-XI2(J))*NODE(LN(I,2),L)+(1+XI2(J))*NODE(LN(I,3),L)+(-(1+XI2(J)))*NODE(LN(I,4),L))
          PX2(L)=0.25*(-(1-XI1(J))*NODE(LN(I,1),L)+(-(1+XI1(J)))*NODE(LN(I,2),L)+(1+XI1(J))*NODE(LN(I,3),L)+(1-XI1(J))*NODE(LN(I,4),L))
        END DO
        NX=PX1(2)*PX2(3)-PX1(3)*PX2(2)
        NY=PX1(1)*PX2(3)-PX1(3)*PX2(1)
        NZ=PX1(1)*PX2(2)-PX1(2)*PX2(1)
        JCB=(NX**2+NY**2+NZ**2)**0.5
        NORM(1)=NX/JCB
        NORM(2)=NY/JCB
        NORM(3)=NZ/JCB

!***** CALCULATE THE  DIFFERENTIAL OF PHI ,XI1,XI2
         RM1=(PX1(1)**2+PX1(2)**2+PX1(3)**2)**0.5
         RM2=(PX2(1)**2+PX2(2)**2+PX2(3)**2)**0.5
         THI=ACOS(1./(RM1*RM2)*(PX1(1)*PX2(1)+PX1(2)*PX2(2)+PX1(3)*PX2(3)))

         DPDX(1)=0.25*(-(1-XI2(J))*PHI(LN(I,1))+(1-XI2(J))*PHI(LN(I,2))+(1+XI2(J))*PHI(LN(I,3))+(-(1+XI2(J)))*PHI(LN(I,4)))
         DPDX(2)=0.25*(-(1-XI1(J))*PHI(LN(I,1))+(-(1+XI1(J)))*PHI(LN(I,2))+(1+XI1(J))*PHI(LN(I,3))+(1-XI1(J))*PHI(LN(I,4)))

         DPDS(1)=DPDX(1)/RM1
         DPDS(2)=-DPDX(1)*COS(THI)/RM1/SIN(THI)+DPDX(2)/RM2/SIN(THI)

		 S1(1)=PX1(1)/RM1
         S1(2)=PX1(2)/RM1
         S1(3)=PX1(3)/RM1

         S2(1)=NORM(2)*S1(3)-NORM(3)*S1(2)
         S2(2)=NORM(3)*S1(1)-NORM(1)*S1(3)
         S2(3)=NORM(1)*S1(2)-NORM(2)*S1(1)

!------THE FOLLOW IS FOR THE BOUNDARYT CONDITION
         DPPDX1=0.25*(-(1-XI2(J))*PPHI(LN(I,1))+(1-XI2(J))*PPHI(LN(I,2))+(1+XI2(J))*PPHI(LN(I,3))+(-(1+XI2(J)))*PPHI(LN(I,4)))
         DPPDX2=0.25*(-(1-XI1(J))*PPHI(LN(I,1))+(-(1+XI1(J)))*PPHI(LN(I,2))+(1+XI1(J))*PPHI(LN(I,3))+(1-XI1(J))*PPHI(LN(I,4)))

         DPPDS1=DPPDX1/RM1
         DPPDS2=-DPPDX1*COS(THI)/RM1/SIN(THI)+DPPDX2/RM2/SIN(THI)

         DPDS1X1=0.25*(-(1-XI2(J))*DPDS(1)+(1-XI2(J))*DPDS(1)+(1+XI2(J))*DPDS(1)+(-(1+XI2(J)))*DPDS(1))
         DPDS2X1=0.25*(-(1-XI2(J))*DPDS(2)+(1-XI2(J))*DPDS(2)+(1+XI2(J))*DPDS(2)+(-(1+XI2(J)))*DPDS(2))
         DPDS2X2=0.25*(-(1-XI1(J))*DPDS(2)+(-(1+XI1(J)))*DPDS(2)+(1+XI1(J))*DPDS(2)+(1-XI1(J))*DPDS(2))

         DPDS1S1=DPDS1X1/RM1
         DPDS2S2=-DPDS2X1*COS(THI)/RM1/SIN(THI)-DPDS2X2/RM2/SIN(THI)

         DPDNN=-(DPDS1S1+DPDS2S2)

         XMOD=DPDNN*NORM(1)+DPPDS1*S1(1)+DPPDS2*S2(1)
         YMOD=DPDNN*NORM(2)+DPPDS1*S1(2)+DPPDS2*S2(2)
         ZMOD=DPDNN*NORM(3)+DPPDS1*S1(3)+DPPDS2*S2(3)

!***** CALCULATE THE VELOCITY OF LATERAL PARTICLE AND PHI OVER T
       VELX(LN(I,J))=NORM(1)*PPHI(LN(I,J))+S1(1)*DPDS(1)+S2(1)*DPDS(2) !VEL
       VELY(LN(I,J))=NORM(2)*PPHI(LN(I,J))+S1(2)*DPDS(1)+S2(2)*DPDS(2) !0.0
       VELZ(LN(I,J))=NORM(3)*PPHI(LN(I,J))+S1(3)*DPDS(1)+S2(3)*DPDS(2) !S1(3)*DPDS(1)+S2(3)*DPDS(2)

		   IF(NODE(LN(I,J),3) .EQ. 0.) THEN
			VELZ(LN(I,J))=0.0
		   END IF

		   IF (ACBC(LN(I,J),1)==1) THEN
			VELZ(LN(I,J))=VELZ(ACBC(LN(I,J),2))
		   END IF

       ACCMO(LN(I,J))=VELX(LN(I,J))*XMOD+VELY(LN(I,J))*YMOD+VELZ(LN(I,J))*ZMOD
      END DO
      END DO

!*********** 6TH PLANE ***********
      DO I=ME(5)+1,ME(6)
       DO J=1,4
!***** CALCULATE THE JACOBIAN OF XI1 AND XI2 
       DO L=1,3
          PX1(L)=0.25*(-(1-XI2(J))*NODE(LN(I,1),L)+(1-XI2(J))*NODE(LN(I,2),L)+(1+XI2(J))*NODE(LN(I,3),L)+(-(1+XI2(J)))*NODE(LN(I,4),L))
          PX2(L)=0.25*(-(1-XI1(J))*NODE(LN(I,1),L)+(-(1+XI1(J)))*NODE(LN(I,2),L)+(1+XI1(J))*NODE(LN(I,3),L)+(1-XI1(J))*NODE(LN(I,4),L))
        END DO
        NX=PX1(2)*PX2(3)-PX1(3)*PX2(2)
        NY=PX1(1)*PX2(3)-PX1(3)*PX2(1)
        NZ=PX1(1)*PX2(2)-PX1(2)*PX2(1)
        JCB=(NX**2+NY**2+NZ**2)**0.5
        NORM(1)=NX/JCB
        NORM(2)=NY/JCB
        NORM(3)=NZ/JCB

!***** CALCULATE THE  DIFFERENTIAL OF PHI ,XI1,XI2
         RM1=(PX1(1)**2+PX1(2)**2+PX1(3)**2)**0.5
         RM2=(PX2(1)**2+PX2(2)**2+PX2(3)**2)**0.5
         THI=ACOS(1./(RM1*RM2)*(PX1(1)*PX2(1)+PX1(2)*PX2(2)+PX1(3)*PX2(3)))

         DPDX(1)=0.25*(-(1-XI2(J))*PHI(LN(I,1))+(1-XI2(J))*PHI(LN(I,2))+(1+XI2(J))*PHI(LN(I,3))+(-(1+XI2(J)))*PHI(LN(I,4)))
         DPDX(2)=0.25*(-(1-XI1(J))*PHI(LN(I,1))+(-(1+XI1(J)))*PHI(LN(I,2))+(1+XI1(J))*PHI(LN(I,3))+(1-XI1(J))*PHI(LN(I,4)))

         DPDS(1)=DPDX(1)/RM1
         DPDS(2)=-DPDX(1)*COS(THI)/RM1/SIN(THI)+DPDX(2)/RM2/SIN(THI)

		 S1(1)=PX1(1)/RM1
         S1(2)=PX1(2)/RM1
         S1(3)=PX1(3)/RM1

         S2(1)=NORM(2)*S1(3)-NORM(3)*S1(2)
         S2(2)=NORM(3)*S1(1)-NORM(1)*S1(3)
         S2(3)=NORM(1)*S1(2)-NORM(2)*S1(1)

!------THE FOLLOW IS FOR THE BOUNDARYT CONDITION
         DPPDX1=0.25*(-(1-XI2(J))*PPHI(LN(I,1))+(1-XI2(J))*PPHI(LN(I,2))+(1+XI2(J))*PPHI(LN(I,3))+(-(1+XI2(J)))*PPHI(LN(I,4)))
         DPPDX2=0.25*(-(1-XI1(J))*PPHI(LN(I,1))+(-(1+XI1(J)))*PPHI(LN(I,2))+(1+XI1(J))*PPHI(LN(I,3))+(1-XI1(J))*PPHI(LN(I,4)))

         DPPDS1=DPPDX1/RM1
         DPPDS2=-DPPDX1*COS(THI)/RM1/SIN(THI)+DPPDX2/RM2/SIN(THI)

         DPDS1X1=0.25*(-(1-XI2(J))*DPDS(1)+(1-XI2(J))*DPDS(1)+(1+XI2(J))*DPDS(1)+(-(1+XI2(J)))*DPDS(1))
         DPDS2X1=0.25*(-(1-XI2(J))*DPDS(2)+(1-XI2(J))*DPDS(2)+(1+XI2(J))*DPDS(2)+(-(1+XI2(J)))*DPDS(2))
         DPDS2X2=0.25*(-(1-XI1(J))*DPDS(2)+(-(1+XI1(J)))*DPDS(2)+(1+XI1(J))*DPDS(2)+(1-XI1(J))*DPDS(2))

         DPDS1S1=DPDS1X1/RM1
         DPDS2S2=-DPDS2X1*COS(THI)/RM1/SIN(THI)-DPDS2X2/RM2/SIN(THI)

         DPDNN=-(DPDS1S1+DPDS2S2)

         XMOD=DPDNN*NORM(1)+DPPDS1*S1(1)+DPPDS2*S2(1)
         YMOD=DPDNN*NORM(2)+DPPDS1*S1(2)+DPPDS2*S2(2)
         ZMOD=DPDNN*NORM(3)+DPPDS1*S1(3)+DPPDS2*S2(3)

!***** CALCULATE THE VELOCITY OF LATERAL PARTICLE AND PHI OVER T
       VELX(LN(I,J))=NORM(1)*PPHI(LN(I,J))+S1(1)*DPDS(1)+S2(1)*DPDS(2)
       VELY(LN(I,J))=NORM(2)*PPHI(LN(I,J))+S1(2)*DPDS(1)+S2(2)*DPDS(2)
       VELZ(LN(I,J))=0.0 !NORM(3)*PPHI(LN(I,J))+S1(3)*DPDS(1)+S2(3)*DPDS(2)

       ACCMO(LN(I,J))=VELX(LN(I,J))*XMOD+VELY(LN(I,J))*YMOD+VELZ(LN(I,J))*ZMOD
      END DO
      END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE BOUNDT(NNODE,NE,PHIT,PPHIT,ACCMO,DPDNN,ACC,WC)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,NNODE,NE(6)
	   REAL    R,ACC,WC
       REAL    NODE(NNODE,3),PHIT(NNODE),PPHIT(NNODE),ACCMO(NNODE),DPDNN(NNODE)

!   FREE SURFACE
       DO I=1,NE(1)
        PHIT(I)=PHIT(I)
       END DO

!   2ND SURFACE
       DO I=NE(1)+1,NE(2)
        PPHIT(I)=-WC*DPDNN(I)
       END DO

!   3RD SURFACE
       DO I=NE(2)+1,NE(3)
        PPHIT(I)=-WC*DPDNN(I)
       END DO

!   4TH SURFACE
       DO I=NE(3)+1,NE(4)
        PPHIT(I)=-WC*DPDNN(I)
       END DO

!   5TH SURFACE WAVE MAKER
       DO I=NE(4)+1,NE(5)
	   	R=SQRT(NODE(I,1)**2+NODE(I,3)**2)
        PPHIT(I)=-ACC*R-ACCMO(I)
       END DO

!   6TH SURFACE
       DO I=NE(5)+1,NNODE
        PPHIT(I)=0.0
      END DO

       RETURN
       END
!**********************************************************************
      SUBROUTINE TAYES2(NNODE,ENUM,LN,ME,NEL,PPHI,PHIT,PPHIT,VELX,VELY,VELZ,&
					   &NORMX,NORMY,NORMZ,S1X,S1Y,S1Z,S2X,S2Y,S2Z,RM1,RM2,THI,&
					   &DPDS1,DPDS2,X2ND,Y2ND,Z2ND,PHI2ND,GRAV)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL(A-H,O-Z)
      INTEGER NNODE,ENUM,LN(ENUM,4),ME(6),NEL(6)
      INTEGER IC(NNODE)
      REAL GRAV,D2X,D2Y,D2Z
      REAL PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),VELX(NNODE),VELY(NNODE),VELZ(NNODE) 
      REAL X2ND(NEL(1)),Y2ND(NEL(1)),Z2ND(NEL(1)),PHI2ND(NEL(1))
	  REAL NORMX(ME(1),4),NORMY(ME(1),4),NORMZ(ME(1),4)
	  REAL RM1(ME(1),4),RM2(ME(1),4),THI(ME(1),4)
      REAL S1X(ME(1),4),S1Y(ME(1),4),S1Z(ME(1),4)
      REAL S2X(ME(1),4),S2Y(ME(1),4),S2Z(ME(1),4)
      REAL DPDS1(ME(1),4),DPDS2(ME(1),4)
      REAL XI1(4),XI2(4)

       IC=0
       X2ND=0.0
       Y2ND=0.0
       Z2ND=0.0

      DO I=1,ME(1)
      DO J=1,4

!***** CALCULATE PARTIAL PHI PARTIAL T OVER PARTIAL S1 AND S2
         DPDTXI1=0.25*(-(1-XI2(J))*PHIT(LN(I,1))+(1-XI2(J))*PHIT(LN(I,2))+(1+XI2(J))*PHIT(LN(I,3))+(-(1+XI2(J)))*PHIT(LN(I,4)))
         DPDTXI2=0.25*(-(1-XI1(J))*PHIT(LN(I,1))+(-(1+XI1(J)))*PHIT(LN(I,2))+(1+XI1(J))*PHIT(LN(I,3))+(1-XI1(J))*PHIT(LN(I,4)))
!--------
         DPDS1XI1=0.25*(-(1-XI2(J))*DPDS1(I,1)+(1-XI2(J))*DPDS1(I,2)+(1+XI2(J))*DPDS1(I,3)+(-(1+XI2(J)))*DPDS1(I,4))
         DPDS1XI2=0.25*(-(1-XI1(J))*DPDS1(I,1)+(-(1+XI1(J)))*DPDS1(I,2)+(1+XI1(J))*DPDS1(I,3)+(1-XI1(J))*DPDS1(I,4))
!----------
         DPDS2XI1=0.25*(-(1-XI2(J))*DPDS2(I,1)+(1-XI2(J))*DPDS2(I,2)+(1+XI2(J))*DPDS2(I,3)+(-(1+XI2(J)))*DPDS2(I,4))
         DPDS2XI2=0.25*(-(1-XI1(J))*DPDS2(I,1)+(-(1+XI1(J)))*DPDS2(I,2)+(1+XI1(J))*DPDS2(I,3)+(1-XI1(J))*DPDS2(I,4))
!----------
         DPDNXI1=0.25*(-(1-XI2(J))*PPHI(LN(I,1))+(1-XI2(J))*PPHI(LN(I,2))+(1+XI2(J))*PPHI(LN(I,3))+(-(1+XI2(J)))*PPHI(LN(I,4)))
         DPDNXI2=0.25*(-(1-XI1(J))*PPHI(LN(I,1))+(-(1+XI1(J)))*PPHI(LN(I,2))+(1+XI1(J))*PPHI(LN(I,3))+(1-XI1(J))*PPHI(LN(I,4)))
!-----------
         DPDTS1=DPDTXI1/RM1(I,J)
         DPDTS2=-DPDTXI1*COS(THI(I,J))/RM1(I,J)/SIN(THI(I,J))-DPDTXI2/RM2(I,J)/SIN(THI(I,J))
!-----------
         DPDS1S1=DPDS1XI1/RM1(I,J)
         DPDS1S2=-DPDS1XI1*COS(THI(I,J))/RM1(I,J)/SIN(THI(I,J))-DPDS1XI2/RM2(I,J)/SIN(THI(I,J))
!-----------
         DPDS2S2=-DPDS2XI1*COS(THI(I,J))/RM1(I,J)/SIN(THI(I,J))-DPDS2XI2/RM2(I,J)/SIN(THI(I,J))
!-----------
         DPDNS1=DPDNXI1/RM1(I,J)
         DPDNS2=-DPDNXI1*COS(THI(I,J))/RM1(I,J)/SIN(THI(I,J))-DPDNXI2/RM2(I,J)/SIN(THI(I,J))
!-----------
         DPDNN=-(DPDS1S1+DPDS2S2)

!--- X DIRECTION
       PXS1=S1X(I,J)*DPDS1S1+S2X(I,J)*DPDS1S2+NORMX(I,J)*DPDNS1
       PXS2=S1X(I,J)*DPDS1S2+S2X(I,J)*DPDS2S2+NORMX(I,J)*DPDNS2
       PXN =S1X(I,J)*DPDNS1+S2X(I,J)*DPDNS2+NORMX(I,J)*DPDNN      

       PTX=S1X(I,J)*DPDTS1+S2X(I,J)*DPDTS2+NORMX(I,J)*PPHIT(LN(I,J))
       PXX=S1X(I,J)*PXS1+S2X(I,J)*PXS2+NORMX(I,J)*PXN
       PYX=S1Y(I,J)*PXS1+S2Y(I,J)*PXS2+NORMY(I,J)*PXN
       PZX=S1Z(I,J)*PXS1+S2Z(I,J)*PXS2+NORMZ(I,J)*PXN

       D2X=PTX+VELX(LN(I,J))*PXX+VELY(LN(I,J))*PYX+VELZ(LN(I,J))*PZX

!--- Y DIRECTION
       PYS1=S1Y(I,J)*DPDS1S1+S2Y(I,J)*DPDS1S2+NORMY(I,J)*DPDNS1
       PYS2=S1Y(I,J)*DPDS1S2+S2Y(I,J)*DPDS2S2+NORMY(I,J)*DPDNS2
       PYN =S1Y(I,J)*DPDNS1+S2Y(I,J)*DPDNS2+NORMY(I,J)*DPDNN      

       PTY=S1Y(I,J)*DPDTS1+S2Y(I,J)*DPDTS2+NORMY(I,J)*PPHIT(LN(I,J))
       PXY=S1X(I,J)*PYS1+S2X(I,J)*PYS2+NORMX(I,J)*PYN
       PYY=S1Y(I,J)*PYS1+S2Y(I,J)*PYS2+NORMY(I,J)*PYN
       PZY=S1Z(I,J)*PYS1+S2Z(I,J)*PYS2+NORMZ(I,J)*PYN

       D2Y=PTY+VELX(LN(I,J))*PXY+VELY(LN(I,J))*PYY+VELZ(LN(I,J))*PZY

!--- Z DIRECTION
       PZS1=S1Z(I,J)*DPDS1S1+S2Z(I,J)*DPDS1S2+NORMZ(I,J)*DPDNS1
       PZS2=S1Z(I,J)*DPDS1S2+S2Z(I,J)*DPDS2S2+NORMZ(I,J)*DPDNS2
       PZN =S1Z(I,J)*DPDNS1+S2Z(I,J)*DPDNS2+NORMZ(I,J)*DPDNN      

       PTZ=S1Z(I,J)*DPDTS1+S2Z(I,J)*DPDTS2+NORMZ(I,J)*PPHIT(LN(I,J))
       PXZ=S1X(I,J)*PZS1+S2X(I,J)*PZS2+NORMX(I,J)*PZN
       PYZ=S1Y(I,J)*PZS1+S2Y(I,J)*PZS2+NORMY(I,J)*PZN
       PZZ=S1Z(I,J)*PZS1+S2Z(I,J)*PZS2+NORMZ(I,J)*PZN

       D2Z=PTZ+VELX(LN(I,J))*PXZ+VELY(LN(I,J))*PYZ+VELZ(LN(I,J))*PZZ
!---------
      IC(LN(I,J))=IC(LN(I,J))+1
      X2ND(LN(I,J))=X2ND(LN(I,J))+D2X
      Y2ND(LN(I,J))=Y2ND(LN(I,J))+D2Y
      Z2ND(LN(I,J))=Z2ND(LN(I,J))+D2Z
      END DO
      END DO

      DO I=1,NEL(1)
       X2ND(I)=X2ND(I)/IC(I)
       Y2ND(I)=Y2ND(I)/IC(I)
       Z2ND(I)=Z2ND(I)/IC(I)
       PHI2ND(I)=VELX(I)*X2ND(I)+VELY(I)*Y2ND(I)+VELZ(I)*Z2ND(I)-GRAV*VELZ(I)
      END DO

      RETURN
      END
!********************************************************************
SUBROUTINE PRESSURE(ICON,THO,GRAV,DEP,NNODE,NE,NODE,PHIT,VELX,VELY,VELZ,PR,P_ATM)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,ICON,NNODE,NE(6)
      REAL DEP,THO,GRAV,P1,P2,P3,P_ATM
      REAL NODE(NNODE,2),PHIT(NNODE),VELX(NNODE),VELY(NNODE),VELZ(NNODE),PR(NNODE)
	  REAL CP1(NE(1)),CP2(NE(1)),CP3(NE(1)),CP(NE(1))

!----ATOM PRESSURE (BERNOULLI CONSTANT) ON THE FREE SURFACE
	DO I=ICON,ICON ! NE(1) !
	CP1(I)=THO*PHIT(I)
	CP2(I)=THO*0.5*(VELX(I)**2+VELY(I)**2+VELZ(I)**2)
	CP3(I)=THO*GRAV*(NODE(I,2)-DEP)
	CP(I)=CP1(I)+CP2(I)+CP3(I)
	ENDDO
	P_ATM=CP(ICON)

!----PRESSURE ON BOUNDARY
	DO I=NE(1)+1,NNODE
	P1=THO*PHIT(I)
	P2=THO*0.5*(VELX(I)**2+VELY(I)**2+VELZ(I)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP)
	PR(I)=P_ATM-(P1+P2+P3)
	END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE DOMAIN(NFIELD,NNODE,ENUM,NELEM,LN,DN,THO,GRAV,DEP,P_ATM,QPI,&
					   &NODE,NORMX,NORMY,NORMZ,JCB,PHI,PPHI,PHIT,PPHIT,&
					   &VELX,VELY,VELZ,PR,SH,W2,SHA1,SHA2,SHA3,SHA4)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,M
	  INTEGER NFIELD,NNODE,ENUM,NELEM(6,2),LN(ENUM,4),DN(2,(NELEM(1,1)-1)*(NELEM(1,2)-1))
	  REAL DX,DY,DZ,TEMP,R,P1,P2,P3
	  REAL THO,GRAV,DEP,P_ATM,QPI
	  REAL NODE(NNODE,3),NORMX(ENUM,7),NORMY(ENUM,7),NORMZ(ENUM,7),JCB(ENUM,7)
	  REAL PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),VELX(NNODE),VELY(NNODE),VELZ(NNODE),PR(NNODE)
	  REAL DNODE(NFIELD,3),DVX(NFIELD),DVY(NFIELD),DVZ(NFIELD),DPHIT(NFIELD),DPR(NFIELD)
	  REAL KER1X(NFIELD,NNODE),KER2X(NFIELD,NNODE),KER1Y(NFIELD,NNODE),KER2Y(NFIELD,NNODE)
	  REAL KER1Z(NFIELD,NNODE),KER2Z(NFIELD,NNODE),KER1T(NFIELD,NNODE),KER2T(NFIELD,NNODE)
	  REAL HX(4),GX(4),HY(4),GY(4),HZ(4),GZ(4),HT(4),GT(4)
      REAL SH(4,7),W2(7),SHA1(7),SHA2(7),SHA3(7),SHA4(7),XFUNC(7),YFUNC(7),ZFUNC(7)

!----CREATE DOMAIN POINT
	K=1
	DO I=1,(NELEM(1,1)-1)*(NELEM(1,2)-1)
	DX=(NODE(DN(1,I),1)-NODE(DN(2,I),1))/NELEM(2,2)
	DY=(NODE(DN(1,I),2)-NODE(DN(2,I),2))/NELEM(2,2)
	DZ=(NODE(DN(1,I),3)-NODE(DN(2,I),3))/NELEM(2,2)
		DO J=1,NELEM(2,2)-1
		DNODE(K,1)=NODE(DN(2,I),1)+DX*J
		DNODE(K,2)=NODE(DN(2,I),2)+DY*J
		DNODE(K,3)=NODE(DN(2,I),3)+DZ*J
		K=K+1
		END DO
	END DO

!*******CALCULATE X, Y, Z VELOCITY AND DPHIT BY BIE******************************************
	KER1X=0.0
	KER2X=0.0
	KER1Y=0.0
	KER2Y=0.0
	KER1Z=0.0
	KER2Z=0.0
	KER1T=0.0
	KER2T=0.0
	DVX=0.D0
	DVY=0.D0
	DVZ=0.D0
	DPHIT=0.D0
!-----THE SURFACE KERNELS
      DO I = 1,NFIELD
       DO J=1,ENUM
		DO M=1,7
		  XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)+SHA3(M)*NODE(LN(J,3),1)+SHA4(M)*NODE(LN(J,4),1)
		  YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)+SHA3(M)*NODE(LN(J,3),2)+SHA4(M)*NODE(LN(J,4),2)
		  ZFUNC(M)=SHA1(M)*NODE(LN(J,1),3)+SHA2(M)*NODE(LN(J,2),3)+SHA3(M)*NODE(LN(J,3),3)+SHA4(M)*NODE(LN(J,4),3)
		END DO
        DO K=1,4
         GX(K)=0.0
         HX(K)=0.0
         GY(K)=0.0
         HY(K)=0.0
         GZ(K)=0.0
         HZ(K)=0.0
         GT(K)=0.0
         HT(K)=0.0

		  DO M=1,7
		  	TEMP=JCB(J,M)*SH(K,M)*W2(M)
			R=SQRT((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2+(ZFUNC(M)-DNODE(I,3))**2)

			HX(K)=HX(K)+QPI*(3.0*(XFUNC(M)-DNODE(I,1))*&
					&((XFUNC(M)-DNODE(I,1))*NORMX(J,M)+(YFUNC(M)-DNODE(I,2))*NORMY(J,M)+(ZFUNC(M)-DNODE(I,3))*NORMZ(J,M))/R**5+&
					&NORMX(J,M)/R**3)*TEMP
			GX(K)=GX(K)+(-QPI*(XFUNC(M)-DNODE(I,1))/R**3)*TEMP

			HY(K)=HY(K)+QPI*(3.0*(YFUNC(M)-DNODE(I,2))*&
					&((XFUNC(M)-DNODE(I,1))*NORMX(J,M)+(YFUNC(M)-DNODE(I,2))*NORMY(J,M)+(ZFUNC(M)-DNODE(I,3))*NORMZ(J,M))/R**5+&
					&NORMY(J,M)/R**3)*TEMP
			GY(K)=GY(K)+(-QPI*(YFUNC(M)-DNODE(I,2))/R**3)*TEMP

			HZ(K)=HZ(K)+QPI*(3.0*(ZFUNC(M)-DNODE(I,3))*&
					&((XFUNC(M)-DNODE(I,1))*NORMX(J,M)+(YFUNC(M)-DNODE(I,2))*NORMY(J,M)+(ZFUNC(M)-DNODE(I,3))*NORMZ(J,M))/R**5+&
					&NORMZ(J,M)/R**3)*TEMP
			GZ(K)=GZ(K)+(-QPI*(ZFUNC(M)-DNODE(I,3))/R**3)*TEMP

            HT(K)=HT(K)+(-QPI/R**3*((XFUNC(M)-DNODE(I,1))*NORMX(J,M)+(YFUNC(M)-DNODE(I,2))*NORMY(J,M)+(ZFUNC(M)-DNODE(I,3))*NORMZ(J,M))*TEMP)
            GT(K)=GT(K)+QPI/R*TEMP
		  END DO
         KER1X(I,LN(J,K))=KER1X(I,LN(J,K))+HX(K)
         KER2X(I,LN(J,K))=KER2X(I,LN(J,K))+GX(K)
         KER1Y(I,LN(J,K))=KER1Y(I,LN(J,K))+HY(K)
         KER2Y(I,LN(J,K))=KER2Y(I,LN(J,K))+GY(K)
         KER1Z(I,LN(J,K))=KER1Z(I,LN(J,K))+HZ(K)
         KER2Z(I,LN(J,K))=KER2Z(I,LN(J,K))+GZ(K)
         KER1T(I,LN(J,K))=KER1T(I,LN(J,K))+HT(K)
         KER2T(I,LN(J,K))=KER2T(I,LN(J,K))+GT(K)
		END DO
	   END DO
      END DO
DVX=MATMUL(KER2X,PPHI)-MATMUL(KER1X,PHI)
DVY=MATMUL(KER2Y,PPHI)-MATMUL(KER1Y,PHI)
DVZ=MATMUL(KER2Z,PPHI)-MATMUL(KER1Z,PHI)
DPHIT=MATMUL(KER2T,PPHIT)-MATMUL(KER1T,PHIT)

!----CALCULATE PRESSURE DISTRIBUTION IN DOMAIN
	DO I=1,NFIELD
	P1=THO*DPHIT(I)
	P2=THO*0.5*(DVX(I)**2+DVY(I)**2+DVZ(I)**2)
	P3=THO*GRAV*(DNODE(I,2)-DEP)
	DPR(I)=P_ATM-(P1+P2+P3)
	END DO

	WRITE(11,'(5000(1X,F15.7))') NODE(:,1),DNODE(:,1)
	WRITE(11,'(5000(1X,F15.7))') NODE(:,2),DNODE(:,2)
	WRITE(11,'(5000(1X,F15.7))') NODE(:,3),DNODE(:,3)
	WRITE(11,'(5000(1X,F15.7))') VELX,DVX
	WRITE(11,'(5000(1X,F15.7))') VELY,DVY
	WRITE(11,'(5000(1X,F15.7))') VELZ,DVZ
	WRITE(11,'(5000(1X,F15.7))') PR,DPR

      RETURN
      END
!**********************************************************************
      SUBROUTINE GAUSS(W2,RT1,RT2)
!**********************************************************************
      REAL   W2(7),RT1(7),RT2(7)
!****** THE GAUSS POINT *****
        W2(1)=8./7
        RT1(1)=0
        RT2(1)=0
        W2(2)=5./9
        RT1(2)=(1./3)**0.5
        RT2(2)=(3./5)**0.5
        W2(3)=W2(2)
        RT1(3)=RT1(2)
        RT2(3)=-RT2(2)
        W2(4)=W2(2)
        RT1(4)=-RT1(2)
        RT2(4)=RT2(2)
        W2(5)=W2(2)
        RT1(5)=-RT1(2)
        RT2(5)=-RT2(2)
        W2(6)=20./63
        RT1(6)=(14./15)**0.5
        RT2(6)=0
        W2(7)=W2(6)
        RT1(7)=-RT1(6)
        RT2(7)=0
      RETURN
      END
!**********************************************************************
      SUBROUTINE SHAFUN(SHAPE,XI,ITA,K)
!**********************************************************************
      INTEGER K
      REAL SHAPE,XI,ITA
        SELECT CASE(K)
          CASE(1)
            SHAPE=0.25*(1-XI)*(1-ITA)
          CASE(2)
            SHAPE=0.25*(1+XI)*(1-ITA)
          CASE(3)
            SHAPE=0.25*(1+XI)*(1+ITA)
          CASE(4)
            SHAPE=0.25*(1-XI)*(1+ITA)
        END SELECT
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
