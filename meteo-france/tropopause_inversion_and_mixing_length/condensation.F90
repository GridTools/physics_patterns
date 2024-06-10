!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
    SUBROUTINE CONDENSATION(D, CST, &
        &PZZ, PT, PRC_IN, PRI_IN,  &
        &OSIGMAS, &
        &PLV, PLS, PCPH)

!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CST,            ONLY: CST_t

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
TYPE(CST_t),                  INTENT(IN)    :: CST
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRC_IN ! grid scale r_c mixing ratio (kg/kg) in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRI_IN ! grid scale r_i (kg/kg) in input

LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                   ! or that from turbulence scheme
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLV    ! Latent heat L_v
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLS
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCPH   ! Specific heat C_ph


!*       0.2   Declarations of local variables :
!
INTEGER :: JIJ, JK, JKP, JKM                    ! loop index
INTEGER :: IKTB, IKTE, IKB, IKE, IKL, IIJB, IIJE
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTLK           ! work arrays for T_l and total water mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT) :: ZL             ! length scale
INTEGER, DIMENSION(D%NIJT)  :: ITPL             ! top levels of troposphere
REAL,    DIMENSION(D%NIJT)  :: ZTMIN            ! minimum Temp. related to ITPL
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZLV, ZLS, ZCPD

REAL :: ZLL, DZZ, ZZZ                           ! used for length scales

REAL, DIMENSION(D%NIJT) :: ZDZ
REAL :: ZPRIFACT


REAL,PARAMETER :: ZL0     = 600.        ! tropospheric length scale


------------------------------------
!
!
!
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
IIJB=D%NIJB
IIJE=D%NIJE

ZPRIFACT = 1.    ! Initialize value

IF(PRESENT(PLV) .AND. PRESENT(PLS)) THEN
    ZLV(:,:)=PLV(:,:)
    ZLS(:,:)=PLS(:,:)

IF(PRESENT(PCPH)) THEN
    ZCPD(:,:)=PCPH(:,:)

! Preliminary calculations needed for computing the "turbulent part" of Sigma_s
IF ( .NOT. OSIGMAS ) THEN
    DO JK=IKTB,IKTE
        DO JIJ=IIJB,IIJE
! store temperature at saturation
            ZTLK(JIJ,JK) = PT(JIJ,JK) - ZLV(JIJ,JK)*PRC_IN(JIJ,JK)/ZCPD(JIJ,JK) &
                 - ZLS(JIJ,JK)*PRI_IN(JIJ,JK)/ZCPD(JIJ,JK)*ZPRIFACT
        END DO
    END DO
! Determine tropopause/inversion  height from minimum temperature
    ITPL(:)  = IKB+IKL
    ZTMIN(:) = 400.
    DO JK = IKTB+1,IKTE-1
        DO JIJ=IIJB,IIJE
            IF ( PT(JIJ,JK) < ZTMIN(JIJ) ) THEN
                ZTMIN(JIJ) = PT(JIJ,JK)
                ITPL(JIJ) = JK
            ENDIF
        END DO
    END DO
    ! Set the mixing length scale
    ZL(:,IKB) = 20.
    DO JK = IKB+IKL,IKE,IKL
        DO JIJ=IIJB,IIJE
            ! free troposphere
            ZL(JIJ,JK) = ZL0
            ZZZ =  PZZ(JIJ,JK) -  PZZ(JIJ,IKB)
            JKP = ITPL(JIJ)
            ! approximate length for boundary-layer
            IF ( ZL0 > ZZZ ) ZL(JIJ,JK) = ZZZ
            ! gradual decrease of length-scale near and above tropopause
            IF ( ZZZ > 0.9*(PZZ(JIJ,JKP)-PZZ(JIJ,IKB)) ) &
                ZL(JIJ,JK) = .6 * ZL(JIJ,JK-IKL)
        END DO
    END DO
END IF