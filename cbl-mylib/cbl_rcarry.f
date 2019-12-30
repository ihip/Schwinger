C     ******************************************************************
      SUBROUTINE cbl_rcarry(RVEC, LENV, ADD)
C     ------------------------------------------------------------------
C     ******************************************************************
C     **                   **                                         **
C     **    R C A R R Y    **     Implemented by I. Hip, Mar 1992     **
C     **                   **     Amended by CBL,        Jan 1994     **
C     **                   **                                         **
C     ******************************************************************
C     Algorithm due to: G. Marsaglia and A. Zaman
C     Implementation based on paper by Fred James: "A Review of pseudo-
C     random number generators", Comput. Phys. Commun. 60 (1990) 329.
C     including the correction due to Luescher: exchange I24, J24 
C     the original version in line UNI=SEEDS(I24)-SEEDS(J24)-CARRY
C     ------------------------------------------------------------------
C     ATTENTION: program does NOT vectorize
C     ------------------------------------------------------------------
C     Program delivers vector r(n) of random numbers in the
C     interval c < r =< c+1.0
C     Time requirement decreases with length of vector and is
C     in microseconds
C     for vector length    DN4000(FPX) DN10000 HP400   MIPS4000
C               1           58.2        9.3      23.0
C              50           10.0        1.4       3.2
C            1000            6.3        1.1       2.7
C           50000            6.8        1.3       2.2  0.32
C     ------------------------------------------------------------------
      LOGICAL NOTINIT
      DIMENSION RVEC(LENV)
      COMMON /R_CARRY/ SEEDS(48), CARRY, NOTINIT
      DATA NOTINIT/.TRUE./
      PARAMETER (TWOP24 = 16777216.0)
      PARAMETER (TWOM24 = 1.0/TWOP24)
 
      EPSADD = ADD + TWOM24
      IF(NOTINIT) CALL cbl_rc_seed(54217137)

      IF(LENV .LT. 25) THEN

        DO 10 I = 1, LENV
          UNI = SEEDS(I + 14) - SEEDS(I) - CARRY
          IF (UNI.LT.0.0) THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
          ELSE
            CARRY = 0
          ENDIF
          SEEDS(I + 24) = UNI
          RVEC(I) = UNI + EPSADD
10      CONTINUE
        DO 20 I = 1, 24
          SEEDS(I) = SEEDS(I + LENV)
20      CONTINUE

      ELSE

        DO 30 I = 1, 10
          UNI = SEEDS(I + 14) - SEEDS(I) - CARRY
          RVEC(I) = UNI
          IF (UNI .LT. 0.0) THEN
            CARRY = TWOM24
            RVEC(I) = UNI + 1.0
          ELSE
            CARRY = 0.0
          ENDIF
30      CONTINUE

        DO 31 I = 11, 24
          UNI = RVEC(I - 10) - SEEDS(I) - CARRY
          RVEC(I) = UNI
          IF (UNI .LT. 0.0) THEN
            CARRY = TWOM24
            RVEC(I) = UNI + 1.0
          ELSE
            CARRY = 0.0
          ENDIF
31      CONTINUE

        DO 40 I = 25, LENV
          UNI = RVEC(I - 10) - RVEC(I - 24) - CARRY
          RVEC(I) = UNI
          IF (UNI .LT. 0.0) THEN
            CARRY = TWOM24
            RVEC(I) = UNI  + 1.0
          ELSE
            CARRY = 0.0
          ENDIF

 40      CONTINUE

c
c       the preceeding verison is somewhat faster than this one:
c        DO 40 I = 25, LENV
c          UNI = RVEC(I - 10) - RVEC(I - 24) - CARRY
c          IF (UNI .LT. 0.0) THEN
c            UNI = UNI + 1.0
c            CARRY = TWOM24
c          ELSE
c            CARRY = 0.0
c          ENDIF
c          RVEC(I) = UNI
c40      CONTINUE
 
        DO 50 I = 1, 24
          SEEDS(I) = RVEC(LENV - 24 + I)
50      CONTINUE
 
        DO 60 I = 1, LENV
          RVEC(I) = RVEC(I) + EPSADD
60      CONTINUE

      ENDIF
 
      RETURN
      END
 
      SUBROUTINE cbl_rc_seed(IJKL)
C     ------------------------------------------------------------------
C     One of the two ways to initialize rcarry - this (a bit strange)
C     random number generator is given in paper by F. James, and is
C     used for initializing SEEDS table. CARRY is set to 0.0.
C     James:
C       "The input value should be in the range: 0 <= IJKL <= 900 000 000"
C     ------------------------------------------------------------------
      LOGICAL NOTINIT
      COMMON /R_CARRY/ SEEDS(48), CARRY, NOTINIT
 
      IF((IJKL .LT. 0) .OR. (IJKL .GT. 900000000)) IJKL = 54217137
C         - if out of range, set to default value
      IJ = IJKL/30082
      KL = IJKL - 30082 * IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177) + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      DO 2 II = 24, 1, -1
C         - this implementation of rcarry uses SEEDS in inverse order,
C           to make output fully compatible to original, initialization
C           is also made inverse (from 24 to 1, step -1)
        S = 0.0
        T = 0.5
        DO 3 JJ = 1, 24
          M = MOD(MOD(I * J, 179) * K, 179)
          I = J
          J = K
          K = M
          L = MOD(53 * L + 1, 169)
          IF (MOD(L * M, 64).GE.32) S = S + T
          T = 0.5 * T
3       CONTINUE
        SEEDS(II) = S
2     CONTINUE
      CARRY = 0.0
C         - James: "the starting value of CARRY must also be
C           initialized, but it can be started with zero."
      NOTINIT = .FALSE.
      RETURN
      END
 
      SUBROUTINE cbl_rc_set(R)
C     ----------------------------------------------------------------
C     SEEDS initialization with R(1:24), CARRY = R(25)
C     ----------------------------------------------------------------
      LOGICAL NOTINIT
      DIMENSION R(25)
      COMMON /R_CARRY/ SEEDS(48), CARRY, NOTINIT
 
      DO 70 I = 1, 24
        SEEDS(I) = R(I)
70    CONTINUE
      CARRY = R(25)
      NOTINIT = .FALSE.
 
      RETURN
      END
 
      SUBROUTINE cbl_rc_get(R)
C     ----------------------------------------------------------------
C     Get last SEEDS in R(1:24) and last CARRY in R(25)
C     ----------------------------------------------------------------
      LOGICAL NOTINIT
      DIMENSION R(25)
      COMMON /R_CARRY/ SEEDS(48), CARRY, NOTINIT
 
      DO 80 I = 1, 24
        R(I) = SEEDS(I)
80    CONTINUE
      R(25) = CARRY
 
      RETURN
      END
 
      SUBROUTINE cbl_rc_log(RVEC, LENV, ADD)
C     ----------------------------------------------------------------
C     Same as rcarry, but it gives the natural logarithms of
C     the random numbers as output
C     ----------------------------------------------------------------
      LOGICAL NOTINIT
      DIMENSION RVEC(LENV)
      COMMON /R_CARRY/ SEEDS(48), CARRY, NOTINIT
      DATA NOTINIT/.TRUE./
      PARAMETER (TWOP24 = 16777216.0)
      PARAMETER (TWOM24 = 1.0/TWOP24)
 
      EPSADD = ADD + TWOM24
      IF(NOTINIT) CALL cbl_rc_seed(54217137)

      IF(LENV .LT. 25) THEN

        DO 10 I = 1, LENV
          UNI = SEEDS(I + 14) - SEEDS(I) - CARRY
          IF (UNI.LT.0.0) THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
          ELSE
            CARRY = 0
          ENDIF
          SEEDS(I + 24) = UNI
          RVEC(I) = ALOG(UNI + EPSADD)
10      CONTINUE
        DO 20 I = 1, 24
          SEEDS(I) = SEEDS(I + LENV)
20      CONTINUE

      ELSE

        DO 30 I = 1, 10
          UNI = SEEDS(I + 14) - SEEDS(I) - CARRY
          RVEC(I) = UNI
          IF (UNI .LT. 0.0) THEN
            CARRY = TWOM24
            RVEC(I) = UNI  + 1.0
          ELSE
            CARRY = 0.0
          ENDIF

30      CONTINUE
        DO 31 I = 11, 24
          UNI =  RVEC(I - 10) - SEEDS(I) - CARRY
          RVEC(I) = UNI
          IF (UNI .LT. 0.0) THEN
            CARRY = TWOM24
            RVEC(I) = UNI  + 1.0
          ELSE
            CARRY = 0.0
          ENDIF
31      CONTINUE
 
        DO 40 I = 25, LENV
          UNI =  RVEC(I - 10) - RVEC(I - 24) - CARRY
          RVEC(I) = UNI
          IF (UNI .LT. 0.0) THEN
            CARRY = TWOM24
            RVEC(I) = UNI  + 1.0
          ELSE
            CARRY = 0.0
          ENDIF
40      CONTINUE
 
        DO 50 I = 1, 24
          SEEDS(I) = RVEC(LENV - 24 + I)
50      CONTINUE
 
        DO 60 I = 1, LENV
          RVEC(I) = ALOG(RVEC(I) + EPSADD)
60      CONTINUE
      ENDIF
 
      RETURN
      END
