       PROGRAM Cs3_JUD
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION R(3)
c       EXTERNAL POTSO2
        INTEGER PERM, P
              DIMENSION Perm(3,6)

       OPEN(1, FILE="Cs3eq.dat", STATUS="UNKNOWN")
!        OPEN(2, FILE="config.dat", STATUS="UNKNOWN")

       THETA = 60/(180.0D0/3.141592653589793D0)!2.0826801557D0!3.14159D0!   ! Ângulo fixo (60 graus = 1.047 radianos)
       COS_THETA = COS(THETA)   ! Calcular o cosseno do ângulo fixo
!6,1400,0.01 demora
!       R(1) = 6.0D0
!       DO 100 I = 1, 1400
!          R(2) = 6.0D0
!          DO 200 J = 1, 1400
       R(1) = 6.0D0
       DO 100 I = 1, 140
          R(2) = 6.0D0
          DO 200 J = 1, 140

             ! Calcular R(3) pela lei dos cossenos
             R3_SQ = R(1)**2 + R(2)**2 - 2.0D0 * R(1) * R(2) * COS_THETA
             
             IF (R3_SQ .GT. 0.0D0) THEN
                R(3) = DSQRT(R3_SQ)  ! Garantir que não tente raiz de número negativo
             IF ((R(1) + R(2) .GT. R(3)) .AND. 
     &              (R(1) + R(3) .GT. R(2)) .AND. 
     &              (R(2) + R(3) .GT. R(1))) THEN
                V = POTen_cm(R)

                ! Verificar se V é um número válido e <= 0.2
                IF (V .EQ. V .AND. V .LE. 0.2D0) THEN
                   WRITE(1,*) R(1), R(2), R(3), V

C                DO P = 1, 6
C        WRITE(2,*) R(Perm(1,P)), R(Perm(2,P)), R(Perm(3,P)), V
C                END DO
                ENDIF
             ENDIF
             ENDIF
             
             R(2) = R(2) + 0.1D0
 200      CONTINUE
          R(1) = R(1) + 0.1D0
 100   CONTINUE
!              R(2) = R(2) + 0.01D0
! 200      CONTINUE
!          R(1) = R(1) + 0.01D0
! 100   CONTINUE

c       CLOSE(1)
       END

      DOUBLE PRECISION FUNCTION poten_cm(R)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(3), Q(3), J(3)
      EXTERNAL VPOT, VHH
      INTEGER I, II, JJ

      DOUBLE PRECISION FDAMP

      Cddd = 71.2D4      ! Constante C9 de Tang (1976), em a.u.
      BETA = 0.4D0       ! Parâmetro de amortecimento (a.u.^-1)

      R1 = R(1)
      R2 = R(2)
      R3 = R(3)

      ! lei dos cossenos
      COST1 = (R1**2 + R2**2 - R3**2) / (2.0D0 * R1 * R2)
      COST2 = (R3**2 + R2**2 - R1**2) / (2.0D0 * R2 * R3)
      COST3 = (R1**2 + R3**2 - R2**2) / (2.0D0 * R3 * R1)

      ! Parte angular do termo de Axilrod-Teller
      VISB0 = COST1 * COST2 * COST3
      VISB1 = 1.0D0 + 3.0D0 * VISB0

      ! Parte radial (produto de distâncias ao cubo)
      VISB2 = (R1 * R2 * R3)**3

      ! Produto das funções de amortecimento
      DAMP = FDAMP(R1, BETA) * FDAMP(R2, BETA) * FDAMP(R3, BETA)

      ! Termo de 3 corpos amortecido
      AT = (VISB1 / VISB2) * DAMP

      ! Energia total = pares + dispersão tripla amortecida
      poten_cm = Cddd * AT + (VPOT(R1) + VPOT(R2) + VPOT(R3))
      
     ! LEPS
c      DO I = 1, 3
c         Q(I) = 0.5D0 * (VPOT(R(I)) + VHH(R(I)))
c         J(I) = 0.5D0 * (VPOT(R(I)) - VHH(R(I)))
c      END DO
c
c      ! para i > j
c      SOMA2 = 0.0D0
c      DO II = 2, 3
c         DO JJ = 1, II - 1
c            DELTAJ = J(II) - J(JJ)
c            SOMA2 = SOMA2 + DELTAJ**2
c         END DO
c      END DO
c
c      ! VABC com sinal de menos (−)
c      VABCM = Q(1) + Q(2) + Q(3) - SQRT(0.5D0 * SOMA2)
c
c      ! Corrected line: Removed extra parenthesis
c      poten_cm =  VABCM ! CALIBRE * AT 

      !====================================================
      ! Função de amortecimento Tang-Toennies (ordem 6)
      !====================================================
      DOUBLE PRECISION FUNCTION FDAMP(R, BETA)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N, MAXN
      DOUBLE PRECISION R, SUM, EXPO

      MAXN = 6
      SUM = 0.D0
      DO N = 0, MAXN
         SUM = SUM + (BETA * R)**N / DFLOAT(FACT(N))
      END DO

      EXPO = EXP(-BETA * R)
      FDAMP = 1.D0 - EXPO * SUM
      RETURN
      END

      !====================================================
      ! Função fatorial inteira
      !====================================================
      INTEGER FUNCTION FACT(N)
      IMPLICIT INTEGER(A-Z)
      INTEGER I
      FACT = 1
      DO I = 2, N
         FACT = FACT * I
      END DO
      RETURN
      END

c
c      CALIBRE = 71.2E4 ! Tang(1976)
c      R1 = R(1)
c      R2 = R(2)
c      R3 = R(3)
c      
c      COST1 = (R1**2 + R2**2 - R3**2) / (2.0D0 * R1 * R2)
c      COST2 = (R3**2 + R2**2 - R1**2) / (2.0D0 * R2 * R3)
c      COST3 = (R1**2 + R3**2 - R2**2) / (2.0D0 * R3 * R1)
c
c      VISB0 = COST1 * COST2 * COST3
c      VISB1 = 1.0D0 + 3.0D0 * VISB0
c      VISB2 = (R1 * R2 * R3)**3
c      AT = VISB1 / VISB2

      DOUBLE PRECISION FUNCTION VPOT(X)
        IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION X, POTENCIAL
      DOUBLE PRECISION ALFA, BETA, GAMA, RE, DE, BE
      DOUBLE PRECISION POT_TRANS
      DOUBLE PRECISION, DIMENSION(2:6) :: C
      DOUBLE PRECISION, PARAMETER :: CONVERT = 4.556335252920D-6
      INTEGER IO, N
      
      POTENCIAL = 0.0D0
      DE = 0.01663000D0
      BE = 0.01174260D0 * CONVERT
      RE = 8.78530000D0
      ALFA = 0.30370000D0
      BETA = 0.10123000D0
      GAMA = 2.00000000D0

!      OPEN (UNIT=IO, FILE="coef.dat", STATUS="OLD", ACTION="READ")
!      READ(IO, *) C(2), C(3), C(4), C(5), C(6)
!      CLOSE(IO)

      C(2) = 4.2825162681146792D-002
      C(3) = 7.0581745889587710D-004
      C(4) = -2.4457218696157912D-002
      C(5) = -1.8062685619215989D-002
      C(6) = 9.3456819932705298D-003

      IF (X .LE. RE) THEN
         DO N = 2, 6
            POT_TRANS = (C(N) * (((1 + EXP(-2.0D0 * BETA * 
     &                  ((X - RE) / RE)))**N) * ((X - RE) / X)**N))
            POTENCIAL = POTENCIAL + POT_TRANS
         END DO
      ELSE
         POTENCIAL = DE * (((1 - EXP(-2.0D0 * ALFA * (X - RE))) / 
     &                     (1 + EXP(-GAMA * ALFA * (X - RE))))**2)
      END IF
      
      VPOT = POTENCIAL - DE
      RETURN
      end

!      DOUBLE PRECISION FUNCTION VHH(Y)
!      IMPLICIT REAL*8(A-H,O-Z)
!      REAL*8 Y, De, Re_au, Cc6, Cc8, Cc10, beta
!      
!      ! --- Parâmetros para Cs₂ ---
!      De = 0.001344D0           ! Profundidade do poço (Hartree)
!      Re_au = 11.91D0         ! Distância de equilíbrio (Bohr)
!      beta = 1.605D0             ! Parâmetro de decaimento
!      
!      ! Coeficientes de van der Waals convertidos (Hartree·Bohr^n)
!      Cc6 = 150.8D0             ! 3.31e7 cm⁻¹ → 150.8 Hartree·Bohr^6
!      Cc8 = 5920.5D0            ! 1.29962e9 cm⁻¹ → 5920.5 Hartree·Bohr^8
!      Cc10 = 234000.0D0         ! 5.136e10 cm⁻¹ → 234000 Hartree·Bohr^10
!      
!      IF (Y > Re_au) THEN
!          ! Termo repulsivo exponencial + van der Waals (atração)
!          VHH = De * EXP(-2.0D0 * beta * (Y - Re_au)) 
!     &      - (Cc6/Y**6 + Cc8/Y**8 + Cc10/Y**10)
!      ELSE
!          ! Termo repulsivo de curto alcance
!          VHH = De * EXP(-2.0D0 * beta * (Y - Re_au))
!      ENDIF
!      
!      VHH = VHH - De   ! Desloca para V(R) → 0 em R → ∞
!      RETURN
!      END

      DOUBLE PRECISION FUNCTION VHH(Y)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Y, De_eV, De, Re_au, we_cm, we, wexe_cm, wexe
      REAL*8 Be, Rc_au, h, c_light, pi, m_atom_u, m_atom_kg
      REAL*8 mu, x, a0, a1, a2, alpha, c, b, beta, aux
      REAL*8 au_to_m

      ! --- Constantes físicas ---
      h = 6.62607015D-34            ! Planck (J.s)
      c_light = 2.99792458D10        ! velocidade da luz (cm/s)
      pi = 3.141592653589793D0
      au_to_m = 5.29177210903D-11    ! 1 Bohr = x m

      ! --- Conversões ---
      m_atom_u = 132.905D0
      m_atom_kg = m_atom_u * 1.66053906660D-27
      mu = m_atom_kg / 2.0D0         ! massa reduzida
      De_eV = 0.03658 !Li(2007), 0.388D0
      De = De_eV * 0.0367493D0       ! eV → Hartree
      we_cm = 11.58d0!Li(2007)!28.2D0
      we = we_cm * 4.556335D-6       ! cm^-1 → Hartree
      wexe_cm = 0.046D0
      wexe = wexe_cm * 4.556335D-6    ! cm^-1 → Hartree
      Re_au = 11.91!Li(2007) 10.64D0                ! raio de equilíbrio em Bohr
      Rc_au = Re_au

      ! --- Cálculo de Be (em Hartree) ---
      Be = h / (8.0D0 * pi**2 * c_light * mu * (Re_au * au_to_m)**2)
      Be = Be * 4.556335D-6          ! cm^-1 → Hartree

      ! --- Coeficientes de Dunham ---
      a0 = we**2 / (4.0D0 * Be)
      aux = Be / (we * wexe)
      alpha = 6.0D0 * Be * wexe * (DSQRT(aux) - aux)
      a1 = -1.0D0 - alpha * (we / (6.0D0 * Be))**2
      a2 = (5.0D0/4.0D0) * a1**2 - (2.0D0/3.0D0) * (we * wexe / Be)
      c = 1.0D0 + a1 * DSQRT(De / a0)
      b = 2.0D0 * ((7.0D0/12.0D0) - (De * a2 / a0)) / c
      beta = we / (2.0D0 * Rc_au * DSQRT(Be * De))

      ! --- Substituição x = beta * (R - Re) ---
      x = beta * (Y - Re_au)

      ! --- Potencial VHH (em Hartree) ---
      VHH = De * ((1.0D0 - DEXP(-x))**2 + 
     &      (1.0D0 + b * x) * c * x**3 * DEXP(-2.0D0 * x)) - De

      RETURN
      END