       PROGRAM Cs3_JUD
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION R(3), Angulos(1000000), ValoresV(1000000)
       DIMENSION R1_vals(1000000), R2_vals(1000000), R3_vals(1000000)
       INTEGER I, J, K, P, II, JJ, M, Count
       INTEGER Perm
       DIMENSION Perm(3,6)
       REAL*8 TEMP_V, TEMP_A, TEMP_R1, TEMP_R2, TEMP_R3, V, A
       REAL*8 poten_cm
       EXTERNAL poten_cm

       OPEN(1, FILE="Cs3.dat", STATUS="UNKNOWN")
       OPEN(2, FILE="config.dat", STATUS="UNKNOWN")
       OPEN(3, FILE="angulo.dat", STATUS="UNKNOWN")

       THETA = 3.141592653589793D0
       COS_THETA = COS(THETA)

       ! Permutacoes possiveis para verificar simetria
       DATA Perm / 1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1 /

       Count = 0

       DO I = 1, 80
          R(1) = 8.0 + I * 0.2D0
          DO J = 1, 80
             R(2) =8.0 + J * 0.2D0
             DO K = 1, 80
                R(3) = 8.0 + K * 0.2D0

                ! Verifica a condição do triângulo (0 A 180 GRAUS)
                IF ((R(1) + R(2) .GE. R(3)) .AND.  
     &              (R(1) + R(3) .GE. R(2)) .AND.  
     &              (R(2) + R(3) .GE. R(1))) THEN

                   V = poten_cm(R)
                   IF (V .EQ. V .AND. V .LE. 0.02D0) THEN
                      WRITE(1,*) R(1), R(2), R(3), V

                      ! Salva todas as permutacoes no arquivo config.dat
                      DO P = 1, 6
        WRITE(2,*) R(Perm(1,P)), R(Perm(2,P)), R(Perm(3,P)), V
                      END DO

                      ! Calcula um dos ângulos
        A = DACOS((R(2)**2 + R(3)**2 - R(1)**2) / (2.0D0 * R(2) * R(3)))
                      
                      Count = Count + 1
                      Angulos(Count) = A
                      ValoresV(Count) = V
                      R1_vals(Count) = R(1)
                      R2_vals(Count) = R(2)
                      R3_vals(Count) = R(3)
                   ENDIF
                ENDIF

             END DO
          END DO
       END DO

       ! Ordena os valores pelo potencial V antes de salvar
       DO II = 1, Count-1
          DO JJ = II+1, Count
             IF (ValoresV(JJ) .LT. ValoresV(II)) THEN
                TEMP_V = ValoresV(II)
                ValoresV(II) = ValoresV(JJ)
                ValoresV(JJ) = TEMP_V

                TEMP_A = Angulos(II)
                Angulos(II) = Angulos(JJ)
                Angulos(JJ) = TEMP_A

                TEMP_R1 = R1_vals(II)
                R1_vals(II) = R1_vals(JJ)
                R1_vals(JJ) = TEMP_R1

                TEMP_R2 = R2_vals(II)
                R2_vals(II) = R2_vals(JJ)
                R2_vals(JJ) = TEMP_R2

                TEMP_R3 = R3_vals(II)
                R3_vals(II) = R3_vals(JJ)
                R3_vals(JJ) = TEMP_R3
             ENDIF
          END DO
       END DO

       ! Escreve os ângulos ordenados por V
       DO I = 1, Count
       WRITE(3,*) ValoresV(I),(Angulos(I)*180.0D0/3.141592653589793D0),
     &   R1_vals(I), R2_vals(I), R3_vals(I)
       END DO

       CLOSE(1)
       CLOSE(2)
       CLOSE(3)
       END

      DOUBLE PRECISION FUNCTION poten_cm(R)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(3), Q(3), J(3)
      EXTERNAL VPOT, VHH
      INTEGER I, II, JJ

      CALIBRE = 71.2E4 ! Tang(1976)
      R1 = R(1)
      R2 = R(2)
      R3 = R(3)
      
      COST1 = (R1**2 + R2**2 - R3**2) / (2.0D0 * R1 * R2)
      COST2 = (R3**2 + R2**2 - R1**2) / (2.0D0 * R2 * R3)
      COST3 = (R1**2 + R3**2 - R2**2) / (2.0D0 * R3 * R1)

      VISB0 = COST1 * COST2 * COST3
      VISB1 = 1.0D0 + 3.0D0 * VISB0
      VISB2 = (R1 * R2 * R3)**3
      AT = VISB1 / VISB2
      
      ! Calcula Q(i) e J(i)
      DO I = 1, 3
         Q(I) = 0.5D0 * (VPOT(R(I)) + VHH(R(I)))
         J(I) = 0.5D0 * (VPOT(R(I)) - VHH(R(I)))
      END DO

      ! Soma (J_i - J_j)^2 para i > j
      SOMA2 = 0.0D0
      DO II = 2, 3
         DO JJ = 1, II - 1
            DELTAJ = J(II) - J(JJ)
            SOMA2 = SOMA2 + DELTAJ**2
         END DO
      END DO

      ! Valor de VABC com sinal de menos (−)
      VABCM = Q(1) + Q(2) + Q(3) - SQRT(0.5D0 * SOMA2)

      ! Corrected line: Removed extra parenthesis
      poten_cm =  VABCM ! CALIBRE * AT 
c      poten_cm = VABCM + (VPOT(R1) + VPOT(R2) + VPOT(R3))
c     poten_cm = VABCM + CALIBRE * AT + (VPOT(R1) + VPOT(R2) + VPOT(R3))
      RETURN
      END

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
         
      
                 
          