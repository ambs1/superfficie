PROGRAM CALCULOPOTENCIAL
IMPLICIT REAL*8(A-H,O-Z)
!============= POTENCIAL DE UM SISTEMA DIATÔMICO =============

INTEGER :: I, io, PONTOS
DOUBLE PRECISION :: PRIMEIRO, ULTIMO, PASSO
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R, V

! Configurações de cálculo
PONTOS = 100000
PRIMEIRO = 0.0D0
ULTIMO = 55.5D0

ALLOCATE(R(PONTOS))
ALLOCATE(V(PONTOS))
PASSO = (ULTIMO - PRIMEIRO) / (PONTOS - 1)

! Calcula o potencial V(R)
DO I = 1, PONTOS
    R(I) = PRIMEIRO + (I - 1) * PASSO
    V(I) = VHH(R(I))
END DO

! Salva em arquivo
OPEN (newunit=io, file="HHsdun.dat", status="replace", action="write")
DO I = 1, PONTOS
    WRITE(io, *) R(I), V(I)
END DO

CONTAINS

FUNCTION VHH(R)
  IMPLICIT REAL*8(A-H,O-Z)
  real(8), intent(in) :: R
  real(8) :: VHH
  real(8) :: De_eV, De, Re_au, we_cm, wexe_cm, f
  real(8) :: we, wexe, mu, Be, Rc_au
  real(8) :: h, c_light, pi, m_atom_u, m_atom_kg, au_to_m
  real(8) :: x, a0, a1, a2, alpha, c, b, beta, aux

  ! --- Constantes físicas ---
  h = 6.62607015d-34            ! Planck (J.s)
  c_light = 2.99792458d10       ! velocidade da luz (cm/s)
  pi = 3.141592653589793d0
  au_to_m = 5.29177210903d-11   ! 1 Bohr = x m

  ! --- Conversões ---
  ! 1 u = 1.66053906660e-27 kg
  ! 1 cm^-1 = 4.556335e-6 Hartree
  ! 1 eV = 0.0367493 Hartree

  ! --- Parâmetros fornecidos ---
  m_atom_u = 132.905d0
  m_atom_kg = m_atom_u * 1.66053906660d-27
  mu = m_atom_kg / 2.0d0                      ! massa reduzida
  De_eV = 0.388d0
  De = De_eV * 0.0367493d0                    ! eV → Hartree
  we_cm = 28.2d0
  we = we_cm * 4.556335d-6                    ! cm^-1 → Hartree
  wexe_cm = 0.046d0
  wexe = wexe_cm * 4.556335d-6                ! cm^-1 → Hartree
  Re_au = 10.64d0                             ! raio de equilíbrio em Bohr
  Rc_au = Re_au
  f = 1.0d0                                   ! parâmetro empírico

  ! --- Cálculo de Be (em Hartree) ---
  Be = h / (8.0d0 * pi**2 * c_light * mu * (Re_au * au_to_m)**2)  ! cm^-1
  Be = Be * 4.556335d-6  ! cm^-1 → Hartree

  ! --- Coeficientes de Dunham ---
  a0 = we**2 / (4.0d0 * Be)

  ! --- Novo cálculo de alpha baseado em Be, we e wexe ---
  aux = Be / (we * wexe)
  alpha = 6.0d0 * Be * wexe * (sqrt(aux) - aux)  ! alpha em Hartree

  a1 = -1.0d0 - alpha * (we / (6.0d0 * Be))**2
  a2 = (5.0d0/4.0d0) * a1**2 - (2.0d0/3.0d0) * (we * wexe / Be)

  c = 1.0d0 + a1 * sqrt(De / a0)
  b = 2.0d0 * ((7.0d0/12.0d0) - (De * a2 / a0)) / c
  beta = we / (2.0d0 * Rc_au * sqrt(Be * De))

  ! --- Substituição x = beta * (R - Re) ---
  x = beta * (R - Re_au)
!  WRITE(*,*) De, 'Be =', Be, 'we =', we, 'wexe =', wexe, 'Re =', Re_au
!  WRITE(*,*) 'alpha =', alpha, 'a0 =', a0, 'a1 =', a1, 'a2 =', a2

  ! --- Potencial VHH (em Hartree) ---
  VHH = De * ((1.0d0 - exp(-x))**2 + (1.0d0 + b * x) * c * x**3 * exp(-2.0d0 * x)) - De

END FUNCTION VHH

END PROGRAM CALCULOPOTENCIAL

