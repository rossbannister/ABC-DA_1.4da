SUBROUTINE Helmholtz_inv (psi, chi, u, v)

! Code to perform the Helmholtz (find psi and chi from u and v)

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  nlongs,                 &
  nlevs,                  &
  dx,                     &
  recipdx,                &
  fourpi2,                &
  fft_worklen_x,          &
  fft_init_x,             &
  fft_wsave_x,            &
  fft_work_x

IMPLICIT NONE

REAL(ZREAL8),   INTENT(INOUT) :: psi(1:nlongs,1:nlevs)
REAL(ZREAL8),   INTENT(INOUT) :: chi(1:nlongs,1:nlevs)
REAL(ZREAL8),   INTENT(IN)    :: u(0:nlongs+1,1:nlevs)
REAL(ZREAL8),   INTENT(IN)    :: v(0:nlongs+1,1:nlevs)

INTEGER                       :: x, z, k
REAL(ZREAL8)                  :: ddx(1:nlongs), facs(2:nlongs/2+1)
REAL(ZREAL8)                  :: factor, km1
INTEGER                       :: ierr, real_index, imag_index

! See Eq (39) of the model paper


! Initialize the FFT in the x-direction if not already done
IF (.NOT.fft_init_x) THEN
  CALL rfft1i (nlongs, fft_wsave_x, fft_worklen_x, ierr)
  fft_init_x = .TRUE.
END IF

! Setup some factors needed multiple times
factor = -1.0 * dx * dx * REAL(nlongs) * REAL(nlongs) / fourpi2
DO k = 2, nlongs/2+1
  km1        = REAL(k-1)
  facs(k)    = factor / (km1 * km1)
END DO

! 1. Find chi
! -----------
DO z = 1, nlevs
  ! Compute du/dx
  DO x = 1, nlongs
    ddx(x) = (u(x,z) - u(x-1,z)) * recipdx
  END DO
  ! Fourier transform this
  CALL rfft1f (nlongs, 1,                                            &
               ddx(1:nlongs), nlongs,                                &
               fft_wsave_x, fft_worklen_x, fft_work_x, nlongs, ierr)
  ! Multiply by -dx^2 nlongs^2 / (k-1)^2 (but set constant term to zero)
  ddx(1) = 0.0
  DO k = 2, nlongs/2
    real_index      = 2*k-2
    imag_index      = 2*k-1
    ddx(real_index) = ddx(real_index) * facs(k)
    ddx(imag_index) = ddx(imag_index) * facs(k)
  END DO
  ! Smallest scale
  ddx(nlongs) = ddx(nlongs) * facs(nlongs/2+1)
  ! Inverse Fourier transform this
  CALL rfft1b (nlongs, 1,                                           &
               ddx(1:nlongs), nlongs,                               &
               fft_wsave_x, fft_worklen_x, fft_work_x, nlongs, ierr)
  ! ddx no longer contains du/dx, but chi
  chi(1:nlongs,z) = ddx(1:nlongs)
END DO



! 2. Find psi
! -----------
DO z = 1, nlevs
  ! Compute dv/dx
  DO x = 1, nlongs
    ddx(x) = (v(x+1,z) - v(x,z)) * recipdx
  END DO
  ! Fourier transform this
  CALL rfft1f (nlongs, 1,                                            &
               ddx(1:nlongs), nlongs,                                &
               fft_wsave_x, fft_worklen_x, fft_work_x, nlongs, ierr)
  ! Multiply by -dx^2 nlongs^2 / (k-1)^2 (but set constant term to zero)
  ddx(1) = 0.0
  DO k = 2, nlongs/2
    real_index      = 2*k-2
    imag_index      = 2*k-1
    ddx(real_index) = ddx(real_index) * facs(k)
    ddx(imag_index) = ddx(imag_index) * facs(k)
  END DO
  ! Smallest scale
  ddx(nlongs) = ddx(nlongs) * facs(nlongs/2+1)
  ! Inverse Fourier transform this
  CALL rfft1b (nlongs, 1,                                           &
               ddx(1:nlongs), nlongs,                               &
               fft_wsave_x, fft_worklen_x, fft_work_x, nlongs, ierr)
  ! ddx no longer contains dv/dx, but psi
  psi(1:nlongs,z) = ddx(1:nlongs)
END DO


END SUBROUTINE Helmholtz_inv
