SUBROUTINE Lscales_from_fft (dataset, nxpoints, nzpoints, xlength, zlength)

! Code to estimate the lengthscale of a field of data by the FFT

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  fft_worklen_x,          &
  fft_init_x,             &
  fft_wsave_x,            &
  fft_work_x,             &
  fft_worklen_z,          &
  fft_init_z,             &
  fft_wsave_z,            &
  fft_work_z

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN)    :: dataset(1:nxpoints,1:nzpoints)
INTEGER,      INTENT(IN)    :: nxpoints, nzpoints
REAL(ZREAL8), INTENT(INOUT) :: xlength, zlength

REAL(ZREAL8)                :: Dummy_x(1:nxpoints), Dummy_z(1:nzpoints)
REAL(ZREAL8)                :: mean_len_x(1:nzpoints), mean_len_z(1:nxpoints), mean_len
REAL(ZREAL8)                :: this_len, nxpoints_r, nzpoints_r, norm, re, im, weight
INTEGER                     :: ierr, k, xx, zz


nxpoints_r = REAL(nxpoints)
nzpoints_r = REAL(nzpoints)


! Deal with the x-direction
! -------------------------

! Initialize the FFT in the x-direction if not already done
IF (.NOT.fft_init_x) THEN
  CALL rfft1i (nxpoints, fft_wsave_x, fft_worklen_x, ierr)
  fft_init_x = .TRUE.
END IF

DO zz = 1, nzpoints
  ! Do fft in x for this z value
  Dummy_x(1:nxpoints) = dataset(1:nxpoints,zz)
  CALL rfft1f (nxpoints, 1,                                            &
               Dummy_x(1:nxpoints), nxpoints,                          &
               fft_wsave_x, fft_worklen_x, fft_work_x, nxpoints, ierr)

  ! Treat fft^2 as a PDF - find the mean lengthscale
  mean_len = 0.0
  norm     = 0.0
  DO k = 2, nxpoints/2 + 1
    this_len = nxpoints_r / REAL(k)
    ! See /home/ross/Info/MyComputerDocs/f90.html.LyXconv/f90.html for information on structure of FFT array
    re       = Dummy_x(2*k-2)
    IF (k == nxpoints/2 + 1) THEN
      im = 0.0
    ELSE
      im = Dummy_x(2*k-1)
    END IF
    weight   = re*re + im*im
    mean_len = mean_len + weight * this_len
    norm     = norm + weight
  END DO

  IF (norm < 0.0000000001) PRINT *, 'WARNING FROM HORIZONTAL FFT ROUTINE - ZERO AT ', zz
  mean_len_x(zz) = mean_len / (2.0 * norm)    ! Half wavelength
END DO

xlength = SUM(mean_len_x(1:nzpoints)) / nzpoints_r


! Deal with the z-direction
! -------------------------

! Initialize the FFT in the z-direction if not already done
IF (.NOT.fft_init_z) THEN
  CALL rfft1i (nzpoints, fft_wsave_z, fft_worklen_z, ierr)
  fft_init_z = .TRUE.
END IF

DO xx = 1, nxpoints
  ! Do fft in z for this x value
  Dummy_z(1:nzpoints) = dataset(xx,1:nzpoints)
  CALL rfft1f (nzpoints, 1,                                            &
               Dummy_z(1:nzpoints), nzpoints,                          &
               fft_wsave_z, fft_worklen_z, fft_work_z, nzpoints, ierr)

  ! Treat fft^2 as a PDF - find the mean lengthscale
  mean_len = 0.0
  norm     = 0.0
  DO k = 2, nzpoints/2 + 1
    this_len = nzpoints_r / REAL(k)
    ! See /home/ross/Info/MyComputerDocs/f90.html.LyXconv/f90.html for information on structure of FFT array
    re       = Dummy_z(2*k-2)
    IF (k == nzpoints/2 + 1) THEN
      im = 0.0
    ELSE
      im = Dummy_z(2*k-1)
    END IF
    weight   = re*re + im*im
    mean_len = mean_len + weight * this_len
    norm     = norm + weight
  END DO

  IF (norm < 0.0000000001) PRINT *, 'WARNING FROM VERTICAL  FFT ROUTINE - ZERO AT ', xx
  mean_len_z(xx) = mean_len / (2.0 * norm)    ! Half wavelength
END DO

zlength = SUM(mean_len_z(1:nxpoints)) / nxpoints_r

END SUBROUTINE Lscales_from_fft
