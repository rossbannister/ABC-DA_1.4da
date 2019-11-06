SUBROUTINE U_h_adj (ControlVar1, ControlVar2, CVT)

! Code to perform the adjoint horizontal cvt: ControlVar1 = U_h* ControlVar2

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  CV_type,                &
  CVT_type,               &
  nlongs,                 &
  nlevs


IMPLICIT NONE

TYPE(CV_type),  INTENT(INOUT) :: ControlVar1
TYPE(CV_type),  INTENT(IN)    :: ControlVar2
TYPE(CVT_type), INTENT(IN)    :: CVT

INTEGER                       :: z, k, real_index, imag_index, last_index
TYPE(CV_type)                 :: Interim


! Do transform from real space to spectral space

CALL fft_real2spec (ControlVar2 % v1(1:nlongs,1:nlevs), Interim % v1(1:nlongs,1:nlevs))
CALL fft_real2spec (ControlVar2 % v2(1:nlongs,1:nlevs), Interim % v2(1:nlongs,1:nlevs))
CALL fft_real2spec (ControlVar2 % v3(1:nlongs,1:nlevs), Interim % v3(1:nlongs,1:nlevs))
CALL fft_real2spec (ControlVar2 % v4(1:nlongs,1:nlevs), Interim % v4(1:nlongs,1:nlevs))
CALL fft_real2spec (ControlVar2 % v5(1:nlongs,1:nlevs), Interim % v5(1:nlongs,1:nlevs))
CALL fft_real2spec (ControlVar2 % v6(1:nlongs,1:nlevs), Interim % v6(1:nlongs,1:nlevs))


! Multiply by the evs (these are actually the square-roots)
! These are factors in spectral space
! Deal with the largest scale first
Interim % v1(1,1:nlevs) = Interim % v1(1,1:nlevs) * CVT % HorizEV1(1,1:nlevs)
Interim % v2(1,1:nlevs) = Interim % v2(1,1:nlevs) * CVT % HorizEV2(1,1:nlevs)
Interim % v3(1,1:nlevs) = Interim % v3(1,1:nlevs) * CVT % HorizEV3(1,1:nlevs)
Interim % v4(1,1:nlevs) = Interim % v4(1,1:nlevs) * CVT % HorizEV4(1,1:nlevs)
Interim % v5(1,1:nlevs) = Interim % v5(1,1:nlevs) * CVT % HorizEV5(1,1:nlevs)
Interim % v6(1,1:nlevs) = Interim % v6(1,1:nlevs) * CVT % HorizEV6(1,1:nlevs)
! Deal with the bulk of the scales
DO k = 2, nlongs/2
  real_index = 2*k-2
  imag_index = 2*k-1
  Interim % v1(real_index,1:nlevs) = Interim % v1(real_index,1:nlevs) * CVT % HorizEV1(k,1:nlevs)
  Interim % v1(imag_index,1:nlevs) = Interim % v1(imag_index,1:nlevs) * CVT % HorizEV1(k,1:nlevs)
  Interim % v2(real_index,1:nlevs) = Interim % v2(real_index,1:nlevs) * CVT % HorizEV2(k,1:nlevs)
  Interim % v2(imag_index,1:nlevs) = Interim % v2(imag_index,1:nlevs) * CVT % HorizEV2(k,1:nlevs)
  Interim % v3(real_index,1:nlevs) = Interim % v3(real_index,1:nlevs) * CVT % HorizEV3(k,1:nlevs)
  Interim % v3(imag_index,1:nlevs) = Interim % v3(imag_index,1:nlevs) * CVT % HorizEV3(k,1:nlevs)
  Interim % v4(real_index,1:nlevs) = Interim % v4(real_index,1:nlevs) * CVT % HorizEV4(k,1:nlevs)
  Interim % v4(imag_index,1:nlevs) = Interim % v4(imag_index,1:nlevs) * CVT % HorizEV4(k,1:nlevs)
  Interim % v5(real_index,1:nlevs) = Interim % v5(real_index,1:nlevs) * CVT % HorizEV5(k,1:nlevs)
  Interim % v5(imag_index,1:nlevs) = Interim % v5(imag_index,1:nlevs) * CVT % HorizEV5(k,1:nlevs)
  Interim % v6(real_index,1:nlevs) = Interim % v6(real_index,1:nlevs) * CVT % HorizEV6(k,1:nlevs)
  Interim % v6(imag_index,1:nlevs) = Interim % v6(imag_index,1:nlevs) * CVT % HorizEV6(k,1:nlevs)
END DO
! Deal with the smallest scale
last_index = nlongs/2+1
Interim % v1(nlongs,1:nlevs) = Interim % v1(nlongs,1:nlevs) * CVT % HorizEV1(last_index,1:nlevs)
Interim % v2(nlongs,1:nlevs) = Interim % v2(nlongs,1:nlevs) * CVT % HorizEV2(last_index,1:nlevs)
Interim % v3(nlongs,1:nlevs) = Interim % v3(nlongs,1:nlevs) * CVT % HorizEV3(last_index,1:nlevs)
Interim % v4(nlongs,1:nlevs) = Interim % v4(nlongs,1:nlevs) * CVT % HorizEV4(last_index,1:nlevs)
Interim % v5(nlongs,1:nlevs) = Interim % v5(nlongs,1:nlevs) * CVT % HorizEV5(last_index,1:nlevs)
Interim % v6(nlongs,1:nlevs) = Interim % v6(nlongs,1:nlevs) * CVT % HorizEV6(last_index,1:nlevs)

ControlVar1 % v1(1:nlongs,1:nlevs) = Interim % v1(1:nlongs,1:nlevs)
ControlVar1 % v2(1:nlongs,1:nlevs) = Interim % v2(1:nlongs,1:nlevs)
ControlVar1 % v3(1:nlongs,1:nlevs) = Interim % v3(1:nlongs,1:nlevs)
ControlVar1 % v4(1:nlongs,1:nlevs) = Interim % v4(1:nlongs,1:nlevs)
ControlVar1 % v5(1:nlongs,1:nlevs) = Interim % v5(1:nlongs,1:nlevs)
ControlVar1 % v6(1:nlongs,1:nlevs) = Interim % v6(1:nlongs,1:nlevs)


END SUBROUTINE U_h_adj
