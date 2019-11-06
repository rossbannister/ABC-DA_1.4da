SUBROUTINE U_v (ControlVar1, ControlVar2, CVT)

! Code to perform the vertical cvt: ControlVar2 = U_v ControlVar1

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  CV_type,                &
  CVT_type,               &
  nlongs,                 &
  nlevs


IMPLICIT NONE

TYPE(CV_type),  INTENT(IN)    :: ControlVar1
TYPE(CV_type),  INTENT(INOUT) :: ControlVar2
TYPE(CVT_type), INTENT(IN)    :: CVT

INTEGER                       :: x, z, k, real_index, imag_index, last_index
TYPE(CV_type)                 :: Interim


IF (CVT % CVT_order == 1) THEN
  ! Original Met Office transform order
  ! -----------------------------------
  IF (CVT % CVT_vert_opt_sym == 2) THEN
    ! Do symmetric transform - transform from levels to vert mode space
    DO x = 1, nlongs
      Interim % v1(x,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode1(1:nlevs,1:nlevs,1)), ControlVar1 % v1(x,1:nlevs))
      Interim % v2(x,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode2(1:nlevs,1:nlevs,1)), ControlVar1 % v2(x,1:nlevs))
      Interim % v3(x,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode3(1:nlevs,1:nlevs,1)), ControlVar1 % v3(x,1:nlevs))
      Interim % v4(x,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode4(1:nlevs,1:nlevs,1)), ControlVar1 % v4(x,1:nlevs))
      Interim % v5(x,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode5(1:nlevs,1:nlevs,1)), ControlVar1 % v5(x,1:nlevs))
      Interim % v6(x,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode6(1:nlevs,1:nlevs,1)), ControlVar1 % v6(x,1:nlevs))
    END DO
  ELSE
    ! Do non-symmetric transform - fields are already in vert mode space
    Interim % v1(1:nlongs,1:nlevs) = ControlVar1 % v1(1:nlongs,1:nlevs)
    Interim % v2(1:nlongs,1:nlevs) = ControlVar1 % v2(1:nlongs,1:nlevs)
    Interim % v3(1:nlongs,1:nlevs) = ControlVar1 % v3(1:nlongs,1:nlevs)
    Interim % v4(1:nlongs,1:nlevs) = ControlVar1 % v4(1:nlongs,1:nlevs)
    Interim % v5(1:nlongs,1:nlevs) = ControlVar1 % v5(1:nlongs,1:nlevs)
    Interim % v6(1:nlongs,1:nlevs) = ControlVar1 % v6(1:nlongs,1:nlevs)
  END IF
  ! Multiply by the evs (these are actually the square-roots)
  ! These are factors in vert mode space
  DO x = 1, nlongs
    Interim % v1(x,1:nlevs) = CVT % VertEV1(1:nlevs,1) * Interim % v1(x,1:nlevs)
    Interim % v2(x,1:nlevs) = CVT % VertEV2(1:nlevs,1) * Interim % v2(x,1:nlevs)
    Interim % v3(x,1:nlevs) = CVT % VertEV3(1:nlevs,1) * Interim % v3(x,1:nlevs)
    Interim % v4(x,1:nlevs) = CVT % VertEV4(1:nlevs,1) * Interim % v4(x,1:nlevs)
    Interim % v5(x,1:nlevs) = CVT % VertEV5(1:nlevs,1) * Interim % v5(x,1:nlevs)
    Interim % v6(x,1:nlevs) = CVT % VertEV6(1:nlevs,1) * Interim % v6(x,1:nlevs)
  END DO
  ! Do transform from vert mode space to level space
  DO x = 1, nlongs
    ControlVar2 % v1(x,1:nlevs) = MATMUL(CVT % VertMode1(1:nlevs,1:nlevs,1), Interim % v1(x,1:nlevs))
    ControlVar2 % v2(x,1:nlevs) = MATMUL(CVT % VertMode2(1:nlevs,1:nlevs,1), Interim % v2(x,1:nlevs))
    ControlVar2 % v3(x,1:nlevs) = MATMUL(CVT % VertMode3(1:nlevs,1:nlevs,1), Interim % v3(x,1:nlevs))
    ControlVar2 % v4(x,1:nlevs) = MATMUL(CVT % VertMode4(1:nlevs,1:nlevs,1), Interim % v4(x,1:nlevs))
    ControlVar2 % v5(x,1:nlevs) = MATMUL(CVT % VertMode5(1:nlevs,1:nlevs,1), Interim % v5(x,1:nlevs))
    ControlVar2 % v6(x,1:nlevs) = MATMUL(CVT % VertMode6(1:nlevs,1:nlevs,1), Interim % v6(x,1:nlevs))
  END DO




ELSE
  ! Reversed transform order (different vert transform for each horz scale)
  ! This is the only other option when vertical transform is called.
  ! ----------------------------------------------------------------------------------
  IF (CVT % CVT_vert_opt_sym == 2) THEN
    ! Do symmetric transform - transform from levels to vert mode space

    ! Deal with the largest scale first
    Interim % v1(1,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode1(1:nlevs,1:nlevs,1)), ControlVar1 % v1(1,1:nlevs))
    Interim % v2(1,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode2(1:nlevs,1:nlevs,1)), ControlVar1 % v2(1,1:nlevs))
    Interim % v3(1,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode3(1:nlevs,1:nlevs,1)), ControlVar1 % v3(1,1:nlevs))
    Interim % v4(1,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode4(1:nlevs,1:nlevs,1)), ControlVar1 % v4(1,1:nlevs))
    Interim % v5(1,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode5(1:nlevs,1:nlevs,1)), ControlVar1 % v5(1,1:nlevs))
    Interim % v6(1,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode6(1:nlevs,1:nlevs,1)), ControlVar1 % v6(1,1:nlevs))
    ! Deal with the bulk of the scales
    DO k = 2, nlongs/2
      real_index = 2*k-2
      imag_index = 2*k-1
      Interim % v1(real_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode1(1:nlevs,1:nlevs,k)), ControlVar1 % v1(real_index,1:nlevs))
      Interim % v1(imag_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode1(1:nlevs,1:nlevs,k)), ControlVar1 % v1(imag_index,1:nlevs))
      Interim % v2(real_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode2(1:nlevs,1:nlevs,k)), ControlVar1 % v2(real_index,1:nlevs))
      Interim % v2(imag_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode2(1:nlevs,1:nlevs,k)), ControlVar1 % v2(imag_index,1:nlevs))
      Interim % v3(real_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode3(1:nlevs,1:nlevs,k)), ControlVar1 % v3(real_index,1:nlevs))
      Interim % v3(imag_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode3(1:nlevs,1:nlevs,k)), ControlVar1 % v3(imag_index,1:nlevs))
      Interim % v4(real_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode4(1:nlevs,1:nlevs,k)), ControlVar1 % v4(real_index,1:nlevs))
      Interim % v4(imag_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode4(1:nlevs,1:nlevs,k)), ControlVar1 % v4(imag_index,1:nlevs))
      Interim % v5(real_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode5(1:nlevs,1:nlevs,k)), ControlVar1 % v5(real_index,1:nlevs))
      Interim % v5(imag_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode5(1:nlevs,1:nlevs,k)), ControlVar1 % v5(imag_index,1:nlevs))
      Interim % v6(real_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode6(1:nlevs,1:nlevs,k)), ControlVar1 % v6(real_index,1:nlevs))
      Interim % v6(imag_index,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode6(1:nlevs,1:nlevs,k)), ControlVar1 % v6(imag_index,1:nlevs))
    END DO
    ! Deal with the smallest scale
    last_index = nlongs/2+1
    Interim % v1(nlongs,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode1(1:nlevs,1:nlevs,last_index)), ControlVar1 % v1(nlongs,1:nlevs))
    Interim % v2(nlongs,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode2(1:nlevs,1:nlevs,last_index)), ControlVar1 % v2(nlongs,1:nlevs))
    Interim % v3(nlongs,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode3(1:nlevs,1:nlevs,last_index)), ControlVar1 % v3(nlongs,1:nlevs))
    Interim % v4(nlongs,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode4(1:nlevs,1:nlevs,last_index)), ControlVar1 % v4(nlongs,1:nlevs))
    Interim % v5(nlongs,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode5(1:nlevs,1:nlevs,last_index)), ControlVar1 % v5(nlongs,1:nlevs))
    Interim % v6(nlongs,1:nlevs) = MATMUL(TRANSPOSE(CVT % VertMode6(1:nlevs,1:nlevs,last_index)), ControlVar1 % v6(nlongs,1:nlevs))

  ELSE
    ! Do non-symmetric transform - fields are already in vert mode space
    Interim % v1(1:nlongs,1:nlevs) = ControlVar1 % v1(1:nlongs,1:nlevs)
    Interim % v2(1:nlongs,1:nlevs) = ControlVar1 % v2(1:nlongs,1:nlevs)
    Interim % v3(1:nlongs,1:nlevs) = ControlVar1 % v3(1:nlongs,1:nlevs)
    Interim % v4(1:nlongs,1:nlevs) = ControlVar1 % v4(1:nlongs,1:nlevs)
    Interim % v5(1:nlongs,1:nlevs) = ControlVar1 % v5(1:nlongs,1:nlevs)
    Interim % v6(1:nlongs,1:nlevs) = ControlVar1 % v6(1:nlongs,1:nlevs)
  END IF

  ! Multiply by the evs (these are actually the square-roots)
  ! These are factors in vert mode space
  ! Deal with the largest scale first
  Interim % v1(1,1:nlevs) = CVT % VertEV1(1:nlevs,1) * Interim % v1(1,1:nlevs)
  Interim % v2(1,1:nlevs) = CVT % VertEV2(1:nlevs,1) * Interim % v2(1,1:nlevs)
  Interim % v3(1,1:nlevs) = CVT % VertEV3(1:nlevs,1) * Interim % v3(1,1:nlevs)
  Interim % v4(1,1:nlevs) = CVT % VertEV4(1:nlevs,1) * Interim % v4(1,1:nlevs)
  Interim % v5(1,1:nlevs) = CVT % VertEV5(1:nlevs,1) * Interim % v5(1,1:nlevs)
  Interim % v6(1,1:nlevs) = CVT % VertEV6(1:nlevs,1) * Interim % v6(1,1:nlevs)
  ! Deal with the bulk of the scales
  DO k = 2, nlongs/2
    real_index = 2*k-2
    imag_index = 2*k-1
    Interim % v1(real_index,1:nlevs) = CVT % VertEV1(1:nlevs,k) * Interim % v1(real_index,1:nlevs)
    Interim % v1(imag_index,1:nlevs) = CVT % VertEV1(1:nlevs,k) * Interim % v1(imag_index,1:nlevs)
    Interim % v2(real_index,1:nlevs) = CVT % VertEV2(1:nlevs,k) * Interim % v2(real_index,1:nlevs)
    Interim % v2(imag_index,1:nlevs) = CVT % VertEV2(1:nlevs,k) * Interim % v2(imag_index,1:nlevs)
    Interim % v3(real_index,1:nlevs) = CVT % VertEV3(1:nlevs,k) * Interim % v3(real_index,1:nlevs)
    Interim % v3(imag_index,1:nlevs) = CVT % VertEV3(1:nlevs,k) * Interim % v3(imag_index,1:nlevs)
    Interim % v4(real_index,1:nlevs) = CVT % VertEV4(1:nlevs,k) * Interim % v4(real_index,1:nlevs)
    Interim % v4(imag_index,1:nlevs) = CVT % VertEV4(1:nlevs,k) * Interim % v4(imag_index,1:nlevs)
    Interim % v5(real_index,1:nlevs) = CVT % VertEV5(1:nlevs,k) * Interim % v5(real_index,1:nlevs)
    Interim % v5(imag_index,1:nlevs) = CVT % VertEV5(1:nlevs,k) * Interim % v5(imag_index,1:nlevs)
    Interim % v6(real_index,1:nlevs) = CVT % VertEV6(1:nlevs,k) * Interim % v6(real_index,1:nlevs)
    Interim % v6(imag_index,1:nlevs) = CVT % VertEV6(1:nlevs,k) * Interim % v6(imag_index,1:nlevs)
  END DO
  ! Deal with the smallest scale
  last_index = nlongs/2+1
  Interim % v1(nlongs,1:nlevs) = CVT % VertEV1(1:nlevs,last_index) * Interim % v1(nlongs,1:nlevs)
  Interim % v2(nlongs,1:nlevs) = CVT % VertEV2(1:nlevs,last_index) * Interim % v2(nlongs,1:nlevs)
  Interim % v3(nlongs,1:nlevs) = CVT % VertEV3(1:nlevs,last_index) * Interim % v3(nlongs,1:nlevs)
  Interim % v4(nlongs,1:nlevs) = CVT % VertEV4(1:nlevs,last_index) * Interim % v4(nlongs,1:nlevs)
  Interim % v5(nlongs,1:nlevs) = CVT % VertEV5(1:nlevs,last_index) * Interim % v5(nlongs,1:nlevs)
  Interim % v6(nlongs,1:nlevs) = CVT % VertEV6(1:nlevs,last_index) * Interim % v6(nlongs,1:nlevs)

  ! Do transform from vert mode space to level space
  ! Deal with the largest scale first
  ControlVar2 % v1(1,1:nlevs) = MATMUL(CVT % VertMode1(1:nlevs,1:nlevs,1), Interim % v1(1,1:nlevs))
  ControlVar2 % v2(1,1:nlevs) = MATMUL(CVT % VertMode2(1:nlevs,1:nlevs,1), Interim % v2(1,1:nlevs))
  ControlVar2 % v3(1,1:nlevs) = MATMUL(CVT % VertMode3(1:nlevs,1:nlevs,1), Interim % v3(1,1:nlevs))
  ControlVar2 % v4(1,1:nlevs) = MATMUL(CVT % VertMode4(1:nlevs,1:nlevs,1), Interim % v4(1,1:nlevs))
  ControlVar2 % v5(1,1:nlevs) = MATMUL(CVT % VertMode5(1:nlevs,1:nlevs,1), Interim % v5(1,1:nlevs))
  ControlVar2 % v6(1,1:nlevs) = MATMUL(CVT % VertMode6(1:nlevs,1:nlevs,1), Interim % v6(1,1:nlevs))
  ! Deal with the bulk of the scales
  DO k = 2, nlongs/2
    real_index = 2*k-2
    imag_index = 2*k-1
    ControlVar2 % v1(real_index,1:nlevs) = MATMUL(CVT % VertMode1(1:nlevs,1:nlevs,k), Interim % v1(real_index,1:nlevs))
    ControlVar2 % v1(imag_index,1:nlevs) = MATMUL(CVT % VertMode1(1:nlevs,1:nlevs,k), Interim % v1(imag_index,1:nlevs))
    ControlVar2 % v2(real_index,1:nlevs) = MATMUL(CVT % VertMode2(1:nlevs,1:nlevs,k), Interim % v2(real_index,1:nlevs))
    ControlVar2 % v2(imag_index,1:nlevs) = MATMUL(CVT % VertMode2(1:nlevs,1:nlevs,k), Interim % v2(imag_index,1:nlevs))
    ControlVar2 % v3(real_index,1:nlevs) = MATMUL(CVT % VertMode3(1:nlevs,1:nlevs,k), Interim % v3(real_index,1:nlevs))
    ControlVar2 % v3(imag_index,1:nlevs) = MATMUL(CVT % VertMode3(1:nlevs,1:nlevs,k), Interim % v3(imag_index,1:nlevs))
    ControlVar2 % v4(real_index,1:nlevs) = MATMUL(CVT % VertMode4(1:nlevs,1:nlevs,k), Interim % v4(real_index,1:nlevs))
    ControlVar2 % v4(imag_index,1:nlevs) = MATMUL(CVT % VertMode4(1:nlevs,1:nlevs,k), Interim % v4(imag_index,1:nlevs))
    ControlVar2 % v5(real_index,1:nlevs) = MATMUL(CVT % VertMode5(1:nlevs,1:nlevs,k), Interim % v5(real_index,1:nlevs))
    ControlVar2 % v5(imag_index,1:nlevs) = MATMUL(CVT % VertMode5(1:nlevs,1:nlevs,k), Interim % v5(imag_index,1:nlevs))
    ControlVar2 % v6(real_index,1:nlevs) = MATMUL(CVT % VertMode6(1:nlevs,1:nlevs,k), Interim % v6(real_index,1:nlevs))
    ControlVar2 % v6(imag_index,1:nlevs) = MATMUL(CVT % VertMode6(1:nlevs,1:nlevs,k), Interim % v6(imag_index,1:nlevs))
  END DO
  ! Deal with the smallest scale
  last_index = nlongs/2+1
  ControlVar2 % v1(nlongs,1:nlevs) = MATMUL(CVT % VertMode1(1:nlevs,1:nlevs,last_index), Interim % v1(nlongs,1:nlevs))
  ControlVar2 % v2(nlongs,1:nlevs) = MATMUL(CVT % VertMode2(1:nlevs,1:nlevs,last_index), Interim % v2(nlongs,1:nlevs))
  ControlVar2 % v3(nlongs,1:nlevs) = MATMUL(CVT % VertMode3(1:nlevs,1:nlevs,last_index), Interim % v3(nlongs,1:nlevs))
  ControlVar2 % v4(nlongs,1:nlevs) = MATMUL(CVT % VertMode4(1:nlevs,1:nlevs,last_index), Interim % v4(nlongs,1:nlevs))
  ControlVar2 % v5(nlongs,1:nlevs) = MATMUL(CVT % VertMode5(1:nlevs,1:nlevs,last_index), Interim % v5(nlongs,1:nlevs))
  ControlVar2 % v6(nlongs,1:nlevs) = MATMUL(CVT % VertMode6(1:nlevs,1:nlevs,last_index), Interim % v6(nlongs,1:nlevs))

END IF

END SUBROUTINE U_v
