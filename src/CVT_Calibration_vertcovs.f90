SUBROUTINE CVT_Calibration_vertcovs (CVT, pert, div_cons)

! If div_cons = 0.0, use the perturbation to contribute to the vertical covariances

! If div_cons > 0.0, divide vertical covs by div_cons
!                    (to compute vertical covariances)

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    CVT_type,                &
    CV_type,                 &
    nlongs, nlevs,           &
    small

IMPLICIT NONE

TYPE(CVT_type), INTENT(INOUT)   :: CVT
TYPE(CV_type),  INTENT(IN)      :: pert
REAL(ZREAL8),   INTENT(IN)      :: div_cons

REAL(ZREAL8)                    :: mul_cons
INTEGER                         :: x, z1, z2, k, lim
INTEGER                         :: real_index, imag_index


IF (DABS(div_cons) < small) THEN
  ! Contribute to the vertical covariances for this perturbation

  SELECT CASE (CVT % CVT_order)

  CASE (1)
    ! Original MetO order (horiz co-ordinate of pert is longitude)
    ! Use the VertMode structure for storage
    DO x = 1, nlongs
      DO z1 = 1, nlevs
        DO z2 = 1, nlevs
          CVT % VertMode1(z1,z2,1) = CVT % VertMode1(z1,z2,1) + &
                                     pert % v1(x, z1) *         &
                                     pert % v1(x, z2)
          CVT % VertMode2(z1,z2,1) = CVT % VertMode2(z1,z2,1) + &
                                     pert % v2(x, z1) *         &
                                     pert % v2(x, z2)
          CVT % VertMode3(z1,z2,1) = CVT % VertMode3(z1,z2,1) + &
                                     pert % v3(x, z1) *         &
                                     pert % v3(x, z2)
          CVT % VertMode4(z1,z2,1) = CVT % VertMode4(z1,z2,1) + &
                                     pert % v4(x, z1) *         &
                                     pert % v4(x, z2)
          CVT % VertMode5(z1,z2,1) = CVT % VertMode5(z1,z2,1) + &
                                     pert % v5(x, z1) *         &
                                     pert % v5(x, z2)
          CVT % VertMode6(z1,z2,1) = CVT % VertMode6(z1,z2,1) + &
                                     pert % v6(x, z1) *         &
                                     pert % v6(x, z2)
        END DO
      END DO
    END DO


  CASE (2)
    ! Reverse MetO order (horiz co-ordinate of pert is wavenumber)
    ! Use the VertMode structure for storage

    ! The data in pert are complex in general
    ! The k=1 component though is real
    k = 1
    DO z1 = 1, nlevs
      DO z2 = 1, nlevs
        CVT % VertMode1(z1,z2,k) = CVT % VertMode1(z1,z2,k) + &
                                   pert % v1(k, z1) *         &
                                   pert % v1(k, z2)
        CVT % VertMode2(z1,z2,k) = CVT % VertMode2(z1,z2,k) + &
                                   pert % v2(k, z1) *         &
                                   pert % v2(k, z2)
        CVT % VertMode3(z1,z2,k) = CVT % VertMode3(z1,z2,k) + &
                                   pert % v3(k, z1) *         &
                                   pert % v3(k, z2)
        CVT % VertMode4(z1,z2,k) = CVT % VertMode4(z1,z2,k) + &
                                   pert % v4(k, z1) *         &
                                   pert % v4(k, z2)
        CVT % VertMode5(z1,z2,k) = CVT % VertMode5(z1,z2,k) + &
                                   pert % v5(k, z1) *         &
                                   pert % v5(k, z2)
        CVT % VertMode6(z1,z2,k) = CVT % VertMode6(z1,z2,k) + &
                                   pert % v6(k, z1) *         &
                                   pert % v6(k, z2)
      END DO
    END DO

    ! Perts contain complex data, but covs itself must be real
    DO k = 2, nlongs/2
      real_index = 2*k-2
      imag_index = 2*k-1
      DO z1 = 1, nlevs
        DO z2 = 1, nlevs
          CVT % VertMode1(z1,z2,k) = CVT % VertMode1(z1,z2,k) + &
                                     pert % v1(real_index,z1) * &
                                     pert % v1(real_index,z2) + &
                                     pert % v1(imag_index,z1) * &
                                     pert % v1(imag_index,z2)
          CVT % VertMode2(z1,z2,k) = CVT % VertMode2(z1,z2,k) + &
                                     pert % v2(real_index,z1) * &
                                     pert % v2(real_index,z2) + &
                                     pert % v2(imag_index,z1) * &
                                     pert % v2(imag_index,z2)
          CVT % VertMode3(z1,z2,k) = CVT % VertMode3(z1,z2,k) + &
                                     pert % v3(real_index,z1) * &
                                     pert % v3(real_index,z2) + &
                                     pert % v3(imag_index,z1) * &
                                     pert % v3(imag_index,z2)
          CVT % VertMode4(z1,z2,k) = CVT % VertMode4(z1,z2,k) + &
                                     pert % v4(real_index,z1) * &
                                     pert % v4(real_index,z2) + &
                                     pert % v4(imag_index,z1) * &
                                     pert % v4(imag_index,z2)
          CVT % VertMode5(z1,z2,k) = CVT % VertMode5(z1,z2,k) + &
                                     pert % v5(real_index,z1) * &
                                     pert % v5(real_index,z2) + &
                                     pert % v5(imag_index,z1) * &
                                     pert % v5(imag_index,z2)
          CVT % VertMode6(z1,z2,k) = CVT % VertMode6(z1,z2,k) + &
                                     pert % v6(real_index,z1) * &
                                     pert % v6(real_index,z2) + &
                                     pert % v6(imag_index,z1) * &
                                     pert % v6(imag_index,z2)
        END DO
      END DO
    END DO

    ! The k=nlongs/2+1 component though is real
    k = nlongs/2 + 1
    DO z1 = 1, nlevs
      DO z2 = 1, nlevs
        CVT % VertMode1(z1,z2,k) = CVT % VertMode1(z1,z2,k) + &
                                   pert % v1(k, z1) *         &
                                   pert % v1(k, z2)
        CVT % VertMode2(z1,z2,k) = CVT % VertMode2(z1,z2,k) + &
                                   pert % v2(k, z1) *         &
                                   pert % v2(k, z2)
        CVT % VertMode3(z1,z2,k) = CVT % VertMode3(z1,z2,k) + &
                                   pert % v3(k, z1) *         &
                                   pert % v3(k, z2)
        CVT % VertMode4(z1,z2,k) = CVT % VertMode4(z1,z2,k) + &
                                   pert % v4(k, z1) *         &
                                   pert % v4(k, z2)
        CVT % VertMode5(z1,z2,k) = CVT % VertMode5(z1,z2,k) + &
                                   pert % v5(k, z1) *         &
                                   pert % v5(k, z2)
        CVT % VertMode6(z1,z2,k) = CVT % VertMode6(z1,z2,k) + &
                                   pert % v6(k, z1) *         &
                                   pert % v6(k, z2)
      END DO
    END DO


  END SELECT

ELSE
  ! Normalize
  mul_cons = 1.0 / div_cons

  SELECT CASE (CVT % CVT_order)

  CASE (1)
    ! Original MetO order (horiz co-ordinate of pert is longitude)
    lim = 1
  CASE (2)
    ! Reverse MetO order (hoiz co-ordinate of pert is wavenumber)
    lim = nlongs/2 + 1
  END SELECT

  DO k = 1, lim
    CVT % VertMode1(1:nlevs,1:nlevs,k) = CVT % VertMode1(1:nlevs,1:nlevs,k) * mul_cons
    CVT % VertMode2(1:nlevs,1:nlevs,k) = CVT % VertMode2(1:nlevs,1:nlevs,k) * mul_cons
    CVT % VertMode3(1:nlevs,1:nlevs,k) = CVT % VertMode3(1:nlevs,1:nlevs,k) * mul_cons
    CVT % VertMode4(1:nlevs,1:nlevs,k) = CVT % VertMode4(1:nlevs,1:nlevs,k) * mul_cons
    CVT % VertMode5(1:nlevs,1:nlevs,k) = CVT % VertMode5(1:nlevs,1:nlevs,k) * mul_cons
    CVT % VertMode6(1:nlevs,1:nlevs,k) = CVT % VertMode6(1:nlevs,1:nlevs,k) * mul_cons
  END DO
END IF

END SUBROUTINE CVT_Calibration_vertcovs
