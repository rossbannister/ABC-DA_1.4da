SUBROUTINE VertEigens (CVT)

! Compute the vertical eigenvectors and eigenvalues for each control variable

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    CVT_type,                &
    nlongs, nlevs,           &
    small, ConditionFudge

IMPLICIT NONE

TYPE(CVT_type), INTENT(INOUT)   :: CVT

INTEGER                         :: k, lim, ierr
REAL(ZREAL8)                    :: Work(1:3*nlevs), increment

  SELECT CASE (CVT % CVT_order)

  CASE (1)
    ! Original MetO order (horiz co-ordinate of pert is longitude)
    lim = 1
  CASE (2)
    ! Reverse MetO order (hoiz co-ordinate of pert is wavenumber)
    lim = nlongs/2 + 1
  END SELECT

  DO k = 1, lim
    IF (CVT % CVT_order == 2) PRINT *, '    Wavenumber', k, ' out of', nlongs/2+1
    CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
               'U',                                  & ! Upper triangular of matrix specified
               nlevs,                                & ! Order of the matrix to be diagonalised
               CVT % VertMode1(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
               nlevs,                                & ! Leading dimension of the matrix
               CVT % VertEV1(1:nlevs,k),             & ! Eigenvalue array
               Work(1:3*nlevs),                      & ! Work array
               3*nlevs,                              & ! Size of work array
               ierr)
    ! Condition for zero values
    increment = SUM(CVT % VertEV1(1:nlevs,k)) * ConditionFudge
    IF (increment < small) increment = small
    CVT % VertEV1(1:nlevs,k) = CVT % VertEV1(1:nlevs,k) + increment
    CVT % VertEV1(1:nlevs,k) = SQRT(CVT % VertEV1(1:nlevs,k))

    CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
               'U',                                  & ! Upper triangular of matrix specified
               nlevs,                                & ! Order of the matrix to be diagonalised
               CVT % VertMode2(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
               nlevs,                                & ! Leading dimension of the matrix
               CVT % VertEV2(1:nlevs,k),             & ! Eigenvalue array
               Work(1:3*nlevs),                      & ! Work array
               3*nlevs,                              & ! Size of work array
               ierr)
    ! Condition for zero values
    increment = SUM(CVT % VertEV2(1:nlevs,k)) * ConditionFudge
    IF (increment < small) increment = small
    CVT % VertEV2(1:nlevs,k) = CVT % VertEV2(1:nlevs,k) + increment
    ! Square-root
    CVT % VertEV2(1:nlevs,k) = SQRT(CVT % VertEV2(1:nlevs,k))

    CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
               'U',                                  & ! Upper triangular of matrix specified
               nlevs,                                & ! Order of the matrix to be diagonalised
               CVT % VertMode3(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
               nlevs,                                & ! Leading dimension of the matrix
               CVT % VertEV3(1:nlevs,k),             & ! Eigenvalue array
               Work(1:3*nlevs),                      & ! Work array
               3*nlevs,                              & ! Size of work array
               ierr)
    ! Condition for zero values
    increment = SUM(CVT % VertEV3(1:nlevs,k)) * ConditionFudge
    IF (increment < small) increment = small
    CVT % VertEV3(1:nlevs,k) = CVT % VertEV3(1:nlevs,k) + increment
    ! Square-root
    CVT % VertEV3(1:nlevs,k) = SQRT(CVT % VertEV3(1:nlevs,k))

    CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
               'U',                                  & ! Upper triangular of matrix specified
               nlevs,                                & ! Order of the matrix to be diagonalised
               CVT % VertMode4(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
               nlevs,                                & ! Leading dimension of the matrix
               CVT % VertEV4(1:nlevs,k),             & ! Eigenvalue array
               Work(1:3*nlevs),                      & ! Work array
               3*nlevs,                              & ! Size of work array
               ierr)
    ! Condition for zero values
    increment = SUM(CVT % VertEV4(1:nlevs,k)) * ConditionFudge
    IF (increment < small) increment = small
    CVT % VertEV4(1:nlevs,k) = CVT % VertEV4(1:nlevs,k) + increment
    ! Square-root
    CVT % VertEV4(1:nlevs,k) = SQRT(CVT % VertEV4(1:nlevs,k))

    CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
               'U',                                  & ! Upper triangular of matrix specified
               nlevs,                                & ! Order of the matrix to be diagonalised
               CVT % VertMode5(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
               nlevs,                                & ! Leading dimension of the matrix
               CVT % VertEV5(1:nlevs,k),             & ! Eigenvalue array
               Work(1:3*nlevs),                      & ! Work array
               3*nlevs,                              & ! Size of work array
               ierr)
    ! Condition for zero values
    increment = SUM(CVT % VertEV5(1:nlevs,k)) * ConditionFudge
    IF (increment < small) increment = small
    CVT % VertEV5(1:nlevs,k) = CVT % VertEV5(1:nlevs,k) + increment
    ! Square-root
    CVT % VertEV5(1:nlevs,k) = SQRT(CVT % VertEV5(1:nlevs,k))

    CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
               'U',                                  & ! Upper triangular of matrix specified
               nlevs,                                & ! Order of the matrix to be diagonalised
               CVT % VertMode6(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
               nlevs,                                & ! Leading dimension of the matrix
               CVT % VertEV6(1:nlevs,k),             & ! Eigenvalue array
               Work(1:3*nlevs),                      & ! Work array
               3*nlevs,                              & ! Size of work array
               ierr)
    ! Condition for zero values
    increment = SUM(CVT % VertEV6(1:nlevs,k)) * ConditionFudge
    IF (increment < small) increment = small
    CVT % VertEV6(1:nlevs,k) = CVT % VertEV6(1:nlevs,k) + increment
    ! Square-root
    CVT % VertEV6(1:nlevs,k) = SQRT(CVT % VertEV6(1:nlevs,k))

  END DO

END SUBROUTINE VertEigens
