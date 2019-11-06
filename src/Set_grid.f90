SUBROUTINE Set_grid (init_um_data, dims, regular_grid)

!**********************************
!* Subroutine to set model grid   *
!*                                *
!*                                *
!* Set regular_grid = 1 to use a  *
!* regular vertical grid          *
!*                                *
!* R. Petrie                      *
!* version 1.0                    *
!* 6/6/2011                       *
!**********************************

USE DefConsTypes, ONLY :   &
    ZREAL8,                &
    UM_type,               &
    dims_type,             &
    nlongs, nlevs,         &
    dx

IMPLICIT NONE

! Declare Top level variables
TYPE(UM_type),   INTENT(IN)     :: init_um_data
TYPE(dims_type), INTENT(INOUT)  :: dims
INTEGER,         INTENT(IN)     :: regular_grid

! Declare local variables
REAL(ZREAL8)                    :: level_reg_fact, half_lev_fact, halfdx
INTEGER                         :: i

! Set vertical dimensions
IF (regular_grid .EQ. 1) THEN
  PRINT*, 'Using regular vertical grid'
  level_reg_fact = init_um_data % full_levs(nlevs) / nlevs
  half_lev_fact  = level_reg_fact * 0.5
  halfdx         = dx * 0.5
  DO i = 0, nlevs
     dims % full_levs(i) = REAL(i) * level_reg_fact
  ENDDO

  DO i = 1, nlevs+1
     dims % half_levs(i) = half_lev_fact + (REAL(i-1) * level_reg_fact)
  ENDDO
  dims % half_levs(0) = -1. * dims % half_levs(1)
ELSE
  PRINT*, 'Using MetO model vertical grid'
  dims % full_levs(0:nlevs)   = init_um_data % full_levs(0:nlevs)
  dims % half_levs(1:nlevs+1) = init_um_data % half_levs(1:nlevs+1)
  dims % half_levs(0)         = -1. * dims % half_levs(1)
END IF

! Set horizontal dimension
DO i = 0,nlongs+1
  dims % longs_v(i) = REAL(i) * dx
ENDDO

DO i = 0,nlongs+1
  dims % longs_u(i) = REAL(i) * dx + halfdx
ENDDO

END SUBROUTINE Set_grid
