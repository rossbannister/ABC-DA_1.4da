SUBROUTINE Set_ht_dep_cons (dims)

!********************************
!* Subroutine to set some height*
!* dependent constants          *
!* R. Bannister                 *
!* version 1.0                  *
!* 20/02/14                     *
!********************************

USE DefConsTypes, ONLY :   &
    dims_type,       &
    nlevs

IMPLICIT NONE

! Declare Top level variables
TYPE(dims_type), INTENT(INOUT)  :: dims

! Declare local variables
INTEGER                               :: z

! For full-to-half level interpolation
DO z = 1, nlevs+1
   dims % a1(z) = (dims % half_levs(z) - dims % full_levs(z-1)) / (dims % full_levs(z) - dims % full_levs(z-1))
   dims % b1(z) = (dims % full_levs(z) - dims % half_levs(z))   / (dims % full_levs(z) - dims % full_levs(z-1))
ENDDO

! For half-to-full level interpolation
DO z = 0, nlevs
   dims % a2(z) = (dims % full_levs(z)   - dims % half_levs(z)) / (dims % half_levs(z+1) - dims % half_levs(z))
   dims % b2(z) = (dims % half_levs(z+1) - dims % full_levs(z)) / (dims % half_levs(z+1) - dims % half_levs(z))
ENDDO

! Reciprocal of some level differences used commonly
DO z = 1, nlevs+1
  dims % recip_half_k_km1(z) = 1.0 / ( dims % half_levs(z) - dims % half_levs(z-1) )
  dims % recip_full_k_km1(z) = 1.0 / ( dims % full_levs(z) - dims % full_levs(z-1) )
END DO
DO z = 0, nlevs
  dims % recip_half_kp1_k(z) = 1.0 / ( dims % half_levs(z+1) - dims % half_levs(z) )
  dims % recip_full_kp1_k(z) = 1.0 / ( dims % full_levs(z+1) - dims % full_levs(z) )
END DO


END SUBROUTINE Set_ht_dep_cons

