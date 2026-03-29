!#**********************************************************************************
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! miscellaneous debuging routines for CBMZ and MOSAIC - ADAPTED FOR MPAS
!**********************************************************************************
module module_peg_util

  implicit none

  private

  ! Public interfaces to be used by other modules
  public :: peg_debugmsg, &
            peg_message, &
            peg_error_fatal

contains

!-----------------------------------------------------------------------
subroutine peg_debugmsg( lun, level, str )
!
! Writes debug message to standard output
!
  implicit none
  ! subr arguments
  integer, intent(in) :: lun, level
  character(len=*), intent(in) :: str
  ! local variables
  integer n

  n = max( 1, len_trim(str) )
  write(*,'(a)') str(1:n)

  return
end subroutine peg_debugmsg


!-----------------------------------------------------------------------
subroutine peg_message( lun, str )
!
! Writes message to standard output
!
  implicit none
  ! subr arguments
  integer, intent(in) :: lun
  character(len=*), intent(in) :: str
  ! local variables
  integer n

  n = max( 1, len_trim(str) )
  write(*,'(a)') str(1:n)

  return
end subroutine peg_message


!-----------------------------------------------------------------------
subroutine peg_error_fatal( lun, str )
!
! Writes fatal error message and stops execution
!
  implicit none
  ! subr arguments
  integer, intent(in) :: lun
  character(len=*), intent(in) :: str
  ! local variables
  integer n

  n = max( 1, len_trim(str) )

  write(*,'(a)') ' '
  write(*,'(a)') '------------- FATAL ERROR in Mie code ------------'
  write(*,'(a)') str(1:n)
  write(*,'(a)') 'MPAS Mie Code abort due to fatal error.'
  
  stop 'FATAL ERROR'
  
  return
end subroutine peg_error_fatal


!-----------------------------------------------------------------------
end module module_peg_util
