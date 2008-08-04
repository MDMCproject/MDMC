module converters_class
use various_constants_class, only : db

implicit none

private
!
! Takes a string and turns it into useful data types
!
!
! NOTE: the code below is a modified version of the m_converters.f90
!       in the xmf90 sax library.
!
! There might be better ways to convert a string an data types but
! this brute force was the easiest and fastest way I could think of
! doing this
!
public :: string_to_int
public :: string_to_db 


private :: token_analysis, is_separator, is_CR_or_LF

CONTAINS

!---------------------------------------------------------------
!subroutine build_data_array_real_dp(str,x,n)
!integer, parameter  :: dp = selected_real_kind(14)
!
!character(len=*), intent(in)                ::  str
!real(kind=dp), dimension(:), intent(inout)  ::    x
!integer, intent(inout)                      ::    n

!integer                            :: ntokens, status, last_pos
!character(len=len(str))  :: s

!s = str
!call token_analysis(s,ntokens,last_pos)
!if (debug) print *, "ntokens, last_pos ", ntokens, last_pos
!if (debug) print *, s
!if ((n + ntokens) > size(x)) STOP "data array full"
!read(unit=s(1:last_pos),fmt=*,iostat=status) x(n+1:n+ntokens)
!if (status /= 0) STOP "real conversion error"
!n = n + ntokens

!end subroutine build_data_array_real_dp
!---------------------------------------------------------------

function string_to_db(str) result(db_out)
  character(len=*), intent(in) :: str

  real(db) :: db_out

  integer                            :: ntokens, status, last_pos
  character(len=len(str))  :: s

  s = str
  call token_analysis(s,ntokens,last_pos)
  if (ntokens /= 1) STOP "ERROR in string_to_db"
  read(unit=s(1:last_pos),fmt=*,iostat=status) db_out
  if (status /= 0) then
    print *, "Cannot convert the following string to double:"
    print *, str
    STOP "db conversion error"
  end if

end function string_to_db


function string_to_int(str) result(int_out)
  character(len=*), intent(in) :: str

  integer :: int_out

  integer                            :: ntokens, status, last_pos
  character(len=len(str))  :: s

  s = str
  call token_analysis(s,ntokens,last_pos)
  if (ntokens /= 1) STOP "ERROR in string_to_int"
  read(unit=s(1:last_pos),fmt=*,iostat=status) int_out
  if (status /= 0) STOP "integer conversion error"

end function string_to_int


!==================================================================

function is_separator(c) result(sep)
character(len=1), intent(in)          :: c
logical                               :: sep

 sep = ((c == char(32)) .or. (c == char(10))             &
         .or. (c == char(9)) .or. (c == char(13)))

end function is_separator
!----------------------------------------------------------------
function is_CR_or_LF(c) result(res)
character(len=1), intent(in)          :: c
logical                               :: res

 res = ((c == char(10)) .or. (c == char(13)))

end function is_CR_or_LF

!==================================================================

subroutine token_analysis(str,ntokens,last_pos)
!
character(len=*), intent(inout)          :: str
integer, intent(out)                     :: ntokens, last_pos
!
!
! Checks the contents of a string and finds the number of tokens it contains
! The standard separator is generalized whitespace (space, tab, CR, or LF)
! It also returns the last useful position in the string (excluding
! separator characters which are not blanks, and thus not caught by the
! (len_)trim fortran intrinsic). This is necessary to perform list-directed
! I/O in the string as an internal file.
! 
! Also, replace on the fly CR and LF by blanks. This is necessary if
! str spans more than one record. In that case, internal reads only 
! look at the first record. 
! -- ** Compiler limits on size of internal record??
!
integer           :: i, str_length
logical           :: in_token
character(len=1)  :: c

in_token = .false.
ntokens = 0
last_pos = 0

str_length = len_trim(str)
!print *, "string length: ", str_length

do i = 1, str_length
      c = str(i:i)

      if (in_token) then
         if (is_separator(c)) then
            in_token = .false.
            if (is_CR_or_LF(c)) str(i:i) = " "
         else
            last_pos = i
         endif

      else   ! not in token
         
         if (is_separator(c)) then
            if (is_CR_or_LF(c)) str(i:i) = " "
            ! do nothing
         else
            in_token = .true.
            last_pos = i
            ntokens = ntokens + 1
         endif
      endif
enddo
!print *, "ntokens, last_pos: ", ntokens, last_pos

end subroutine token_analysis


end module converters_class







