!    The module prints the Date,Month,Year,Hours,Minutes and seonds.
!
!    Programmer :
!                 Shahid Maqbool
!
!    Modified :
!                January 20, 2025

module timestamp_module
  use precision_module
  implicit none
  
  character (20) ampm
  integer (sp) day
  integer (sp) hour
  integer (sp) minute
  integer (sp) mm
  character (len=9), parameter,dimension (12) :: month = (/ &
  'January  ', 'February ', 'March    ', 'April    ', &
  'May      ', 'June     ', 'July     ', 'August   ', &
  'September', 'October  ', 'November ', 'December ' /)
  integer (sp) n
  integer (sp) sec
  integer (sp) values (8)
  integer (sp) year
  
  contains
  
  subroutine time_stamp()
    
    call date_and_time (values=values)
    
    year = values (1)
    minute = values (2)
    day = values (3)
    hour = values (5)
    n = values (6)
    sec = values (7)
    mm = values (8)
    
    if ( hour<12 ) then
    ampm = 'AM'
    else if (hour==12) then
      if (n==0 .and. sec==0) then
      ampm = 'Noon'
      else
        ampm = 'PM'
      end if
    
    else
      hour = hour - 12
      if ( hour<12 ) then
      ampm = 'PM'
      else if (hour==12) then
        if ( n==0 .and. sec==0 ) then
        ampm = 'Midnight'
        else
          ampm = 'AM'
        end if
      end if
    end if
    
    write (*,'(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    day,trim(month(minute)),year,hour,':',n,':',sec,'.',mm,trim(ampm)
    
    return
    
  end subroutine time_stamp

end module timestamp_module