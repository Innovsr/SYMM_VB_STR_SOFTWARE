!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine get the arguments from the python input GUI interface!
!!! and create the required flags and update the values in commondat for !
!!! global use                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cisvb_inp (ctrl_keys, values, n)
!use commondat
!use commondat1

implicit none

integer:: i, j, k, l, m, n, values(n)
CHARACTER(LEN=*), DIMENSION(n) :: ctrl_keys

do i=1,n
print*,'ctrl_key:',ctrl_keys(i),' = ',values(i) 
enddo

stop


end subroutine cisvb_inp
