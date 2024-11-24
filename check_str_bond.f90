!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! checking the bonds of the structures how many are presented there which!
! helps to identify most effective equally distributed sets.             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_str_bond(nl,strn,str3,ncqs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,j,j1,j2,nl,strn,ncqs,str3(15000,20),num_orb(20),num_orb1(20),&
totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),least
common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb


print*,'enter check_str_bond'

i1=0
do i=1,nae-nl*2
  i1=i1+1
  num_orb(i1)=str3(1,i+nl*2)
enddo

! searchinh lowest orbital number
least=1000
do i2=1,i1
  if(least.gt.num_orb(i2))least=num_orb(i2)
enddo

i3 = 0
do
    do i2 = 1, i1
        if (least == num_orb(i2)) then
            i3 = i3 + 1
            num_orb1(i3) = least
            exit  ! Exit the inner loop once the condition is met
        endif
    end do
    least = least + 1
    if (i3 >= i1) exit  ! Exit the outer loop when the condition is satisfied
end do


!i3=0
!320 do i2=1,i1
!  if(least.eq.num_orb(i2))then
!    i3=i3+1
!    num_orb1(i3)=least
!    goto 321
!  endif
!enddo
!321 least=least+1
!if (i3.lt.i1)goto 320

do i=1,i1
num_orb(i)=num_orb1(i)
enddo

i2=0
do j1=1,i1
do j2=j1+1,i1
i2=i2+1
bond_sl(i2,1)=num_orb(j1)
bond_sl(i2,2)=num_orb(j2)
enddo
enddo

totbnd=i2
totseorb=i1

do i=1,15000
do i1=1,1000
str_bnd(i,i1)=0
enddo
enddo

do i=1,ncqs
do i1=nl*2+1,nae-nlast,2
do j=1,totbnd
j1=0
do i2=i1,i1+1
do i3=1,2
if(str3(i,i2).eq.bond_sl(j,i3))then
j1=j1+1
endif
enddo
enddo
if(j1.eq.2)then
str_bnd(i,j)=1
goto 200
endif
enddo
200 enddo
enddo

do i=1,500
do i1=1,20
str_rad(i,i2)=0
enddo
enddo

if(nlast.ne.0)then
do i=1,ncqs
do i1=nae-nlast+1,nae
do i2=1,totseorb
if(str3(i,i1).eq.num_orb(i2))then
str_rad(i,i2)=1
goto 201
endif
enddo
201 enddo
enddo
endif

print*,'exit check_str_bond'
return
end subroutine check_str_bond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
