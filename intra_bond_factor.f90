!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intra_bond_factor(nl,str,tonstruc,str_quality_1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::nl,i2,m8,l1,l2,l3,k1,k2,m13,m14,str(15000,20),str_quality_1(15000),tonstruc

print*,'enter intra_bond_factor'
do m8=1,15000
str_quality_1(m8)=1
enddo

l3=1
do m8=1,tonstruc
!print*,'ssttr',(str(m8,i2),i2=1,nae)
l2=1
do k2=1+nl*2,nae-nlast,2
do m13=1,atom
l1=0
do m14=1,atn(active_atoms(m13))
if(atn(active_atoms(m13)).eq.1)goto 505
do k1=k2,k2+1
if(str(m8,k1).eq.atoset(active_atoms(m13),m14))then
l1=l1+1
endif
enddo
505 enddo
if(l1.eq.2) then
l2=l2+l3
goto 506
endif
enddo
506 enddo
str_quality_1(m8)=l2
!print*,'intra_bond_factor',str_quality_1(m8)
enddo


print*,'exit intra_bond_factor'
return
end subroutine intra_bond_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
