!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nnat_bond_cal(nl,str1,ncqs,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! This subroutine calculates the number of nearest neighbour bonds !!!!!

use commondat
implicit none

integer::m19,m18,i1,i2,i3,i4,i5,i6,i7,ii,iii,iiii,nl,ncqs
integer::nn(10),str1(15000,20),bondq(15000)

print*,'enter nnat_bond_cal'
do i1=1,15000
bondq(i1)=0
enddo

if(nnnatom.eq.0) goto 900
231 format (30I3)
i4=0
do i1=1,ncqs
write(*,231)(str1(i1,i2),i2=1,nae)
iiii=1+(nae-nl*2-nlast)/2
do i3=1+nl*2,(nae-nlast),2
iii=0
nn(1)=0
nn(2)=0
do i5=i3,i3+1
do i4=1,atom
do i7=1,atn(active_atoms(i4))
if(str1(i1,i5).eq.atoset(active_atoms(i4),i7))then
iii=iii+1
if(mod(i5,2).eq.0) nn(2)=i4
if(mod(i5,2).eq.1) nn(1)=i4
goto 517
endif
enddo

enddo
517 enddo
if(nn(2).eq.nn(1))goto 400
if(iii.ne.2) goto 400
do i6=1,nnnatom
ii=0
do i7=1,2
do i4=1,2
if(nn(i4).eq.nnat_bond(i6,i7))then
ii=ii+1
endif
enddo
enddo
if(ii.eq.2)then
iiii=iiii-1
goto 400
endif
enddo

400 enddo
bondq(i1)=iiii
!print*,'nnbd:iiii',iiii
enddo

print*,'exit nnat_bond_cal'
900 return
end subroutine nnat_bond_cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
