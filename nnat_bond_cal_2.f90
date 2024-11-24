!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nnat_bond_cal_2(nl,str1,ncqs,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::m19,m18,i,i1,i2,i3,i4,i5,i6,i7,iii,iiii,nl,ncqs,k
integer::nn(10),str1(15000,20),bondq(15000),sl(5000)
real::least,bondq1_dist(15000)
double precision::ii,bondq_dist(15000)

print*,'enter nnat_bond_cal_2'

open(unit=13,file='nnbd.temp',status='unknown')
do i1=1,15000
bondq(i1)=0
bondq_dist(i1)=0.0
enddo

231 format (30I3)
i4=0
do i1=1,ncqs
ii=0.0
!write(*,231)i1,(str1(i1,i2),i2=1,nae)
do i3=1+nl*2,(nae-nlast),2
iii=0
nn(1)=0
nn(2)=0
do i5=i3,i3+1
do i4=1,atom
do i7=1,atn(active_atoms(i4))
!print*,'atn(i4)',atn(active_atoms(i4)),active_atoms(i4)
!print*,'atoset',str1(i1,i5),atoset(active_atoms(i4),i7)
if(str1(i1,i5).eq.atoset(active_atoms(i4),i7))then
iii=iii+1
if(mod(i5,2).eq.0) nn(2)=active_atoms(i4)
if(mod(i5,2).eq.1) nn(1)=active_atoms(i4)
!print*,'------------',nn(1),nn(2)
goto 517
endif
enddo

enddo
517 enddo
print*,'nn(2),nn(1)',nn(2),nn(1),iii
if(iii.ne.2) goto 400
if(nn(2).eq.nn(1))then

!!! changes done for EDEN
iab_length=1.00
!!!!!!!!!!!!!!!!!!!!!!!!!

if (iab_length.eq.0.0)ii=ii+100.0
if (iab_length.ne.0.0)ii=ii+iab_length
!print*,'iab_length',iab_length
goto 400
endif
ii=ii+dist_act_rel_mat(nn(1),nn(2))
!print*,'dist_act_rel_mat(nn(1),nn(2))',dist_act_rel_mat(nn(1),nn(2)),ii
400 enddo
bondq_dist(i1)=ii
!print*,'bondq_dist(i1)',bondq_dist(i1),ii
write(13,*)bondq_dist(i1)
enddo
rewind(13)
do i1=1,ncqs
read(13,105)bondq1_dist(i1)
!write(*,*),i1,bondq1_dist(i1)
enddo
105 format (F8.4)
iii=0
iiii=0
211 least=10000.0
do i=1,ncqs
do i1=1,iii
if(i.eq.sl(i1))goto 210
enddo
if(bondq1_dist(i).le.least)least=bondq1_dist(i)
210 enddo
!print*,'least',least
if(least.ne.100.0.or.least.ne.0.0)iiii=iiii+1
do i=1,ncqs
if(bondq1_dist(i).eq.least)then
iii=iii+1
bondq(i)=iiii
sl(iii)=i
endif
enddo
if(iii.ne.ncqs)goto 211

!do i=1,ncqs
!write(*,231)i,(str1(i,i2),i2=1,nae),bondq(i)
!!print*,'bondq(i)',i,bondq(i)
!enddo
CALL SYSTEM ("rm nnbd.temp")
print*,'exit nnat_bond_cal_2'
900 return
end subroutine nnat_bond_cal_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
