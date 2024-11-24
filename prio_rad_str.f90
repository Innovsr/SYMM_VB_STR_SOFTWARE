!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prio_rad_str(nl,str1,ncqs,pref_radical)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,j,i3,i4,i5,i6,i7,ii,m16,m17,m18,m19,m23,nl,ncqs,qual2(15000),qul1(15000),&
qul2(15000),jj,nqul1,nqul2,qt,qual3(15000),qual4(15000),k6,k7,k8,pref_radical(15000)
integer::str1(15000,20),qual1(15000),qcs_serial(15000),q_fac(15000)

!!!!! starting the arrangement of the structures according to priority radicles!!!


print*,'enter prio_rad_str'
do i=1,15000
!pref_radical(i)=2
pref_radical(i)=nlast+1
enddo

!print*,'nlpset',nlpset
if(nlpset.eq.0)goto 102
do i5=1,nlpset
ii=0
do i=1,nl*2,2
!print*,'i',i
do j=1,nl
if(plpair(i5,j).eq.str1(1,i))then
ii=ii+1
endif
enddo
enddo
!print*,'ii',ii
if(ii.ne.nl) goto 103
if(ii.eq.nl)then
jj=i5
if(plpair(i5,nl+1).ne.0.)goto 100
if(plpair(i5,nl+1).eq.0.)goto 102
endif
103 enddo
goto 101
!102 print*,'i51',jj,plpair(i5,nl+1)


!102 print*,'sourav'
102 do i=1,ncqs
!print*,'sourav1',(str1(i,i4),i4=1,nae)
do j=1,prad
!print*,'norad',norad(j)
ii=0
do i3=1,norad(j)
do i4=nae-nlast+1,nae
if(str1(i,i4).eq.prio_rad(j,i3))then
ii=ii+1
endif
enddo
enddo
if(ii.eq.norad(j))then
pref_radical(i)=pref_radical(i)-1
!print*,'pref_radical(i)',pref_radical(i)
endif
enddo
enddo

!print*,'sourav2'
if(nlpset.eq.0)goto 101
100 if(plpair(jj,nl+1).eq.0) goto 101


do i=1,ncqs
do j=nl+1,lp(jj)
i7=plpair(jj,j)
ii=0
do i3=1,norad(i7)
do i4=nae-nlast+1,nae
if(str1(i,i4).eq.prio_rad(i7,i3))then
ii=ii+1
endif
enddo
enddo
!print*,'norad(i7)',norad(i7)
if(ii.eq.norad(i7))then
pref_radical(i)=pref_radical(i)-1
endif
enddo
enddo

!101 enddo


!101 do m19=1,ncqs
!!print*,'pref_rad',pref_radical(m19)
!write(*,231),m19,(str1(m19,i4),i4=1,nae)
!enddo
231 format(30I3)

print*,'exit prio_rad_str'

101 return
end subroutine prio_rad_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
