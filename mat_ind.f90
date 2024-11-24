!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mat_ind(nl,numstr,totstr,strno,Ifail,det_inv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use commondat
use commondat1
implicit none

!common/str/str5
common/fail/faiil
!common/sdt/strdet,detmnt,det_sign
integer::nl,str1(10000,15),Ifail,totstr,tndet,vec(10000,1000),numstr,strno(1000),detmnt1(10000,15)&
,i,j,k,l,m,n,i1,comdet(10000),mat(10000,10000),sdet(10000),mat_sign(10000,10000),S00(10000,10000)&
,det_sign1(10000),faiil,det_inv
!,strdet(100000),detmnt(100000,15),det_sign(100000)
integer::col(1000)
DOUBLE PRECISION :: SPL0,SPL,D1(100)
character(len=3)::ind

!write(9,*),'enter mat_ind'
!do l=1,totstr*ndet
!print*,'strdet(l)',l,'<',strdet(l)
!enddo
!totstr=6
!call vector_rep(nl,str1,totstr,vec,detmnt,strdet,det_sign,tndet)
tndet=0
!write(9,*),'totstr,numstr',totstr,numstr
!write(9,103),'mat_ind_strno',ndet*totstr,(strno(i),i=1,numstr)

!print*,'numstr',numstr,totstr,ndet
do i=1,numstr
j=strno(i)
!write(9,*),'structures',i,j,numstr
do l=1,totstr*ndet
!write(9,*),'lll',l
if(strdet(l).ne.j)goto 200
!write(9,*),'strdet(l)',strdet(l),j,totstr*ndet
tndet=tndet+1
do m=1,nae
detmnt1(tndet,m)=detmnt(l,m)
enddo
!write(9,103),'detmnt(l,m)',(detmnt(l,m),m=1,nae),det_sign(l)
det_sign1(tndet)=det_sign(l)
200 enddo
enddo
!do i=1,tndet
!write(9,103),'detmnt1***',(detmnt1(i,m),m=1,nae)
!enddo
!print*,'totstr*ndet,tndet',totstr,ndet,tndet

l=0
n=0
do i=1,tndet
!print*,'det',i,(detmnt(i,j),j=1,nae)
!print*,'lll',l
do i1=1,l
if(comdet(i1).eq.i)goto 101
enddo
n=n+1
m=0
do j=i,tndet
do k=1,nae
if(detmnt1(i,k).ne.detmnt1(j,k))goto 100
enddo
l=l+1
m=m+1
mat(n,m)=strdet(j)
mat_sign(n,m)=det_sign1(j)
comdet(l)=j
!print*,'iii',l,comdet(l),n,m,mat(n,m),mat_sign(n,m)
100 enddo
sdet(n)=m
!print*,'sdet(n)',n,sdet(n)
101 enddo

!do i=1,n
!print*,'**',i,sdet(i)
!write(*,102),(mat(i,j),j=1,sdet(i))
!enddo
!do i=1,n
!write(*,102),(mat_sign(i,j),j=1,sdet(i))
!enddo
102 format(50I5)
103 format(a,2x,50I5)

do i=1,totstr
do j=1,totstr
ind_mat(i,j)=0.0
enddo
enddo

do i=1,totstr
ind_mat(i,i)=ndet
!write(9,*),'mat_ind::ndet',ndet
enddo

do i=1,numstr
!print*,i,'****',i
do j=i+1,numstr
!print*,'j',j
do k=1,n
!print*,'k',k
do l=1,sdet(k)
!print*,'l,sdet(k),mat(k,l),i',l,sdet(k),mat(k,l),i
if(mat(k,l).eq.i)then
do m=1,sdet(k)
!print*,'m,sdet(k),mat(k,m),j',m,sdet(k),mat(k,m),j
if(mat(k,m).eq.j)then
ind_mat(i,j)=ind_mat(i,j)+mat_sign(k,m)*mat_sign(k,l)
!print*,'i,j,k,l,m,ind_mat(i,j)',i,j,k,l,m,ind_mat(i,j),mat_sign(k,m),mat_sign(k,l)
endif
enddo
endif
enddo
enddo
enddo
enddo
!      Do k=1,n
!print*,'sdet(k)',sdet(k)
!       Do i=1,sdet(k)
!        IIs=Nsd(Is,Jdet)
!        Spl0=Dble(mat_sign(k,i))
!print*,'Spl0',k,i,Spl0
!        Do j=i,sdet(k)
!         JJs=Nsd(Js,Jdet)
!         Spl=Spl0*Dble(mat_sign(k,j))
!         S00(j,i)=S00(j,i)+Spl
!        Enddo
!       Enddo
!      Enddo
      Do i=1,numstr
       Do j=i+1,numstr
        ind_mat(j,i)=ind_mat(i,j)
       Enddo
      Enddo

!Do i=1,numstr
!write(*,102),(ind_mat(i,j),j=1,numstr)
!
!!write(*,102),(S00(i,j),j=1,totstr)
!Enddo


!call MatLDR('ind',col,numstr,D1)
!det_inv=1.0
!do i=1,numstr
!if (dabs(D1(i)).lt.1.D-10)then
!Ifail=1
!goto 111
!endif
!
!det_inv=det_inv*D1(i)
!enddo
!print*,'mat_ind',(D1(i),i=1,numstr),'det_inv=',det_inv
!
!if (dabs(det_inv).lt.1.D-10)then
!Ifail=1
!else
!Ifail=0
!endif

call Invmat(numstr,Ifail)
111 faiil=Ifail
!if(numstr.eq.14.and.Ifail.eq.0)then
!print*,(strno(i),i=1,numstr)
!if(Ifail.eq.0)then
!write(9,*),'Ifail:invmat',Ifail
!do i=1,numstr
!write(9,103),'mat_ind_str',(str5(strno(i),j),j=1,nae)
!enddo
!Do i=1,numstr
!write(9,102),(ind_mat(i,j),j=1,numstr)
!!!!!
!!!!!!write(*,102),(S00(i,j),j=1,totstr)
!Enddo
!endif
!endif
!if(Ifail.eq.0.and.numstr.eq.5)serial=serial+1
!print*,'serial=',serial
!print*,'exit mat_ind'

!if(numstr.eq.7)stop
return
end subroutine mat_ind
