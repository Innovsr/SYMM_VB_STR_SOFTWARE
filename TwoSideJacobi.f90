!----------------------------------------------------------------
      SUBROUTINE TwoSideJacobi(i,j,Aii,Ajj,Aij,Aji,U,VT,A,N)
!----------------------------------------------------------------
!
!     This SUBROUTINE calculate two sides Jacobi rotation to following
!     matrix
!     | Aii Aij |       | Aii~  0  |
!     |         |  -->  |          |
!     | Aji Ajj |       |  0  Ajj~ |
!
      Implicit none
      Double Precision A(N,N),U(N,N),D(N),VT(N,N)
      Double Precision t1,c1,s1
      Double Precision t2,c2,s2
      Double Precision tp,cp,sp
      Double Precision tm,cm,sm
      Double Precision Aii,Ajj,Aij,Aji,tt,Aod,Adif,tmp
      integer N
      integer I,J,k
      Double Precision Zer,One,Eps
      Data Zer,One,Eps/0.D0,1.D0,1.0D-20/
!     write(0,*)"Size 2 = ",Size(VT,1),size(vt,2),size(vt)
!print*,'A11,A12,A21,A22',A(1,1),A(1,2),A(2,1),A(2,2)

!     angle alpha + beta
!     write(0,*)'I J = ',I,J
!     write(0,*)'aii = ',aii
!     write(0,*)'ajj = ',ajj
!     write(0,*)'aij = ',aij
!     write(0,*)'aji = ',aji
      Adif = Aii - Ajj
      Aod = Aij + Aji
      if(dabs(Adif) .lt. dabs(Aod)) then
!       tt is cot(2x)
!       if(tt > Zer) then
!         tp = dsqrt(One + tt*tt) - tt
!       else
!         tp =-dsqrt(One + tt*tt) - tt
!       endif
        tt = Adif / Aod
!       write(0,*)'alpha + beta adif < aod, tt = ',tt
        tp = dsqrt(One + tt*tt)
        if(tt .lt. Zer) tp = -tp
        tp = tp - tt
      else if (dabs(Adif).lt.Eps) then
        tp = Zer
      else
!       tt is tan(2x)
        tt = Aod / Adif
!       write(0,*)'alpha + beta adif >= aod, tt = ',tt
        if(dabs(tt) .gt. 1.d-7) then
          tp = (dsqrt(One + tt*tt) - One) / tt
        else
!         numerical stable:::Taylor expansion
          tp = 0.5d0*tt - 0.125d0*tt*tt + 0.0625d0*tt*tt*tt 
        endif
      endif
      cp = One / dsqrt(One + tp*tp)
      sp = tp * cp
!     angle alpha - beta
      Adif = Aii + Ajj
      Aod = Aji - Aij
      if(dabs(Adif) .lt. dabs(Aod)) then
        tt = Adif / Aod
!       write(0,*)'alpha - beta adif < aod, tt = ',tt
        tm = dsqrt(One + tt*tt)
        if(tt .lt. Zer) tm = -tm
        tm = tm - tt
      else if (dabs(Adif).lt.Eps) then
        tm = Zer
      else
        tt = Aod / Adif
!       write(0,*)'alpha - beta adif >= aod, tt = ',tt
        if(dabs(tt) .gt. 1.d-7) then
          tm = (dsqrt(One + tt*tt) - One) / tt
        else
!         numerical stable:::Taylor expansion
          tm = 0.5d0*tt - 0.125d0*tt*tt + 0.0625d0*tt*tt*tt 
        endif
      endif
      cm = One / dsqrt(One + tm*tm)
      sm = tm * cm
      c1 = cp*cm - sp*sm
      s1 = sp*cm + cp*sm
      c2 = cp*cm + sp*sm
      s2 = sp*cm - cp*sm
!     update U and A
      Do k = 1, N
        tmp=c1*U(k,j)-s1*U(k,i)
        U(k,i)=s1*U(k,j)+c1*U(k,i)
        U(k,j)=tmp
        tmp=c1*A(j,k)-s1*A(i,k)
!print*,A(i,k),A(j,k),s1,c1
        A(i,k)=s1*A(j,k)+c1*A(i,k)
        A(j,k)=tmp
!print*,'A_two',A(i,k),A(j,k),s1,c1,k
      enddo
!     update VT and A
      Do k = 1, N
        tmp=c2*VT(k,j)-s2*VT(k,i)
        VT(k,i)=s2*VT(k,j)+c2*VT(k,i)
        VT(k,j)=tmp
        tmp=c2*A(k,j)-s2*A(k,i)
        A(k,i)=s2*A(k,j)+c2*A(k,i)
        A(k,j)=tmp
!print*,'A_two_2',A(k,i),A(k,j)
      enddo

      A(j,i) = Zer
      A(i,j) = Zer

      return
      End
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
