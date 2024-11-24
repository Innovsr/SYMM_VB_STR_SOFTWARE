!---------------------------------------------
      SUBROUTINE Invmat(Ndim,Ifail)
!---------------------------------------------
!
!
!
use commondat
use commondat1
      Implicit DOUBLE PRECISION(A-H,O-Z)
      Dimension A(1000,1000),B(1000,1000)
      integer::ndim,ifail



!print*,'enter Invmat'
do i=1,ndim
do j=1,ndim
A(i,j)=ind_mat(i,j)
!print*,A(i,j)
enddo
enddo
!do i=1,ndim
!write(9,111),'ind_mat',(int(A(i,j)),j=1,ndim)
!enddo
111 format(a,30I5)

      N=Ndim
      Ifail = 0
      Xdet = 1.D0
      Do I = 1,N
       Do J = 1,N
        B(J,I) = 0.D0
       Enddo
       B(I,I) = 1.D0
      Enddo
!do i=1,ndim
!print*,'A(i,j)',(A(i,j),j=1,ndim)
!enddo

      Do I = 1,N-1
!--Find the Max
       Amax = Dabs(A(I,I))
!print*,'AmaxII',I,Amax
       Nmax = I

       Do J = I+1,N
        Aij = Dabs(A(I,J))
!print*,'Aij',I,J,Aij
        If(Aij.Gt.Amax) Then
         Amax = Aij
         Nmax = J
        Endif
       Enddo

!print*,'Amaxij',Amax,Nmax
       If(Amax.Lt.1.D-11) Then
!       If(Amax.Lt.1.D-1) Then
        Ifail = 1
!print*,'Amax',Amax
        Xdet = 0.D0
        Return
       Endif

       If(Nmax.Ne.I) Then
        Xdet = -Xdet
        Do J = 1,N
         Atmp = A(J,I)
         A(J,I) = A(J,Nmax)
         A(J,Nmax) = Atmp
         Atmp = B(J,I)
         B(J,I) = B(J,Nmax)
         B(J,Nmax) = Atmp
        Enddo
       Endif

!-- Vanish Lower
       Do J = I+1,N
        Atmp = A(I,J)/A(I,I)
        Do K = 1,N
!print*,'A(K,J)',I,J,K,A(K,J),Atmp,A(K,I)
         A(K,J) = A(K,J)-Atmp*A(K,I)
         B(K,J) = B(K,J)-Atmp*B(K,I)
!print*,'A(K,J)**',I,J,K,A(K,J),Atmp,A(K,I)
!do i1=1,ndim
!print*,'A(i,j)',(A(i1,j1),j1=1,ndim)
!enddo
        Enddo
       Enddo
        Do K = 1,N
!print*,'A(K,J)',k,(A(K,J),J=1,N)
        Enddo
      Enddo
!do i=1,ndim
!print*,'A(i,j)',(A(i,j),j=1,ndim)
!enddo

      If(Dabs(A(N,N)).Lt.1.D-11) Then
!      If(Dabs(A(N,N)).lt.1.D-1) Then
       Ifail = 1
!print*,'Dabs(A(N,N))',N,Dabs(A(N,N))
       Xdet = 0.D0
       Return
      Endif

      Do I = 1,N
       Xdet = Xdet*A(I,I)
      Enddo

!-- Vanish Upper
      Do I = N,2,-1
       Do J = I-1,1,-1
        Atmp = A(I,J)/A(I,I)
        Do K = 1,N
         A(K,J) = A(K,J)-Atmp*A(K,I)
         B(K,J) = B(K,J)-Atmp*B(K,I)
!print*,'A(K,J)****',I,J,K,A(K,J),Atmp,A(K,I)
        Enddo
       Enddo
      Enddo

!do i=1,ndim
!print*,'A(i,j)',(A(i,j),j=1,ndim)
!enddo
!-- Normalized
      Do I = 1,N
       Atmp = A(I,I)
       Do J = 1,N
        B(J,I) = B(J,I)/Atmp
       Enddo
      Enddo
      End

