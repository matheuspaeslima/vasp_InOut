subroutine LU(A,B,x,N)
implicit none
integer N
real*8 A(N,N),B(N),x(N)
!internal variables
real*8 L(N,N), U(N,N), P(N,N)
real*8 B2(N)
real*8 s, y(N)
real*8  piv
integer ipiv
logical long
integer i, j, k
long=.false.
! inicio fatoracao LU
!iniciando matriz de permutacao
P=0.0
do i=1,N
    P(i,i)=1.d0
enddo
if (long) then
print*, "============================================="
print*, "        inicio da Eliminacao de Gauss        "
endif
do k=1,N-1 ! indice coluna
if (long) then
print*, "                        "
print*, " => eliminacao coluna", k
endif
    !begin PIVOTEAMENTO =====
    ipiv=k
    piv=A(k,k)
    do j=k+1,N
        if (abs(A(j,k)).gt.abs(piv)) then 
            ipiv=j
            piv=A(j,k)
        endif
    end do
    if (long) print*, "piv, ipiv=", piv, ipiv
    ! permutar linha j e ipiv
    call permut(A,k,ipiv,N)
    call permut(P,k,ipiv,N)
    !end PIVOTEAMENTO =====
    if (long) print*, "A(k,k)=", A(k,k)
    do i=k+1,N ! indice linha
        A(i,k)=A(i,k)/A(k,k)
        do j=k+1,N
        A(i,j)=A(i,j)-A(i,k)*A(k,j)
        end do
    enddo
    if (long) then 
    print*, "matriz A permutada e eliminada:"
    do i=1,N
        print 15, A(i,:)
    enddo
    endif
enddo
if (long) then
print*, "        inicio da Eliminacao de Gauss        "
print*, "============================================="
endif
! calculo de L e U
U=A
L=0.
do i=1,N
    L(i,i)=1.d0
enddo

do i=2,N
    do j=1,i-1
        U(i,j)=0.d0
        L(i,j)=A(i,j)
    end do
end do

if (long) then
print*, "============================================="
print*, "          fatores L e U                      "
print*," L"
do i=1,N
    print 15, L(i,:)
enddo
print*, " "
print*," U"
do i=1,N
    print 15, U(i,:)
enddo
!   fim fatoracao LU
print*, "============================================="
! calculo de B2=PB
print*, "=========== triangular inferior ============="
endif
B2=0.d0
do i=1,N
do j=1,N
B2(i)=B2(i)+P(i,j)*B(j)
end do
end do
if (long) then
print*, " matriz de permutacao:"
print*," P"
do i=1,N
    print 15, P(i,:)
enddo
print*, "  "
print'(a4,20f10.6)',"B=", B(:)
print*, " "
print'(a4,20f10.6)',"PB=", B2(:)
print*, " "

print*, " solving: b=Ly"
endif
do i=1,N
    s=0.
    do j=1,i-1
        s=s+y(j)*L(i,j)
    enddo
    y(i)=B2(i)-s
enddo
if (long) print'(a4,20f10.6)',"y=", y(:)


if (long) print*, " solving: y=Ux"
do i=N,1,-1
    s=0.
    do j=i,N-1
        s=s+x(j+1)*U(i,j+1)
    enddo
    x(i)=(y(i)-s)/U(i,i)
enddo

15 format(20f10.6)
end subroutine LU

subroutine permut(A,i,j,N)
implicit none
integer i,j,N
real*8 A(N,N)
real*8 dummy

integer k, l
if (i==k) return
do k=1,N
    dummy=A(i,k)
    A(i,k)=A(j,K)
    A(j,k)=dummy
enddo
return
end


subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
