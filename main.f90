program Uebungsaufgabe
  implicit none
  real :: determinantA,determinantB,determinantC
  real, allocatable :: A(:,:)        ! Matrix, für die die Determinante berechnet wird
  integer           ::  n             ! Größe der Matrix A
  integer           ::  j,i           ! Laufindizes

  open (10,file='matrixC.txt')
  read(10,*) n
  allocate(A(n,n))

  do i=1,n
      read(10,*) (A(i,j), j=1,n)
  enddo
  close(10)

  call determinanteLaplace(A, n, determinantA)
  call determinanteLaplaceMitAbsuchen(A, n, determinantB)
  call determinanteDreiecksform(A, n, determinantC)

  write(*,*) 'Determinante A=',determinantA
  write(*,*) 'Determinante B=',determinantB
  write(*,*) 'Determinante C=',determinantC

  read (*,*)

end program Uebungsaufgabe
