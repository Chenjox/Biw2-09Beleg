program Uebungsaufgabe
  implicit none
  real :: determinant
  real, allocatable :: A(:,:)        ! Matrix, für die die Determinante berechnet wird
  integer           ::  n             ! Größe der Matrix A
  integer           ::  j,i           ! Laufindizes

  open (10,file='matrixA.txt')
  read(10,*) n
  allocate(A(n,n))

  do i=1,n
      read(10,*) (A(i,j), j=1,n)
  enddo
  close(10)

  call determinanteLaplace(A, n, determinant)

  write(*,*) 'Determinante=',determinant

  read (*,*)

end program Uebungsaufgabe
