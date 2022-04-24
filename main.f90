program Uebungsaufgabe
  implicit none
  real :: determinantA,determinantB,determinantC
  real :: start, finish
  real, allocatable :: A(:,:)         ! Matrix, für die die Determinante berechnet wird
  integer           ::  n             ! Größe der Matrix A
  integer           ::  j,i           ! Laufindizes

  open (10,file='matrixX.txt')
  read(10,*) n
  allocate(A(n,n))

  do i=1,n
    read(10,*) (A(i,j), j=1,n)
  enddo
  close(10)

  call cpu_time(start)
  call determinanteLaplace(A, n, determinantA)
  call cpu_time(finish)
  write(*,*) 'Determinante A; t=',finish-start

  call cpu_time(start)
  call determinanteLaplaceMitAbsuchen(A, n, determinantB)
  call cpu_time(finish)
  write(*,*) 'Determinante B; t=',finish-start

  call cpu_time(start)
  call determinanteDreiecksform(A, n, determinantC)
  call cpu_time(finish)
  write(*,*) 'Determinante C; t=',finish-start

  write(*,*) 'Determinante A=',determinantA
  write(*,*) 'Determinante B=',determinantB
  write(*,*) 'Determinante C=',determinantC

  read (*,*)

end program Uebungsaufgabe
