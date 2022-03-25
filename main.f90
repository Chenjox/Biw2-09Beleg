program Uebungsaufgabe
  implicit none
  real :: determinant
  real, allocatable :: A(:,:)        ! Matrix, für die die Determinante berechnet wird
  integer           ::  n             ! Größe der Matrix A
  integer           ::  j,i           ! Laufindizes
  ! Zuerst den simplen LAPLACE'schen Entwicklungssatz
  ! Mittels Rekursion, iterativ geht das glaube ich nicht.

  open (10,file='matrixB.txt')
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

subroutine determinanteLaplace(A,n,det)

  implicit none

  real    :: A(n,n)          !
  real    :: det             !
  real    :: laplace         !
  integer :: n               !
  integer :: i,j,k           !

  ! Vorbelegung der Variablen
  !**************************************************************************!

  det = laplace(A,n)

end subroutine

recursive function laplace(A, n) result(d)
  integer:: n          ! Groesse
  real   :: A(n,n)     ! Die Matrix an sich.
  real   :: d         ! Determinante
  real                :: B(n-1,n-1) ! Die Untermatrix
  integer             :: i,j,k      ! Zeilen und Spaltenindezes
  interface
    function UnterMatrix(A,n,z,s)
      integer :: n,z,s
      real    :: A(n,n)
      real    :: UnterMatrix(n-1,n-1)
    end function UnterMatrix
  end interface

  d = 0.0

  if(n.eq.1) then !Basisfall ist abbruchbedingung
    d = A(1,1)
    return
  end if

  !Hier machen wir mal etwas
  do i = 1, n
    d = A(i,1)
    ! Dann holen wir uns die Untermatrix
    B = UnterMatrix(A,n,i,1)

    !write(*,*) 'A()=',A(i,1)
    !do k = 1, n-1
    !  write(*,*) (B(k,j), j=1,n-1)
    !end do

    d = -d + A(i,1)*laplace(B,n-1)
  end do
end function laplace

function UnterMatrix(A,n,z,s)
  integer :: n,z,s ! groesse, zeile und spalte die weggelassen werden
  real    :: A(n,n)
  real   :: UnterMatrix(n-1,n-1)
  integer :: i,j

  do i = 1, n
    if ( i.eq.z ) then ! Die Zeile mit Index i wird ignoriert
      continue
    end if
    do j = 1, n
      if ( j.eq.s ) then ! Die Spalte mit index j wird ignoriert
        continue
      else if ( j.gt.s ) then ! Ist j größer als m dann müssen wir immer eins von j abziehen
        if(i.gt.z) then
          UnterMatrix(i-1,j-1) = A(i,j)
        else
          UnterMatrix(i,j-1) = A(i,j)
        end if
      else ! j ist kleiner als s
        if(i.gt.z) then
          UnterMatrix(i-1,j) = A(i,j)
        else
          UnterMatrix(i,j) = A(i,j)
        end if
      end if
    end do
  end do
end function UnterMatrix
