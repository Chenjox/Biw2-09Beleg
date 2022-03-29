program Uebungsaufgabe
  implicit none
  real :: determinant
  real, allocatable :: A(:,:)        ! Matrix, für die die Determinante berechnet wird
  integer           ::  n             ! Größe der Matrix A
  integer           ::  j,i           ! Laufindizes
  ! Zuerst den simplen LAPLACE'schen Entwicklungssatz
  ! Mittels Rekursion, iterativ geht das glaube ich nicht.

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
  real   :: d,h        ! Determinante und Platzhalter
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

  !Hier wird nach der ersten Spalte entwickelt!
  do i = 1, n
    h = A(i,1)
    if(h.eq.0.0) then ! Sollte das Element 0.0 sein, dann sparen wir uns den Aufwand
      d = -d
      cycle ! Daher überspringen wir diese Iteration
    endif
    ! Sonst holen wir uns die Untermatrix
    B = UnterMatrix(A,n,i,1)

    ! Und davon letztlich die Determinante
    d = -d + h*laplace(B,n-1)
  end do
  d = -d
end function laplace

function UnterMatrix(A,n,z,s)
  integer :: n,z,s ! groesse, zeile und spalte die weggelassen werden
  real    :: A(n,n)
  real   :: UnterMatrix(n-1,n-1)
  integer :: i,j

  !write(*,*) 'z=',z,',s=',s
  !do i = 1, n
  !  write(*,*) (A(i,j), j=1,n)
  !end do

  zeilen: do i = 1, n
    if ( i.eq.z ) then ! Die Zeile mit Index i wird ignoriert
      cycle zeilen
    end if
    spalten: do j = 1, n
      if ( j.eq.s ) then ! Die Spalte mit index s wird ignoriert
        cycle spalten
      else if ( j.gt.s ) then ! Ist j größer als s dann müssen wir immer eins von j abziehen
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
    end do spalten
  end do zeilen

  !do i = 1, n-1
  !  write(*,*) (UnterMatrix(i,j), j=1,n-1)
  !end do
end function UnterMatrix
