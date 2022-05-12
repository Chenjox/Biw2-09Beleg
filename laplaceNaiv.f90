
subroutine determinanteLaplace(A,n,det)

  implicit none

  real    :: A(n,n)          !
  real    :: det             !
  real    :: laplace         !
  integer :: n               !

  ! Vorbelegung der Variablen
  !**************************************************************************!

  det = laplace(A,n)

end subroutine

recursive function laplace(A, n) result(d)
  integer :: n          ! Groesse
  real    :: A(n,n)     ! Die Matrix an sich.
  real    :: d,h        ! Determinante und Platzhalter
  real                :: B(n-1,n-1) ! Die Untermatrix
  integer             :: i,j     ! Zeilen und Spaltenindezes
  interface
    function UnterMatrix(A,n,z,s)
      integer :: n,z,s                !Größe Matrix
      real    :: A(n,n)               !Matrix an sich
      real    :: UnterMatrix(n-1,n-1) !Untermatrix
    end function UnterMatrix
  end interface

  d = 0.0 !Determinante Startwert 0

  if(n.eq.1) then !Basisfall ist abbruchbedingung
    d = A(1,1) ! wenn Matrix Größe 1 hat, ist der Wert die Determinante und wir sind fertig
    return
  end if

  !Hier wird nach der ersten Spalte entwickelt!
  do i = 1, n !Zeilen werden von 1 bis n hochgezählt
    h = A(i,1) !Wert an Stelle (i,1)
    if(h.eq.0.0) then ! Sollte das Element 0.0 sein, dann sparen wir uns den Aufwand
      d = -d
      cycle ! Daher überspringen wir diese Iteration
    endif
    ! Sonst holen wir uns die Untermatrix
    B = UnterMatrix(A,n,i,1)

    ! Und davon letztlich die Determinante
    d = -d + h*laplace(B,n-1)
    ! Wert bei 2. Spalte 1. Zeile * Determinante der Untermatrix (2. Spalte und 1. Zeile weggestrichen)
  end do
  d = -d
end function laplace
