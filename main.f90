program Uebungsaufgabe
  implicit none
  integer :: groesse
  real :: determinant, determinante_laplace
  real, allocatable:: A(:,:)        ! Matrix, für die die Determinante berechnet wird
  integer             n             ! Größe der Matrix A
  integer             j,i           ! Laufindizes
  ! Zuerst den simplen LAPLACE'schen Entwicklungssatz
  ! Mittels Rekursion, iterativ geht das glaube ich nicht.

  open (10,file='matrixB.txt')
  read(10,*) n
  allocate(A(n,n))

  ! Einlesen der Matrix und Schlie�en der Eingabedatei
  !**************************************************************************!
  do i=1,n
  		read(10,*) (A(i,j), j=1,n)
  enddo
  close(10)

  determinant = determinante_laplace(A, n)
  write(*,*) 'Determinante=',determinant
  read (*,*)

end program Uebungsaufgabe

recursive function determinante_laplace(A, n) result(d)
  integer, intent(in) :: n   ! Groesse
  real, intent(in) :: A(n,n) ! Die Matrix an sich.
  real :: d                  ! Determinante
  real :: B(n-1,n-1)         ! Die Untermatrix
  integer :: i,j             ! Zeilen und Spaltenindezes

  d=1

  if(n.eq.1) then !Basisfall ist abbruchbedingung
    d = A(1,1)
    return
  end if
  !Es wird nach der ersten Spalte entwickelt, koste es was es wolle
  spalten: do i=1,n
    if(A(1,i).eq.0.0) exit spalten
    !Zuerst Belegen wir B mit den notwendigen elementen
    unterspalten: do j=1,n
      if (j.eq.i) then
        exit unterspalten !Wenn i = j ist, dann überspringen wir das ganze
      else if (j.gt.i) then !Wenn j > i ist dann müssen wir in B einen Index verschieben
        B(i,j-1) = A(i,j)
      else !Wenn j < i ist, dann wird 1 zu 1 übertragen
        B(i,j) = A(i,j)
      end if
    end do unterspalten
    ! Am Ende kommt das hier raus
    d = -d + A(1,i) * determinante_laplace(B, n-1)
    !d=A(1,i)
  end do spalten
end function determinante_laplace
