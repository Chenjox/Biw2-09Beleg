
! Siehe Vorlesung 6

subroutine determinanteDreiecksform(A,n,det)

  implicit none

  real    A(n,n)          ! Matrix, für die die Determinante berechnet wird
  real    A_diag(n,n)     ! in Dreieckform Überführte Matrix A
  real    det             ! Determinante der Matrix A bzw. A_diag
  real    m,merk          ! Hilfsvariablen
  integer n               ! Größe der Matrix A
  integer i,j,k           ! Laufindizes
  logical det_exist       ! Abbruchvariable

  ! Vorbelegung der Variablen
  !**************************************************************************!
  det_exist=.true.
  det=1
  A_diag=A

  ! Überführung der Matrix A in Dreieckform (ohne Pivotisierung)
  !**************************************************************************!
  do k=1,n-1
    if (A_diag(k,k).eq.0) then                                       ! wenn das aktuelle Diagonalelement = 0,
                                                                 ! dann wird Existenz der Determinante
      det_exist=.false.                                          ! zunächst angezweifelt und die
      do i=k+1,n                                                 ! nächste Zeile gesucht, in der das
        if(A_diag(i,k).ne.0) then                          ! entsprechende Element |= 0 ist,
          do j=1,n                                       ! um diese Zeile mit Zeile k zu tauschen
            merk=A_diag(i,j)                           !
            A_diag(i,j)=A_diag(k,j)                    ! Tausch von Zeile i und k
            A_diag(k,j)=merk                           !
          enddo
          det_exist=.true.                               ! Zweifel an der Existenz der
                                                         ! Determinante ist beseitigt
          det=-det                                       ! Vorzeichenwechsel der Determinante
                                                         ! bei Zeilentausch
          exit                                           ! vorzeitiges Beenden der Schleife
        endif
      enddo
      if(det_exist.eqv..false.) then                             ! alle Elemente der Spalte k sind Null,
        det=0.0                                            ! das heißt det(A)=0
        return                                             ! Rückkehr ins aufrufende (Unter-)Programm
      endif
    endif

    do j=k+1,n                                                       ! Subtrahieren der k-ten Zeile (mit
      m=A_diag(j,k)/A_diag(k,k)                                  ! Faktor m) von den restlichen Zeilen
      do i=k+1,n                                                 ! = Diagonalisieren
        A_diag(j,i)=A_diag(j,i)-m*A_diag(k,i)
      enddo
    enddo
  enddo

  ! Berechnung der Determinante als Produkt der Hauptdiagonalelemente
  !**************************************************************************!
  do i=1,n
   det=det*A_diag(i,i)
  enddo

end subroutine
