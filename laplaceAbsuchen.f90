
subroutine determinanteLaplaceMitAbsuchen(A,n,det)

  implicit none

  real    :: A(n,n)          ! A(Spalte, Zeile)
  real    :: det             !
  real    :: laplaceMitAbsuchen            !
  integer :: n

  ! Vorbelegung der Variablen
  !**************************************************************************!

  det = laplaceMitAbsuchen(A,n)

end subroutine

recursive function laplaceMitAbsuchen(A, n) result(d)

  implicit none
  integer:: n          ! Groesse und Spalte der Matrix nach der entwickelt wird.
  real   :: A(n,n)     ! Die Matrix an sich.
  real   :: T(n,n)     ! Die eventuell transponierte Matrix
  real   :: d,h        ! Determinante und Platzhalter
  integer :: spaltenNullen,zeilenNullen ! Spaltenabsuche
  integer :: spaltenKandidatNullen,zeilenKandidatNullen ! Zeilenabsuche
  integer :: spaltenKandidat,zeilenKandidat,s ! die Tatsächliche Zeile oder Spalte
  real                :: B(n-1,n-1) ! Die Untermatrix
  integer             :: i,j        ! Zeilen und Spaltenindezes
  interface ! Interfaces damit der gfortran Kompiler nicht meckert
    function UnterMatrix(A,n,z,s)
      integer :: n,z,s
      real    :: A(n,n)
      real    :: UnterMatrix(n-1,n-1)
    end function UnterMatrix
    function TranspoMatrix(A,n)
      integer :: n
      real    :: A(n,n)
      real    :: TranspoMatrix(n,n)
    end function TranspoMatrix
  end interface

  d = 0.0
  if(n.eq.1) then !Basisfall ist abbruchbedingung
    d = A(1,1)    !wenn Matrixgröße=1, Det=Wert der Matrix
    return
  end if

  spaltenNullen = 0   !Startwert Spaltennullen = 0
  spaltenKandidat = 1 !Spalte mit aktuell meisten Nullen
  spaltenKandidatNullen = 0 !Wie viele Nullen hat aktueller Kandidat
  zeilenNullen = 0
  zeilenKandidat = 1
  zeilenKandidatNullen = 0

  ! Durchsuchen der Spalten
  !
  do i = 1, n ! Spalten werden von 1 bis n hochgezählt
    do j = 1, n !Zeilen für aktuelle Spalte werden hochgezählt, bis alle spalten abgesucht sind
      if ( A(i,j).eq.0.0 ) then
        spaltenNullen = spaltenNullen + 1
      end if
    end do
    if ( spaltenNullen.gt.spaltenKandidatNullen ) then ! Haben wir eine Spalten gefunden mit mehr nullen, dann wechseln wir den Kandidaten
      spaltenKandidat = i !Spalte in der wir uns gerade befinden ist neuer Spaltenkndidat
      spaltenKandidatNullen = spaltenNullen !Nullen in Kandidat sind Spaltennullen
    end if
    spaltenNullen = 0
  end do

  ! Durchsuchen der Zeilen
  !
  do j = 1, n ! Jede Zeile
    do i = 1, n
      if ( A(i,j).eq.0 ) then
        zeilenNullen = zeilenNullen + 1
      end if
    end do
    if ( zeilenNullen.gt.zeilenKandidatNullen ) then ! Haben wir eine Zeile gefunden mit mehr nullen, dann wechseln wir den Kandidaten
      zeilenKandidat = j
      zeilenKandidatNullen = zeilenNullen
    end if
    zeilenNullen = 0
  end do

  ! Abprüfen des Sonderfalls n = kandidatNullen, also hat eine Zeile genauso viele Nullen wie sie Spalten oder Zeilen hat
  if ( zeilenKandidatNullen.eq.n.or.spaltenKandidatNullen.eq.n ) then
    d = 0
    return
  end if
  ! Eventuelles Transponieren der Matrix
  ! Wenn eine Zeile mehr nullen hat, als eine Spalte, dann müssen wir transponieren
  !s ist zeile oder spalte nach der entwicklt wird
  if(zeilenKandidatNullen.gt.spaltenKandidatNullen) then
    T = TranspoMatrix(A,n)
    s = zeilenKandidat !wenn zeilenkandidat größer als spaltenkandidat ist s zeilenkandidat
  else
    T = A
    s = spaltenKandidat !wenn spaltenkandidat größer als zeilenkandidat ist s spaltenkandidat
  end if
  !Hier wird endlich nach der Spalte s entwickelt!
  do i = 1, n
    h = T(i,s)
    if(h.eq.0.0) then ! Sollte das Element 0.0 sein, dann sparen wir uns den Aufwand
      d = -d
      cycle ! Daher überspringen wir diese Iteration
    endif
    ! Sonst holen wir uns die Untermatrix
    B = UnterMatrix(T,n,i,s)

    ! Und davon letztlich die Determinante
    d = -d + h*laplaceMitAbsuchen(B,n-1)
  end do
  if(mod(s,2).eq.1) then ! Sollten wir nach einer ungeraden Spalte entwickeln (s / 2 = x Rest 1 ergeben) dann tausche das Vorzeichen
    d = -d
  endif
end function laplaceMitAbsuchen
