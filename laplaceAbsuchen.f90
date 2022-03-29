
subroutine determinanteLaplaceMitAbsuchen(A,n,det)

  implicit none

  real    :: A(n,n)          ! A(Spalte, Zeile)
  real    :: det             !
  real    :: laplaceNachSpalte            !
  integer :: n,spaltenNullen,zeilenNullen !
  integer :: spaltenKandidatNullen,zeilenKandidatNullen !
  integer :: spaltenKandidat,zeilenKandidat !
  integer :: i,j,k           !

  interface
    function TranspoMatrix(A,n)
      integer :: n
      real    :: A(n,n)
      real    :: TranspoMatrix(n-1,n-1)
    end function TranspoMatrix
  end interface
  ! Vorbelegung der Variablen
  !**************************************************************************!

  spaltenNullen = 0
  spaltenKandidatNullen = 0
  zeilenNullen = 0
  zeilenKandidatNullen = 0

  ! Durchsuchen der Spalten
  !
  do i = 1, n ! Jede Spalte
    do j = 1, n
      if ( A(i,j).eq.0 ) then
        spaltenNullen = spaltenNullen + 1
      end if
    end do
    if ( spaltenNullen.gt.spaltenKandidatNullen ) then ! Haben wir eine Spalten gefunden mit mehr nullen, dann wechseln wir den Kandidaten
      spaltenKandidat = i
      spaltenKandidatNullen = spaltenNullen
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
    if ( zeilenNullen.gt.zeilenKandidatNullen ) then ! Haben wir eine Spalten gefunden mit mehr nullen, dann wechseln wir den Kandidaten
      zeilenKandidat = j
      zeilenKandidatNullen = zeilenNullen
    end if
    zeilenNullen = 0
  end do

  ! Abprüfen des Sonderfalls n = kandidatNullen, also hat eine Zeile genauso viele Nullen wie sie Spalten oder Zeilen hat
  if ( zeilenKandidatNullen.eq.n.or.spaltenKandidatNullen.eq.n ) then
    det = 0
    return
  end if
  ! Eventuelles Transponieren der Matrix
  ! Wenn eine Zeile mehr nullen hat, als eine Spalte, dann müssen wir transponieren
  if(zeilenKandidatNullen.gt.spaltenKandidatNullen) then
    A = TranspoMatrix(A,n)
    det = laplaceNachSpalte(A,n,zeilenKandidat)
  else
    det = laplaceNachSpalte(A,n,spaltenKandidat)
  end if
end subroutine

recursive function laplaceNachSpalte(A, n, s) result(d)
  integer:: n,s        ! Groesse und Spalte der Matrix nach der entwickelt wird.
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

  !Hier wird nach der Spalte s entwickelt!
  do i = 1, n
    h = A(i,s)
    if(h.eq.0.0) then ! Sollte das Element 0.0 sein, dann sparen wir uns den Aufwand
      d = -d
      cycle ! Daher überspringen wir diese Iteration
    endif
    ! Sonst holen wir uns die Untermatrix
    B = UnterMatrix(A,n,i,s)

    ! Und davon letztlich die Determinante
    d = -d + h*laplace(B,n-1)
  end do
  if(mod(s,2).eq.1) then ! Sollten wir nach einer ungeraden Spalte entwickeln (s / 2 = x Rest 1 ergeben) dann tausche das Vorzeichen
    d = -d
  endif
end function laplaceNachSpalte
