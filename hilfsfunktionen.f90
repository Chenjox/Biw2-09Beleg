
function UnterMatrix(A,n,z,s)
  implicit none
  integer :: n,z,s ! groesse, zeile und spalte die weggelassen werden
  real    :: A(n,n)
  real   :: UnterMatrix(n-1,n-1)
  integer :: i,j

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

end function UnterMatrix

function TranspoMatrix(A,n)
  implicit none
  integer :: n,i,j
  real :: A(n,n)
  real :: TranspoMatrix(n,n)

  do i = 1, n
    do j = 1, n
      TranspoMatrix(i,j) = A(j,i)
    end do
  end do
end function TranspoMatrix

subroutine zufallMatrix(A,n,dichte)
implicit none
integer n,i,k
real u,v
real dichte

  Z=0
  do i=1,n
      do k=1,n
          call random_number(u)    !Funktion erzeugt Pseudozufallszahl von 0 bis 1
          v = FLOOR(101*u)        !Zahlen von 0 bis 100, die Zufallszahl mal 101 erzeugt Zahl von 1-100, FLOOR() rundet ab
          if(v.lt.dichte) then    !Falls die Zufallszahl kleiner als unsere geforderte Dichte, dann Erzeuge in der Zelle eine Null
              Z(i,k)=0
          else                    !Sonst fülle die Zelle mit einer Zufallszahl
              call random_number(u)
              v = FLOOR(9*u) + 1    !Rechen die Zufallszahl +1 damit keine Null erzeugt wird
              Z(i,k)=v
          endif
      enddo
  enddo

end subroutine
