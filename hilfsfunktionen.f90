
function UnterMatrix(A,n,z,s)
  implicit none
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
