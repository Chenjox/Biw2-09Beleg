
subroutine performance(maxGroesse, maxDurchlauf, dichteInkrement)
  implicit none
  integer :: maxGroesse
  integer :: maxDurchlauf
  integer :: groesse
  integer :: dichteInkrement
  real, allocatable :: matrix(:,:)
  integer :: durchlauf, maxDichte, dichte
  real :: zeiten(3)
  real :: startTime, finishTime, timeTaken
  real :: determinant

  timeTaken = 0.0
  startTime = 0.0
  finishTime = 0.0
  determinant = 0.0
  durchlauf = 0
  groesse = 0

  ! maxDichte
  maxDichte = 100
  zeiten = 0.0
  !if(maxDichte.eq.0) maxDichte = 100 / dichteInkrement

  open(unit=30,file='performance.csv')
  write(unit=30,fmt='(A34)') 'groesse,dichte,zeit1,zeit2,zeit3'

  do dichte = dichteInkrement, maxDichte, dichteInkrement
    write(*, *) dichte
    do groesse = 1, maxGroesse
      write(*,*) groesse
      allocate(matrix(groesse,groesse)) ! Die Matrix ist groesse x groesse
      do durchlauf=1, maxDurchlauf
        call zufallMatrix(matrix,groesse,dichte) ! Zufallige belegung der Matrix
        call cpu_time(startTime)
        call determinanteLaplace(matrix, groesse, determinant)
        call cpu_time(finishTime)
        ! Wir berechnen es inplace
        zeiten(1) = zeiten(1) + finishTime - startTime

        call cpu_time(startTime)
        call determinanteLaplaceMitAbsuchen(matrix, groesse, determinant)
        call cpu_time(finishTime)

        zeiten(2) = zeiten(2) + finishTime - startTime

        call cpu_time(startTime)
        call determinanteDreiecksform(matrix, groesse, determinant)
        call cpu_time(finishTime)

        zeiten(3) = zeiten(3) + finishTime - startTime

      end do

      zeiten(1) = zeiten(1)/maxDurchlauf
      zeiten(2) = zeiten(2)/maxDurchlauf
      zeiten(3) = zeiten(3)/maxDurchlauf

      write (30, '(I2,",",I3,3(",", F16.8))') groesse, dichte, zeiten(1), zeiten(2), zeiten(3)

      zeiten = 0.0
      ! Die Matrix muss freigegeben werden.
      deallocate(matrix)
    end do
  end do

  close(unit=30)

end subroutine performance
