
subroutine performance(maxGroesse, maxDurchlauf, dichteInkrement)
  implicit none
  integer :: maxGroesse
  integer :: maxDurchlauf
  integer :: groesse
  integer :: dichteInkrement
  real, allocatable :: matrix(:,:)
  integer :: durchlauf, maxDichte, dichte
  real :: zeiten(4)
  real :: minZeiten(4)
  real :: maxZeiten(4)
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
  maxZeiten = 0.0
  minZeiten = 1.5e+20
  !if(maxDichte.eq.0) maxDichte = 100 / dichteInkrement

  open(unit=30,file='performance.csv')
  write(unit=30,fmt='(A120)') 'groesse,dichte,Zeit1,maxZeit1,minZeit1,&
                                              Zeit2,maxZeit2,minZeit2,&
                                              Zeit3,maxZeit3,minZeit3,&
                                              Zeit4,maxZeit4,minZeit4'

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
        maxZeiten(1) = max(maxZeiten(1), finishTime - startTime)
        minZeiten(1) = min(minZeiten(1), finishTime - startTime)
        zeiten(1) = zeiten(1) + finishTime - startTime

        call cpu_time(startTime)
        call determinanteLaplaceMitAbsuchen(matrix, groesse, determinant)
        call cpu_time(finishTime)

        maxZeiten(2) = max(maxZeiten(2), finishTime - startTime)
        minZeiten(2) = min(minZeiten(2), finishTime - startTime)
        zeiten(2) = zeiten(2) + finishTime - startTime

        call cpu_time(startTime)
        call determinanteDreiecksform(matrix, groesse, determinant)
        call cpu_time(finishTime)

        maxZeiten(3) = max(maxZeiten(3), finishTime - startTime)
        minZeiten(3) = min(minZeiten(3), finishTime - startTime)
        zeiten(3) = zeiten(3) + finishTime - startTime

        call cpu_time(startTime)
        call determinanteDreiecksformPivotisierung(matrix, groesse, determinant)
        call cpu_time(finishTime)

        maxZeiten(4) = max(maxZeiten(4), finishTime - startTime)
        minZeiten(4) = min(minZeiten(4), finishTime - startTime)
        zeiten(4) = zeiten(4) + finishTime - startTime

      end do

      zeiten(1) = zeiten(1)/maxDurchlauf
      zeiten(2) = zeiten(2)/maxDurchlauf
      zeiten(3) = zeiten(3)/maxDurchlauf
      zeiten(4) = zeiten(4)/maxDurchlauf

      write (30, '(I2,",",I3,12(",", e32.16))') groesse, dichte, &
                                                zeiten(1), maxZeiten(1), minZeiten(1), &
                                                zeiten(2), maxZeiten(2), minZeiten(2), &
                                                zeiten(3), maxZeiten(3), minZeiten(3), &
                                                zeiten(4), maxZeiten(4), minZeiten(4)
      ! Zurücksetzen der Variablen
      zeiten = 0.0
      minZeiten = 1.5e+20
      maxZeiten = 0.0
      ! Die Matrix muss freigegeben werden.
      deallocate(matrix)
    end do
  end do

  close(unit=30)

end subroutine performance
