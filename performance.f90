
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
  maxDichte = mod(100,dichteInkrement) * 10

  do groesse = 1, maxGroesse
    allocate(matrix(groesse,groesse)) ! Die Matrix ist groesse x groesse
    do dichte = dichteInkrement, maxDichte, dichteInkrement
      do durchlauf=1, maxDurchlauf
        call zufallMatrix(matrix,groesse,20) ! Zufallige belegung der Matrix
        call cpu_time(startTime)
        call determinanteLaplace(matrix, groesse, determinant)
        call cpu_time(finishTime)
        timeTaken = timeTaken + finishTime - startTime
      end do
      timeTaken = timeTaken / maxDurchlauf
      zeiten(1) = timeTaken ! Schreiben in ein Array
      ! TODO Hier kommt noch das schreiben in eine Datei
      timeTaken = 0.0 ! Alle Variablen zurücksetzen!
      do durchlauf=1, maxDurchlauf
        call zufallMatrix(matrix,groesse,20) ! Zufallige belegung der Matrix
        call cpu_time(startTime)
        call determinanteLaplaceMitAbsuchen(matrix, groesse, determinant)
        call cpu_time(finishTime)
        timeTaken = timeTaken + finishTime - startTime
      end do
      timeTaken = timeTaken / maxDurchlauf
      zeiten(2) = timeTaken
      ! TODO Hier kommt noch das schreiben in eine Datei
      timeTaken = 0.0 ! Alle Variablen zurücksetzen!
      do durchlauf=1, maxDurchlauf
        call zufallMatrix(matrix,groesse,20) ! Zufallige belegung der Matrix
        call cpu_time(startTime)
        call determinanteDreiecksform(matrix, groesse, determinant)
        call cpu_time(finishTime)
        timeTaken = timeTaken + finishTime - startTime
      end do
      timeTaken = timeTaken / maxDurchlauf
      zeiten(3) = timeTaken
      timeTaken = 0.0 ! Alle Variablen zurücksetzen!

      write (1, '(1x, F, 3(",", F))') groesse, dichteInkrement, zeiten(1), zeiten(2), zeiten(3)
    end do
    deallocate(matrix) ! Die Matrix muss freigegeben werden.
  end do
  !write (1, '(1x, F, 3(",", F))')
  !sehr dünn besetzt



  !dünn besetzt

  !halb besetzt

  !dicht besetzt

  !nahezu voll besetzt

end subroutine performance
