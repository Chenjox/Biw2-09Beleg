
program hauptzwei

  implicit none
  ! maxGroesse, durchläufe, dichteInkrement
  call performance(10, 100, 20)

  call system('texify -p --tex-option="-interaction=nonstopmode" --engine=pdftex --tex-option="--synctex=1" "auswertung.tex"')

end program hauptzwei
