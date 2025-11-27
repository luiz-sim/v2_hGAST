set xr[-150:150]
set yr[-150:150]

set title 'Modeshape 1 | Tower Fore-aft'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR001.dat' u 1:2:3 w p t 'Deformed'
pause -1

set title 'Modeshape 2 | Tower Side-to-Side'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR002.dat' u 1:2:3 w p t 'Deformed'
pause -2

set title 'Modeshape 3 | Edgewise Symmetric'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR003.dat' u 1:2:3 w p t 'Deformed'
pause -3

set title 'Modeshape 4 | Flapwise Asymmetric Yaw'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR004.dat' u 1:2:3 w p t 'Deformed'
pause -4

set title 'Modeshape 5 | Flapwise Asymmetric Tilt'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR005.dat' u 1:2:3 w p t 'Deformed'
pause -5

set title 'Modeshape 6 | Flapwise Symmetric'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR006.dat' u 1:2:3 w p t 'Deformed'
pause -6

set title 'Modeshape 7 | Edgewise Asymmetric Vertical'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR007.dat' u 1:2:3 w p t 'Deformed'
pause -7

set title 'Modeshape 8 | Edgewise Asymmetric Horizontal'
sp 'modeR000.dat' u 1:2:3 w p t 'Undeformed', 'modeR008.dat' u 1:2:3 w p t 'Deformed'
pause -8
