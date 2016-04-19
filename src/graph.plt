set grid
set terminal png
set output "../Runs/Run0/Graphics/dndt1.png"
set title "Vol = 480"
plot "../Runs/Run0/data480.dat" using 1:4 title "m" w l
