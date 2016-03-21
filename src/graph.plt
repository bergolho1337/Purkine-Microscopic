set grid
set terminal png
set output "../Runs/Run5/Graphics/dndt1.png"
set title "Vol = 42"
plot "../Runs/Run5/data42.dat" using 1:4 title "m" w l
