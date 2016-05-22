set grid
set terminal png
set output "../Runs/Run0/Graphics/dndt1.png"
set title "Vol = 1122"
plot "../Runs/Run0/data1122.dat" using 1:4 title "m" w l
