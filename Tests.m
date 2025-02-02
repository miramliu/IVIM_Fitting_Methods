
b = [0,50,100,150,200,250,300,500,700,900]
f = 0.07
D = 0.0004
Dstar=0.009
signal = f*exp(-b.*Dstar) + (1-f)*exp(-b.*D)


%test
Algorithm1(b,signal)
Algorithm2(b,signal)
Algorithm3(b,signal)
Algorithm4(b,signal)