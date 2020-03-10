Mtg=28.085;  %Si
Mpr=12.01;   %C
Ztg=14;   %Si
Zpr=6;    %C
%surface binding energy, eV
Esp=0.0;  %for reflection (chemical effects)
Esb=4.70; %for particle emission

%energy param, C->Si data not available ; use N -> Si parameters 
%(b here  = a in Eckstein)
b1=-0.9449e-4;
b2=-1.398;
b3=3.847;
b4=0.1785;

%angular dependence for E0=150eV does not exist; use 200eV Si -> Si (KrC)
c1=0.8545;
c2=0.8531;
c3=2.073;
c4=2.926; %changed sign to avoid divergence