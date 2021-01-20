Mtg=28.085;  %Si
Mpr=12.01;   %C
Ztg=14;   %Si
Zpr=6;    %C
%surface binding energy, eV
Esp=4.70; %for erosion (chemical erosion, used in angular dep.)
Esb=4.70; %for particle emission

%energy dependence ; C->Si data not available ; use N -> Si parameters  
lambda=0.4888;
qeros=1.4367;
mu_eros=1.7970;
Eth=16.6977;

%angular dependence for E0=150eV does not exist;
%option 1: use 200eV Si -> Si (KrC) 
%f=11.0479;
%b=6.1405;
%c=0.6299;
%option 2: use 4500eV Ar -> Si ; 
f=3.1385;
b=0.8328;
c=0.9440;