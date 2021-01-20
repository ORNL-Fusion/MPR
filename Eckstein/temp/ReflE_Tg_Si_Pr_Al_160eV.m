Mtg=28.085;  %Si
Mpr=12.01;   %C
Ztg=14;   %Si
Zpr=6;    %C
%not used in RE or RN, so I'm commenting them out here
%Esp=4.70;    %for erosion (chemical erosion, used in angular dep.)
%Esb=4.70;   %for particle emission

%energy param, C->Si data not available ; use N -> Si parameters
a1=-0.05303;
a2=-5.410;
a3=5.943;
a4=0.1234;

%angular dependence for E0=150eV does not exist;(d in here = c in Eckstein)
%option 1: use 200eV Si -> Si (KrC)
d1=0.9802;
d2=0.9801;
d3=2.861;
d4=-4.471; 
%option 2: use %option 2: use 4500eV Ar -> Si ; does not exist
