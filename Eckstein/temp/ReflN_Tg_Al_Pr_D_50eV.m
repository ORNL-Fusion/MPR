Mtg=26.981;  %Al
Mpr=2.014;   %D
Ztg=13;   %Al
Zpr=1;    %D
%not used in RE or RN, so I'm commenting them out here
%Esp=1.0; %for erosion (chemical erosion, used in angular dep.)
%Esb=3.36; %for energy of particle emission  (thompson distribution)

%energy param, for D -> Al  (b here  = a in Eckstein)
b1=0.2319;
b2=-0.1777;
b3=0.1523;
b4=1.523;

%%angular param for D->Al at appropriate E UNAVAILABLE: 
%option 1, for E0=13333eV D->Si  
%c1=1.177;
%c2=1.162;
%c3=2.026;
%c4=-3.532;

%option 2, for E0=200eV, He->Si (no params for D->Al at low E)
c1=5.273;
c2=5.018;
c3=1.238;
c4=-3.091; 