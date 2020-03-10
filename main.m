%%MPR Trench, for micro-trench erosion studies
%%A. Lasa & J. Coburn

%%Start with a trench surface in 2D based on 
%%hyperbolic tangent function Z = A/2*(tanh(S(x-bx))+tanh(S(y-by)))
%%parameterized by A, S, bx, & by.

%%launch particles from a fixed z0 close to max(Z)
%%and (x0,y0) randomly chosen or forming an uniform grid

%%straight trajectories, specified angular distribution & energy


clear all values
clear variables
%close all

%%%% output folder
currentfolder = pwd;
outfolder=[currentfolder,'/output_files'];
mkdir(outfolder);
addpath(outfolder);


%%%%%%%%%%     DEFINE CASE     %%%%%%%%%%

%%1-DEFINE SURFACE
%amplitude of trench height, (e.g., in um)
A = 3;     %2

%slope parameters; 
% each side follows the formula: 
%       y = A*tanh(S*x) + B*tanh(R*x - C)
% that defines a trench, with sides that have 2 slopes: the main (S) and 
% a secondary slope, close to the edge
% A is the height of trench; B<A; C is the offset (where the 2nd slope starts)
% params depend on FiB & material (S should range from 1 to ~4 or 5) 
  
% downstream in X
SX1 = 8.0 ; %4.0 primary slope      
RX1 = 1.0 ; %1.0 2nd slope
BX1 = 1.0 ; %0.8 2nd depth 
CX1 = 2.0 ; %2.0 2nd length
% upstream in X ;  
SX2 = 8.0 ; %2.0
RX2 = 0.0 ; 
BX2 = 0.0 ;
CX2 = 0.0 ;
% downstream in Y
SY1 = 8.0; 
RY1 = 1.0; 
BY1 = 1.0; 
CY1 = 2.0; 
% upstream in Y 
SY2 = 8.0 ; 
RY2 = 0.0 ; 
BY2 = 0.0 ; 
CY2 = 0.0 ;

%trench length:
bx = 15;      %1/2 average length in x, (e.g., in um)
by = 15;      %1/2 average length in y, (e.g., in um)
surfxmin=-20; %surface x-min
surfxmax=24;  %surface x-max
surfymin=-20; %surface y-min
surfymax=24;  %surface y-max


%%2-DEFINE INITIAL CONDITIONS OF PARTICLES 

%trajectories
phi=5/36.0;       %phi = angle wrt x-axis 0 < phi < pi/2 
dlt=1.5;        %delta = angle wrt -z axis, (pointing to surface); 0<delta<pi/2
th=pi-dlt ;     %theta =angle wrt +z axis ; pi/2 < theta < pi
distr='Chrobk';      %'Boro88'; %Curr85, Boro85, Boro88, Boro89, Chrobak or blank. Replaces dlt & th if not blank

%launching area
initxmin=-20;   %-9; %x-min of initializing ('launching') particles
initxmax=16;     %x-max of initializing particles
initymin=-20;   %y-min of initializing particles
initymax=16;     %y-max of initializing particles
z0=0.01; %specified height of initializing particles

%number of 'particles'
%For micro-trench studies, Np=300000, nsteps=1500 looks good
NP=360000; 
%nsteps = average #impacts per cell, in a flat surface: i.e., represents statistics
nsteps=1800; 
%resolution = number of surface grids;
npoints=floor(NP/nsteps); 


%%3-DEFINE MATERIALS (for Eckstein's fit formula)

Tg='C'; %target material ; e.g. 'W'
Pr='D';  %projectile ; e.g. 'D'
E0=100.0; %impact energy, eV

%%load erosion and reflection parameters
Ecksteinfolder=[currentfolder,'/Eckstein/temp'];
addpath(Ecksteinfolder);

ErosParamFile=[Ecksteinfolder,'/Eros_','Tg_',Tg,'_Pr_',Pr, '_', num2str(E0), 'eV' ];
run(ErosParamFile)
%S0=['running ', ErosParamFile];
%disp(S0)

ReflParamFile=[Ecksteinfolder,'/ReflN_','Tg_',Tg,'_Pr_',Pr, '_', num2str(E0), 'eV'];
run(ReflParamFile)
%S1=['running ', ReflParamFile];
%disp(S1)

ReflParamFile=[Ecksteinfolder,'/ReflE_','Tg_',Tg,'_Pr_',Pr, '_', num2str(E0), 'eV'];
run(ReflParamFile)

%%some constants:
%aB=5.291772e-11; %Bohr radius [m]
%ec=1.602177e-19; %e in [C]
aB=0.529177; % in [A]
ec2=14.399651; % in [eV A]

%Lindh. screening length
aL=(9*pi^2/128)^(1/3)*aB*(Zpr^(2/3)+Ztg^(2/3))^(-1/2);

%reduced energies; to be placed inside loop if not monoenergetic
eps_L=E0*(aL/(Zpr*Ztg*ec2))*(Mtg/(Mtg+Mpr));

%4-OTHER DEFINITIONS
%coefficients for the cosine-like distribution of theta angle
%f=r1*(cos(x))^n1+r2*(cos(x))^n2;
%applied to sputtered particles ; phi is random
%currently also applied to reflected particles: see 'particle_emission'
%for plans to use specular reflection for reflected particles 
r1=1.5;
r2=-1.0;
n1=0.1;
n2=0.9;


%%%%%%%%%%     RUN CASE     %%%%%%%%%%

%create output file and save input data
sfile=[outfolder,'/A',num2str(A),'_bx',num2str(bx),'_by',num2str(by),'_phi',num2str(phi*180/pi),'_delta',num2str(dlt*180/pi),'.mat'];
filename=sfile;

save(filename,'A','bx','by','NP','nsteps', 'phi','th', 'surfxmin','surfxmax','surfymin','surfymax','Tg','Pr');

run('zs_zp_intersect')

%run('surface_n_angles')

run('impact_output_and_distributions')

run('cellarea')

run('localflux')

run('plot_incoming_impacts')

run('erosion')

run('reflection_RN')

run('reflection_RE')

run('particle_emission')

run('plot_sputtering_and_reflection')

run('NSp_alongDiag')

movefile('*.png','outfile/')
