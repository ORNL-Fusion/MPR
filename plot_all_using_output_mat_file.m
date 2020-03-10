%%MPR Trench, for micro-trench erosion studies
%%A. Lasa & J. Coburn

%%This is a simplified version of the 'main' file, that only runs the
%%parts needed for plotting, 
%%Uses the output-mat file from previous runs

%clear all values
clear variables
close all

%load matlab file with all outputs:
filepath='~/work/GA_D3D/DIIID_exp_2019_CSkinner/MPR_output/v8_post-experiment_geometry_and_angle/PrD_TgSi_E50eV_A3um_fitChrobak_phi45_postExpGeom/output_files';
filename=[filepath,'/A3_bx15_by15_phi7.9577_delta85.9437.mat'];
load(filename)

%%%% output folder
currentfolder = pwd;
outfolder=[currentfolder,'/output_files'];
mkdir(outfolder);
addpath(outfolder);


%%%%%%%%%%     DEFINE VARIABLES THAT ARE NOT SAVED IN MAT FILE     %%%%%%%%%%

%%MATERIALS (for Eckstein's fit formula)
E0=100.0; %impact energy, eV

%%GEOMETRY
  
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

%%INITIAL CONDITIONS OF PARTICLES 
%launching area
initxmin=-20;   %-9; %x-min of initializing ('launching') particles
initxmax=16;     %x-max of initializing particles
initymin=-20;   %y-min of initializing particles
initymax=16;     %y-max of initializing particles
z0=0.01; %specified height of initializing particles

%number of 'particles' 
%resolution = number of surface grids;
npoints=floor(NP/nsteps);

%%OTHER DEFINITIONS
%coefficients for the cosine-like distribution of theta angle
%f=r1*(cos(x))^n1+r2*(cos(x))^n2;
r1=1.5;
r2=-1.0;
n1=0.1;
n2=0.9;


%%%%%%%%%%     RUN PLOTTING     %%%%%%%%%%

run('plot_incoming_impacts')

run('plot_sputtering_and_reflection')

movefile('*.png','outfile/')
