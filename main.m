%% MPR Main, for studying erosion of rough surfaces
%% A. Lasa & J. Coburn

%%start with a sinusoidal surface in 2D (aligned valleys or mountains);
%%parametrized by ratio Ax:fx, Ay:fy (amplitude):(1/frequency)

%%launch particles from a fixed z0 (~1.5-2 max z-surf)
%%and (x0,y0) randomly chosen or forming an uniform grid

%% straight trajectories, specified angular distribution & energy
%%angle theta is wrt z-axis, phi wrt x-axis
%%E.g., an initial set of parameter scan was performed: 
%%   A:f=0.1, 1.0, 10.0
%%   theta = 0, 30, 60
%%   phi = 0, (30,45,60?), 90
%% theta can also be given by a distribution (e.g., from PIC)

%%display lines in this file from S1-S9 (if needed)

%clear all values
clear variables
close all

%%%% output folder
currentfolder = pwd;
outfolder=[currentfolder,'/output'];
mkdir(outfolder);
addpath(outfolder);



%%%%%%%%%%     DEFINE CASE     %%%%%%%%%%


%%1-define surface
Ay=0.006;      %oscilation amplitude in y
fy=1.0;      %oscilation frequency in y
Ax=0.1*Ay;      %oscilation amplitude in x
fx=1.0;      %oscilation frequency in x 

surfxmin=-5.0*pi; %surface x-min
surfxmax=5.0*pi;  %surface x-max
surfymin=-5.0*pi; %surface y-min
surfymax=5.0*pi;  %surface y-max



%%2-define particles

%trajectories
phi=pi/4.0;    %phi = angle wrt x-axis 0 < phi < pi/2
dlt=0.0; %delta = angle wrt -z axis, (pointing to surface); 0<delta<pi/2
th=pi-dlt ; %theta =angle wrt +z axis ; pi/2 < theta < pi
distr='Boro88'; %distr = Curr85, Boro85, Boro88, Boro89, or blank. 
% A different name (not implemented in zs_zp_intersect) or blank uses dlt & th

%launching area
initxmin=-5.0*pi; %x-min of initializing ('launching') particles
initxmax=0.0;  %x-max of initializing particles
initymin=-5.0*pi; %y-min of initializing particles
initymax=0.0;  %y-max of initializing particles

%number of 'particles'
%nsteps = average #impacts per cell, in a flat surface: i.e., represents statistics
%npoints = number of surface grids: i.e., represents resolution
NP=200000;
nsteps=2000;
npoints=floor(NP/nsteps);




%%3-define materials (for Eckstein's fit formula)
Tg='W'; %target material
Pr='D'; %projectile 
E0=250.0; %impact energy

%%load erosion and reflection parameters
Ecksteinfolder=[currentfolder,'/Eckstein'];
addpath(Ecksteinfolder);

ErosParamFile=[Ecksteinfolder,'/Eros_','Tg_',Tg,'_Pr_',Pr, '_', num2str(E0), 'eV' ];
run(ErosParamFile)
%S0=['running ', ErosParamFile];
%disp(S0)

ReflNParamFile=[Ecksteinfolder,'/ReflN_','Tg_',Tg,'_Pr_',Pr, '_', num2str(E0), 'eV'];
run(ReflNParamFile)
%S1=['running ', ReflParamFile];
%disp(S1)

ReflEParamFile=[Ecksteinfolder,'/ReflE_','Tg_',Tg,'_Pr_',Pr, '_', num2str(E0), 'eV'];
run(ReflEParamFile)

%%some constants:
%aB=5.291772e-11; %Bohr radius [m]
%ec=1.602177e-19; %e in [C]
aB=0.529177; % in [A]
ec2=14.399651; % in [eV A]

%Lindh. screening length
aL=(9*pi^2/128)^(1/3)*aB*(Zpr^(2/3)+Ztg^(2/3))^(-1/2);

%reduced energies; to be placed inside loop if not monoenergetic
eps_L=E0*(aL/(Zpr*Ztg*ec2))*(Mtg/(Mtg+Mpr));

%coefficients for the cosine-like distribution of emitted particles' angle
%f=r1*(cos(x))^n1+r2*(cos(x)^n2;
r1=1.5;
r2=-1.0;
n1=0.1;
n2=0.9;



%%%%%%%%%%     RUN CASE     %%%%%%%%%%

%create output file and save input data
sfile=[outfolder,'/Ax',num2str(Ax),'_Ay',num2str(Ay),'_phi',num2str(phi*180/pi),'_delta',num2str(dlt*180/pi),'.mat'];
filename=sfile;

save(filename,'Ax','fx','Ay','fy','NP','nsteps', 'phi', 'th', 'surfxmin','surfxmax','surfymin','surfymax','Tg','Pr');

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

movefile('*.png',outfolder)
