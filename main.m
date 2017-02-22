%%JDC RESEARCH PLAN:

%%start with a trench surface in 2D based on sigmoid function, say
%%hyperbolic tangent function Z = (A*tanh(x-bx))
%%parameterized by A, bx, by.

%%launch particles from a fixed z0 (~1.5-2 max z-surf)
%%and (x0,y0) randomly chosen or forming an uniform grid

%%straight trajectories, forming an angle theta wrt z-axis, phi
%%initial set of parameter scan: 
%   A:f=0.1, 1.0, 10.0
%   theta = 0, 30, 60
%   phi = 0, (30,45,60?), 90

%%display lines in this file from S1-S9 (if needed)



%clear all values
clear variables
close all

%%%% output folder
currentfolder = pwd;
outfolder=[currentfolder,'/output_testtrench'];
mkdir(outfolder);
addpath(outfolder);




%%%%%%%%%%     DEFINE CASE     %%%%%%%%%%


%%1-define surface
A = 2;     %amplitude of trench height, in um
bx = 2;      %1/2 average trench length in x, in um
by = 1;      %1/2 average trench length in y, in um
surfxmin=-10; %surface x-min
surfxmax=10;  %surface x-max
surfymin=-10; %surface y-min
surfymax=10;  %surface y-max




%%2-define particles

%trajectories
phi=0.0;  %phi = angle wrt x-axis 0 < phi < pi/2
dlt=pi/3; %delta = angle wrt -z axis, (pointing to surface); 0<delta<pi/2
th=pi-dlt ; %theta =angle wrt +z axis ; pi/2 < theta < pi

%launching area
initxmin=-2; %x-min of initializing ('launching') particles
initxmax=2;  %x-max of initializing particles
initymin=-2; %y-min of initializing particles
initymax=2;  %y-max of initializing particles
z0=1; %start with specified height for now

%number of 'particles'
NP=60000;
nsteps=1000; %resolution: NP/nsteps = npoints = number of surface grids
npoints=floor(NP/nsteps);


%%3-define materials (for Eckstein's fit formula)
Tg='W'; %target material
Pr='N'; %projectile 
E0=100.0; %impact energy

%%load erosion and reflection parameters
Ecksteinfolder=[currentfolder,'/Eckstein'];
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

%coefficients for the cosine-like distribution of emitted particles' angle
%f=r1*(cos(x))^n1+r2*(cos(x)^n2;
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

movefile('*.png','outfile/')
