%use disp S4 - S6
S4='Calculate intersection of surface and trajectory ';
disp(S4)

%initialize; increase range of dim2 if adding output to 'local_angle' 
partglobal(1:NP,1:5)=0.0;
partlocal(1:NP,1:4)=0.0;


%%------------- HERE FOR-LOOP ON PARTICLES:
for p = 1:NP
        
    %S=['part  ',num2str(p)];
    %disp(S)

    %%z0=1.5*max(zs)
    [x,fxmin]=fminbnd(@(x) zs(x,x,Ax,fx,0,fy),surfxmin,surfxmax);
    [y,fymin]=fminbnd(@(y) zs(y,y,0,fx,Ay,fy),surfymin,surfymax);
    z0=-1.5*min(fxmin,fymin);
    
    %%different options for (x0,y0)
    %%a) random initial position (x0,y0), between (initmin,initmax)
    %rng(0,'twister');
    x0=0.5*(initxmax-initxmin)*(1-2*rand); %test print
    y0=0.5*(initymax-initymin)*(1-2*rand); %test print
    %%b) uniform grid of (x0,y0)
    %%NOT IMPLEMENTED YET
    
    %%so global (input) values of the particle are:
    partglobal(p,1)=x0;
    partglobal(p,2)=y0;
    partglobal(p,3)=z0;

%initiate angular distribution of th for each particle, based on a given analytical distribution
currentfolder = pwd;
distribFolder=[currentfolder, '/AngularDistributionFits'];

if strcmp(distr,'Curr85')
    %Curreli data set is in degrees, from theta 0 - 90
    %theta in Curreli data corresponds to delta in this code
    fitParamFile=[distribFolder,'/Curreli85'];
    run(fitParamFile) %load double gaussian parametrization
    run('vonNeumanAngle') %Use Von Neumann rejection method to generate values for delta
elseif strcmp(distr,'Boro85')
    %Borodkina data set is in degrees, from theta 0 - 90
    %theta in Borodkina data corresponds to delta in this code
    fitParamFile=[distribFolder,'/Borodkina85'];
    run(fitParamFile) %load double gaussian parametrization
    run('vonNeumanAngle') %Use Von Neumann rejection method to generate values for delta
elseif strcmp(distr,'Boro88')
    %Borodkina data set is in degrees, from theta 0 - 90
    %theta in Borodkina data corresponds to delta in this code
    fitParamFile=[distribFolder,'/Borodkina88'];
    run(fitParamFile) %load double gaussian parametrization
    run('vonNeumanAngle') %Use Von Neumann rejection method to generate values for delta
elseif strcmp(distr,'Boro89')
    %Borodkina data set is in degrees, from theta 0 - 90
    %theta in Borodkina data corresponds to delta in this code
    fitParamFile=[distribFolder,'/Borodkina89'];
    run(fitParamFile) %load double gaussian parametrization
    run('vonNeumanAngle') %Use Von Neumann rejection method to generate values for delta
else
    %fixed angle
    partglobal(p,4)=dlt;
    partglobal(p,5)=th;
end
    
    

    %%given an analytical surface and particle trajectory, find the intersection
    %%point and angle wrt the surface normal
    %%save all output (local particle's values) as:
    %first index (p) = particle index
    %%component 1:3 = impact point; 4 = angle wrt surface normal
    partlocal(p,:) = part_surf_local_angle(p,x0,y0,z0,phi,partglobal(p,5),Ax,fx,Ay,fy);

    %run('local_angle')
    
    
    if (mod(10*p/NP,1)==0)
        S5=['   ...', num2str(100*p/NP),'% done'];
        disp(S5)
    end
    
end
%%------------- LOOP DONE

save(filename,'partlocal','-append');
