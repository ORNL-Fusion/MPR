nfig=1;

%%display in this file from S30 on
S30='plotting incoming particles output...';
disp(S30)

%due to the way 3D surface plots auto-read axes from matrices, need to
%transpose each [x,y] matrix so that horizontal axis = x and vertical = y.
%Do within plot to not double variable names.

%1-background: contour plot
S31='      ... 1: surface';
disp(S31)

S31a='      ...   a: morphology';
disp(S31a)
figure(nfig)
nfig=nfig+1;

surf(xg,yg,transpose(zg));
hold on;
%set(p1,'EdgeColor','none');
colorbar;


title('1a: surface morphology')
xlabel('x')
ylabel('y')
zlabel('surface height')
axis equal;
hold off;

print('1a_morphology','-dpng')

%1b-cell area
S31b='      ...   b: normalized cell area (wrt flat surf)';
disp(S31b)
figure(nfig)
nfig=nfig+1;


contour(xg,yg,transpose(zg));
colorbar('off');
hold on;

surf(xg,yg,transpose(normal_cell_area));
colorbar;
caxis([1,max(max(normal_cell_area))]);

axis equal;
title('1b: cell area')
xlabel('x')
ylabel('y')
zlabel('area/S0')
hold off;

print('1b_cellarea','-dpng')


%1c-cells theta angle
S31c='      ...   c: theta angle of cells normal ';
disp(S31c)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(sg_theta*180/pi)); %in degrees
colorbar;
caxis([0.0,90.0]);
%caxis([-max(max(-180/pi*sg_theta)),max(max(180/pi*sg_theta))]);


axis equal;
title('1c: theta of cells normal')
xlabel('x')
ylabel('y')
zlabel('cell theta')
hold off;

print('1c_cellntheta','-dpng')


%1d-cells phi angle
S31d='      ...   d: phi angle of cells normal ';
disp(S31d)
figure(nfig)
nfig=nfig+1;


contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(sg_phi*180/pi)); %in degrees
colorbar;
%caxis([-90,90]);
caxis([-180,180]);
%caxis([-max(max(-180/pi*sg_phi)),max(max(180/pi*sg_phi))]);


axis equal;
title('1d: phi of cells normal')
xlabel('x')
ylabel('y')
zlabel('cell phi')
hold off;

print('1d_cellnphi','-dpng')


%2-angle
S32='      ... 2: angle';
disp(S32)

S32a='      ...   a: at each impact point';
disp(S32a)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

scatter(xa,ya,[],pangle*180/pi); 
colorbar;
caxis([0,90]); %180?
%caxis([-max(max(-180/pi*pangle)),max(max(180/pi*pangle))]);

axis equal;
title('2a: angle at each impact point')
xlabel('x')
ylabel('y')
zlabel('angle')
hold off;

print('2a_impactangle','-dpng')


%2b-average impact angle of each cell
S32b='      ...   b: average impact angle of each cell';
disp(S32b)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(avecellangle*180/pi));
colorbar;
caxis([0,90]);
%caxis([-max(max(-180/pi*avecellangle)),max(max(180/pi*avecellangle))]);


axis equal;
title('2b:average impact angle of each cell')
xlabel('x')
ylabel('y')
zlabel('ave angle')
hold off;

print('2b_aveimpactangle','-dpng')


%2c - angle profile (at x=5+/-0.5)
S32c=['      ...   c: profile along x (at y=',num2str(yprofile), '+/-', num2str(0.5), ')'];
disp(S32c)
figure(nfig)
nfig=nfig+1;

scatter(xangleprof,angleprofx*180/pi)
axis ([surfxmin surfxmax 0 90]); %180?
title('2c: profile of the angle along x')
xlabel('x')
ylabel('angle profile')
hold off;

print('2c_angleprofilex','-dpng')



%2d - angle profile (at y=5+/-0.5)
S32d=['      ...   d: profile along y (at x=',num2str(xprofile), '+/-', num2str(0.5), ')'];
disp(S32d)
figure(nfig)
nfig=nfig+1;

scatter(yangleprof,angleprofy*180/pi)
axis ([surfymin surfymax 0 90]);    %180?
title('2d: profile of the angle along y')
xlabel('y')
ylabel('angle profile')
hold off;

print('2d_angleprofiley','-dpng')


%2e-cumulative distribution of angles
S32e='      ...   e: cumulative distribution of angles';
disp(S32e)
figure(nfig)
nfig=nfig+1;

plot(ai,adistr);
title('2e: cumulative distribution of angles')
xlabel('number of cells w/ loc angle')
ylabel('local angle')

print('2e_anglehistogram','-dpng')



%3-flux
S33='      ... 3: flux';
disp(S33)

S33a='      ...   a: number of impacts per cell';
disp(S33a)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p3a=surf(xg,yg,transpose(Ncounts));
set(p3a,'EdgeColor','none');colorbar;
caxis([0,max(max(Ncounts))]);

axis equal;
title('3a: number of impacts per unit cell')
xlabel('x')
ylabel('y')
zlabel('# impacts')
hold off;

print('3a_Nimpact','-dpng')

%pause %troubleshooting

%10-local flux normalized
S33b='      ...   b: normalized local flux';
disp(S33b)
figure(nfig)
nfig=nfig+1;


contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p3b=surf(xg,yg,transpose(normal_local_flux));
set(p3b,'EdgeColor','none');colorbar;
caxis([0,max(max(normal_local_flux))]);

axis equal;
title('3b: normalized local flux (wrt flux at z0)')
xlabel('x')
ylabel('y')
zlabel('flux/flux0')
hold off;

print('3b_localflux','-dpng')


%3c: histogram along x (integrated over y)
S33c='      ...   c: profile of number of impacts along x (integrated along y, normalized)';
disp(S33c)
figure(nfig)
nfig=nfig+1;

plot(xg,profncountsx)
axis normal;
title('3c: profile of number of impacts along x')
xlabel('x')
ylabel('#impacts')
hold off;

print('3c_Nimpactprofilex','-dpng')



%3d: histogram along y (integrated over x)
S33d='      ...   d: profile of number of impacts along y (integrated along x, normalized)';
disp(S33d)
figure(nfig)
nfig=nfig+1;

plot(yg,profncountsy)
axis normal;
title('3d: profile of number of impacts along y')
xlabel('y')
ylabel('#impacts')
hold off;

print('3d_Nimpactprofiley','-dpng')


%3e-cumulative distribution of histogram
S33e='      ...   e: cumulative distribution of histogram';
disp(S33e)
figure(nfig)
nfig=nfig+1;

plot(inc, ncdistr);
title('3e: cumulative distribution of histogram')
xlabel('number of cells w/ #impacts')
ylabel('#impacts')

print('3e_Nimpacthistogram','-dpng')



S50=' done!';
disp(S50)
