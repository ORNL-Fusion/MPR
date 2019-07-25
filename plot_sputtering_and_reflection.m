%due to the way 3D surface plots auto-read axes from matrices, need to
%transpose each [x,y] matrix. Do within plot to not double variable names

%%display in this file from S30 on
S30='plotting reflection and sputtering...';
disp(S30)

%4-sputtering
S34='      ... 4: sputtering';
disp(S34)

%a-angular contribution to sputtering
S34a='      ...   a: angular contribution to sputtering (wrt peak in angle)';
disp(S34a)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(Y_angle));
colorbar;
caxis([0,1]);

axis equal;
title('4a: angular contribution to sputtering (wrt max)')
xlabel('x')
ylabel('y')
zlabel('sputt Y_angle')
hold off;

print('4a_Yangle','-dpng')


%4b-average sputtering yield of each cell
S34b='      ...   b: average sputtering yield of each cell';
disp(S34b)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p4b=surf(xg,yg,transpose(Y_loc));
set(p4b,'EdgeColor','none'); colorbar;
caxis([0,max(max(Y_loc))]);

axis equal;
title('4b: average sputtering yield of each cell')
xlabel('x')
ylabel('y')
zlabel('ave sputt Y')
hold off;

print('4b_avesputtering','-dpng')


%4c-total erosion of each cell, or Ncounts hitting the surface
S34c='      ...   c: total erosion of each cell (N particles)';
disp(S34c)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p4c=surf(xg,yg,transpose(cellNeros));
set(p4c,'EdgeColor','none');colorbar;
caxis([0,max(max(cellNeros))]);

axis equal;
title('4c: total erosion of each cell (N particles)')
xlabel('x')
ylabel('y')
zlabel('erosion (N particles)')
hold off;

print('4c_Nerosion','-dpng')


S34c_abs=['      ... gross erosion (sput yield * N impacts) = ', num2str(sum(sum(cellNeros)))];
disp(S34c_abs)

S34c_frac=['      ... average erosion yield (wrt #impacts) = ', num2str(sum(sum(cellNeros))/NP)];
disp(S34c_frac)

%4d-total erosion of each cell
S34d='      ...   c: total erosion (flux) of each cell';
disp(S34d)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(cellFeros));
colorbar;
caxis([0,max(max(cellFeros))]);

axis equal;
title('4d: total erosion (flux) of each cell')
xlabel('x')
ylabel('y')
zlabel('total erosion (flux)')
hold off;

print('4d_erosionflux','-dpng')



%5: reflection

S35='      ... 5: reflection';
disp(S35)

%12a-angular contribution of reflection
S35a='      ...   a: angular contribution to reflection (wrt peak in angle)';
disp(S35a)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(RN_angle));
colorbar;
caxis([0,max(max(RN_angle))]);


axis equal;
title('5a: angular contribution to reflection (wrt max)')
xlabel('x')
ylabel('y')
zlabel('refl yield, RN_angle')
hold off;

print('5a_Rangle','-dpng')



%5b-average reflection yield of each cell
S35b='      ...   b: average reflection yield of each cell';
disp(S35b)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(RN_loc));
colorbar;
caxis([0,1.0]); %caxis([0,max(max(RN_loc))]);

axis equal;
title('5b: average reflection yield of each cell')
xlabel('x')
ylabel('y')
zlabel('refl yield, RN')
hold off;

print('5b_avereflection','-dpng')


%5c-total reflection from each cell (N particles)
S35c='      ...   c: total reflection from each cell (N particles)';
disp(S35c)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p5c=surf(xg,yg,transpose(cellNrefl));
set(p5c,'EdgeColor','none');colorbar; colormap jet;
caxis([0,max(max(cellNrefl))]);

axis equal;
title('5c: total reflection from each cell (N particles)')
xlabel('x')
ylabel('y')
zlabel('total reflection (N particles)')
hold off;

print('5c_Nreflection','-dpng')

S35c_abs=['      ... reflection (refl yield * N impacts) = ', num2str(sum(sum(cellNrefl)))];
disp(S35c_abs)

S35c_frac=['      ... average reflection yiled (wrt #impacts) = ', num2str(sum(sum(cellNrefl))/NP)];
disp(S35c_frac)

%5d-total reflection from each cell (flux)
S35d='      ...   d: total particle reflection (flux) from each cell';
disp(S35d)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p5d=surf(xg,yg,transpose(cellFrefl));
set(p5d,'EdgeColor','none');colorbar; colormap jet;
caxis([0,max(max(cellFrefl))]);

axis equal;
title('5d: total particle reflection (flux) from each cell')
xlabel('x')
ylabel('y')
zlabel('total particle reflection (flux)')
hold off;

print('5d_reflectionflux','-dpng')



%6: emission by sputtering

S36='      ... 6: particles emitted by sputtering';
disp(S36)

%6a-number of particles emitted by sputtering
S36a='      ...   a: number of particles';
disp(S36a)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p6a=surf(xg,yg,transpose(nSppart));
set(p6a,'EdgeColor','none');colorbar; colormap jet;
caxis([0,max(max(nSppart))]);

axis equal;
title('6a: number of particles emitted by sputtering')
xlabel('x')
ylabel('y')
zlabel('number of particles ')
hold off;

print('6a_Nemittsputt','-dpng')


%6b-energy of particles emitted by sputtering
S36b='      ...   b: average energy out (unused)';
disp(S36b)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(sum(SpEout,3)/size(SpEout,3)));
colorbar;
caxis([0,max(max(max(SpEout)))]);

axis equal;
title('6b: average energy of particles emitted by sputtering')
xlabel('x')
ylabel('y')
zlabel('energy ')
hold off;

print('6b_aveenergyemittsputt','-dpng')


%6c-theta angle of particles emitted by sputtering
S36c='      ...   c: Average Theta out';
disp(S36c)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveSpThout));
colorbar;
caxis([0,180]);
%caxis([-max(max(-180/pi*aveSpThout)),max(max(180/pi*aveSpThout))]);


axis equal;
title('6c: average theta of particles emitted by sputtering')
xlabel('x')
ylabel('y')
zlabel('theta')
hold off;

print('6c_avethetaemittsputt','-dpng')


%6d-average theta angle of emitted particle, wrt surf normal

S36d='      ...   d: ave theta wrt surf normal';
disp(S36d)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveSpThout_loc)); 
colorbar;
caxis([0,90]); %180?
%caxis([-max(max(-180/pi*aveSpThout_loc)),max(max(180/pi*aveSpThout_loc))]);

axis equal;
title('6d: ave theta  wrt surf normal of particles emitted by sputtering')
xlabel('x')
ylabel('y')
zlabel('theta')
hold off;

print('6d_avethetaemittsputt_nsurf','-dpng')


%6e-phi angle of particles emitted by sputtering
S36e='      ...   e: Average Phi out';
disp(S36e)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveSpPhiout));
colorbar;
%caxis([-180,180]);
caxis([-max(max(-180/pi*aveSpPhiout)),max(max(180/pi*aveSpPhiout))]);

axis equal;
title('6e: average phi of particles emitted by sputtering')
xlabel('x')
ylabel('y')
zlabel('phi')
hold off;

print('6e_avephiemittsputt','-dpng')

%6f-average phi angle of emitted particle, wrt surf normal

S36f='      ...   f: phi wrt surf normal';
disp(S36f)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveSpPhiout_loc)); 
colorbar;
%caxis([-180,180]); %180?
caxis([-max(max(-180/pi*aveSpPhiout_loc)),max(max(180/pi*aveSpPhiout_loc))]);


axis equal;
title('6f: ave phi wrt surf normal of particles emitted by sputtering')
xlabel('x')
ylabel('y')
zlabel('phi')
hold off;

print('6f_avephiemittsputt_nsurf','-dpng')



%7: emission by reflection

S36='      ... 7: particles emitted by reflection';
disp(S36)

%7a-number of particles emitted by sputtering
S37a='      ...   a: number of particles';
disp(S37a)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p7a=surf(xg,yg,transpose(nRpart));
set(p7a,'EdgeColor','none');colorbar; colormap jet;
caxis([0,max(max(nRpart))]);

axis equal;
title('7a: number of particles emitted by reflection')
xlabel('x')
ylabel('y')
zlabel('number of particles ')
hold off;

print('7a_Nemittrefl','-dpng')


%7b-energy of particles emitted by reflection
S37b='      ...   b: energy out (unused)';
disp(S37b)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(REout));
colorbar;
caxis([0,max(max(REout))]);

axis equal;
title('7b: energy of particles emitted by reflection')
xlabel('x')
ylabel('y')
zlabel('energy')
hold off;

print('7b_aveenergyemittrefl','-dpng')

%7c-theta of particles emitted by reflection
S37c='      ...   c: Average Theta out';
disp(S37c)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveRThout));
colorbar;
caxis([0,180]);
%caxis([-max(max(-180/pi*aveRThout)),max(max(180/pi*aveRThout))]);


axis equal;
title('7c: average theta of particles emitted by reflection')
xlabel('x')
ylabel('y')
zlabel('average theta')
hold off;

print('7c_avethetaemittrefl','-dpng')


%7d-average theta angle of emitted particle, wrt surf normal

S37d='      ...   d: ave theta wrt surf normal';
disp(S37d)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveRThout_loc)); 
colorbar;
caxis([0,90]); %180?
%caxis([-max(max(-180/pi*aveRThout_loc)),max(max(180/pi*aveRThout_loc))]);



axis equal;
title('7d: ave theta wrt surf normal of particles emitted by reflection')
xlabel('x')
ylabel('y')
zlabel('theta')
hold off;

print('7d_avethetaemittrefl_nsurf','-dpng')

%7e-phi of particles emitted by reflection
S37e='      ...   e: Average Phi out';
disp(S37e)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveRPhiout));
colorbar;
caxis([-180,180]);

axis equal;
title('7e: average phi of particles emitted by reflection')
xlabel('x')
ylabel('y')
zlabel('average phi')
hold off;

print('7e_avephiemittrefl','-dpng')


%7f-mean/average phi angle of emitted particle, wrt surf normal

S37f='      ...   f: ave phi wrt surf normal';
disp(S37f)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

surf(xg,yg,transpose(180/pi*aveRPhiout_loc)); 
colorbar;
caxis([-180,180]); %180?
%caxis([-max(max(-180/pi*aveRThout_loc)),max(max(180/pi*aveRThout_loc))]);

axis equal;
title('7f: ave phi wrt surf normal of particles emitted by reflection')
xlabel('x')
ylabel('y')
zlabel('phi')
hold off;

print('7f_avephiemittrefl_nsurf','-dpng')



%8-redeposition
S38='      ... 8: redeposition';
disp(S38)

%8a-number of particles redeposited after sputtering
S38a='      ...   a: number of particles redeposited after sputtering';
disp(S38a)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p8a=surf(xg,yg,transpose(NredepSp));
set(p8a,'EdgeColor','none');colorbar; colormap jet;
caxis([0,max(max(NredepSp))]);
%caxis([0,10]); %temporarily add specified z axis length, for direct comparison purposes

axis equal;
title('8a: number of redeposited particles (from sput)')
xlabel('x')
ylabel('y')
zlabel('# impacts')
hold off;

print('8a_redepossputt','-dpng')

S38a_abs=['      ... total redeposition (from sput) = ', num2str(sum(sum(NredepSp)))];
disp(S38a_abs)

S38a_frac=['      ... fraction redeposited (from sput) = ', num2str(sum(sum(NredepSp))/sum(sum(nSppart)))];
disp(S38a_frac)


%8b-number of particles redeposited after reflection
S38b='      ...   b: number of particles redeposited after reflection';
disp(S38b)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p8b=surf(xg,yg,transpose(NredepR));
set(p8b,'EdgeColor','none');colorbar; colormap jet;
caxis([0,max(max(NredepR))]);

axis equal;
title('8b: number of redeposited particles (from refl)')
xlabel('x')
ylabel('y')
zlabel('# impacts')
hold off;

print('8b_redeposrefl','-dpng')

S38b_abs=['      ... total redeposition (from refl) = ', num2str(sum(sum(NredepR)))];
disp(S38b_abs)

S38b_frac=['      ... fraction redeposited (from refl) = ', num2str(sum(sum(NredepR))/sum(sum(nRpart)))];
disp(S38b_frac)


%8c-total number of particles redeposited
S38c='      ...   c: total number of particles redeposited';
disp(S38c)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p8c=surf(xg,yg,transpose(NredepR+NredepSp));
set(p8c,'EdgeColor','none');colorbar; colormap jet;
caxis([0,max(max(NredepR+NredepSp))]);

axis equal;
title('8c: total number of redeposited particles')
xlabel('x')
ylabel('y')
zlabel('# impacts')
hold off;

print('8c_totalNredepos','-dpng')


%8d-net erosion (sput - redep by sput)
S38d='      ...   d: net erosion (sput - redep by sput)';
disp(S38d)
figure(nfig)
nfig=nfig+1;

contour(xg,yg,transpose(zg),1); %PLOT?
colorbar('off');
hold on;

p8d=surf(xg,yg,transpose(nSppart-NredepSp));
set(p8d,'EdgeColor','none');colorbar; colormap jet;
caxis([-max(max(-nSppart+NredepSp)),max(max(nSppart-NredepSp))]);

axis equal;
title('8d: net erosion (sput - redep by sput)')
xlabel('x')
ylabel('y')
zlabel('# net erosion')
hold off;

print('8d_neteros','-dpng')

S38d_abs=['      ... total net erosion (gross eros - redep by sput) = ', num2str(sum(sum(nSppart-NredepSp)))];
disp(S38d_abs)

S38d_frac=['      ... net erosion yield (wrt #impacts) = ', num2str(sum(sum(nSppart-NredepSp))/NP)];
disp(S38d_frac)

S38d_frac2=['      ... net erosion yield (gross eros yield * (1- redep fract) = ', num2str(sum(sum(cellNeros))/NP*(1-sum(sum(NredepSp))/sum(sum(nSppart))))];
disp(S38d_frac2)


S50=' done!';
disp(S50)
