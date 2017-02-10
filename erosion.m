%%Erosion of Tg(target) by Pr(projectile)

%%display in this file from S20-S29
S20='erosion:';
disp(S20)

%%use distributions of #particles/unit cell and particle-surface angle
%%to calculate the local erosion using the Eckstein formula
%%eros(local)/eros(0) = Y_loc*Flux(local)/Y_flat*Flux(ave)
%%and the average flux (i.e. that for a flat surface) = #traj/Area

%%xg,yg(1:npoint),zg(1:npoints,1:npoints) = gridded version of the surface; 

%%xa, ya = surface-traject intersection point
%%pangle = local angle (that between surface and traject) at [xa,ya]
%%histrogram: Ncounts(1:npoints,1:npoints)


%0- load parametrization for Ecksteins fit-formula now in main file

w_eps=eps_L+0.1728*(eps_L)^(0.5)+0.008*(eps_L)^(0.1504);
sn=0.5*log(1+1.2288*eps_L)/w_eps;
Y_e=qeros*sn*((E0/Eth-1)^(mu_eros))/((lambda/w_eps)+(E0/Eth-1)^(mu_eros));


save(filename,'Y_e', '-append');   

%%1-initialize
Y_loc(1:npoints,1:npoints)=0.0;
Y_angle(1:npoints,1:npoints)=0.0;
cellNeros(1:npoints,1:npoints)=0.0;
cellFeros(1:npoints,1:npoints)=0.0;

%2-Calculate erosion of each cell

alpha0=pi-acos(sqrt(1/(1+E0/Esp)));

for i = 1:npoints
    for j = 1:npoints
    if (Ncounts(i,j)>0)
        %a-local angle
        locang=avecellangle(i,j);
        %locang=pi/4;%fix angle for tests:
    
        %%b-Calculate the erosion yield for the local angle and energy
        relang=locang/alpha0;
        term1=(cos((relang*pi/2)^c))^(-f);
        term2=exp(b*(1-1/cos((relang*pi/2)^c)));
        tf1=isreal(term1);
        tf2=isreal(term2);
    
        if (tf1==0)
            s21=[' term1 not real', num2str(term1), ' for loc_angle ', num2str(locang), ' rel_angle ', num2str(relang)];
            disp(s21)
        end
        if (tf2==0)
            s22=[' term2 not real', num2str(term2), ' for loc_angle ', num2str(locang),' rel_angle ', num2str(relang)];
            disp(s22)
        end
    
        Y_angle(i,j)=term1*term2;
    else %if not impacts in cell -> erosion yield = 0.0
        Y_angle(i,j)=0.0;
    end
        
    Y_loc(i,j)=Y_angle(i,j)*Y_e;
    %c-Local erosion
    cellNeros(i,j)=Y_loc(i,j)*Ncounts(i,j); %%N particles eroded, for Ncounts hitting the surface
    cellFeros(i,j)=Y_loc(i,j)*local_flux(i,j); %%N particles eroded, for Ncounts hitting the surface
    end
end

save(filename,'Y_angle', '-append');
save(filename,'Y_loc', '-append');
save(filename,'cellNeros', '-append');
save(filename,'cellFeros', '-append');