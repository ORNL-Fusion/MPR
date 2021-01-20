%%Reflection of Pr(projectile) from Tg(target)
%%display in this file from S20-S29
S21='reflection:';
disp(S21)

%%use distributions of #particles/unit cell and particle-surface angle
%%to calculate the local reflection using the Eckstein formula
%%refl(local)/refl(0) = RN_loc*Flux(local)/RN_flat*Flux(ave)
%%and the average flux (i.e. that for a flat surface) = #traj/Area

%%xg,yg(1:npoint),zg(1:npoints,1:npoints) = gridded version of the surface; 

%%xa, ya = surface-traject intersection point
%%pangle = local angle (that between surface and traject) at [xa,ya]
%%histrogram: Ncounts(1:npoints,1:npoints)

%%0-initialize
RN_loc(1:npoints,1:npoints)=0.0;
RN_angle(1:npoints,1:npoints)=0.0;
cellNrefl(1:npoints,1:npoints)=0.0;
cellFrefl(1:npoints,1:npoints)=0.0;


%1-load parametrization for Ecksteins fit-formula in main file


%2-Calculate the reflection yield for E0 energy
tf=strcmp(Pr,Tg); %compage target and projectile components

if (tf~=1) %NOT self-bombardment
    RN_0=b1*(eps_L)^b2/(1+b3*(eps_L)^b4);
else
    RN_0=exp(b1*(eps_L)^b2)/(1+exp(b3*(eps_L)^b4));
end

for i = 1:npoints
    for j = 1:npoints
        if (Ncounts(i,j)>0)

            %a-local angle -- not needed for now
            locang=avecellangle(i,j);

            %tanh works better at small angles; atan at large angles
            %but atan values range in (-pi/2, pi/2) --> normalize
            
            %NOTE: atan argument (c3*locang-c4) should be near zero not to diverge -> 
            %as c3>0 and c4<0 in tabulated values -> use +c4 in our implementation
            if (locang<pi/4)
                RN_angle(i,j)=c1+c2*tanh(c3*locang+c4);
            else
                RN_angle(i,j)=c1+c2*atan(c3*locang+c4)/(pi/2.);
            end
            
        else %if not impacts in cell -> reflection yield = 0.0
            RN_angle(i,j)=0.0;
        end
        RN_loc(i,j)=RN_angle(i,j)*RN_0;
        %c-Total reflection from cell
        cellNrefl(i,j)=RN_loc(i,j)*Ncounts(i,j); %%Nparticles reflected relative to those hitting the cell (Ncounts)
        cellFrefl(i,j)=RN_loc(i,j)*local_flux(i,j); %reflected flux of particles, including local cell area and Ncounts
        
    end
end

save(filename,'RN_0', '-append');   
save(filename,'RN_angle','-append');
save(filename,'RN_loc','-append');
save(filename,'cellNrefl', '-append');
save(filename,'cellFrefl', '-append');



