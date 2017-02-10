%%display in this file from S20-S29
S21='Energy reflection:';
disp(S21)

%%average particle-surface angle of each cell to calculate the local 
%energy reflection yield using the Eckstein formula

%%xg,yg(1:npoint),zg(1:npoints,1:npoints) = gridded version of the surface;

%%0-initialize
RE_angle(1:npoints,1:npoints)=0.0;
RE_loc(1:npoints,1:npoints)=0.0;

%1-load parametrization for Ecksteins fit-formula in main file


%2-Calculate the energy reflection yield for E0 energy
tf=strcmp(Pr,Tg); %compage target and projectile components

if (tf~=1) %NOT self-bombardment
    RE_0=b1*(eps_L)^a2/(1+a3*(eps_L)^b4);
else
    RE_0=exp(a1*(eps_L)^a2)/(1+exp(a3*(eps_L)^a4));
end
       

for i = 1:npoints
    for j = 1:npoints
        if (Ncounts(i,j)>0)

            %a-local angle
            locang=avecellangle(i,j);

            %tanh works better at small angles; atan at large angles
            
            %NOTE: atan argument (d3*locang+d4) should be near zero not to diverge -> 
            %as d3>0 and d4<0 in tabulated values -> use -d4 in our implementation

            if (locang<pi/4)
                RE_angle(i,j)=d1+d2*tanh(d3*locang-d4);
            else
                RE_angle(i,j)=d1+d2*atan(d3*locang-d4);
            end
        
        else %if not impacts in cell -> reflection yield = 0.0
            RE_angle(i,j)=0.0;
        end
        RE_loc(i,j)=RE_angle(i,j)*RE_0;        
    end
end

save(filename,'RE_0', '-append');   
save(filename,'RE_angle','-append');
save(filename,'RE_loc','-append');




