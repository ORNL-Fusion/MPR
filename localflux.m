
%flux for a flat surface: (# particles / (initial area)
Flux0=NP/((surfxmax-surfxmin)*(surfymax-surfymin));

%%initialize
local_flux(1:npoints,1:npoints)=0.0;
normal_local_flux(1:npoints,1:npoints)=0.0;


%local_flux=Ncounts/cell_area;

for i = 1:npoints 
    for j = 1:npoints 
        local_flux(i,j)=Ncounts(i,j)/cell_area(i,j);
        normal_local_flux(i,j)=local_flux(i,j)/Flux0;
    end
end



save(filename,'Flux0', '-append');
save(filename, 'local_flux', '-append');
save(filename, 'normal_local_flux', '-append');

