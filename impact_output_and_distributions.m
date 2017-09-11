%%display in this file from S10-S19
S10='distributions and output:';
disp(S10)

%%derive distributions: #particles/unit cell and particle-surface angle
%%as well as cumulative distributions of both 2D plots

%%Horizontal axis is x; vertical axis is y; axes need to match
%%right-hand-rule!

%%surface-traject intersection point
xa=transpose(partlocal(:,1));
ya=transpose(partlocal(:,2));
za=transpose(partlocal(:,3));
%%and local angle (that between surface and traject) at [xa,ya], in rad
pangle=transpose(partlocal(:,4));

run('surfacegrid')

save(filename,'xg','yg','zg','-append');
save(filename,'xa','ya','za','pangle','-append');

run('histogram_and_distr')

profncountsx=sum(Ncounts,2)/max(sum(Ncounts,2));
save(filename,'profncountsx', '-append');

profncountsy=sum(Ncounts)/max(sum(Ncounts));
save(filename,'profncountsy', '-append');

run('angle_profile_ave_and_distr')

%run('redeposition_output')