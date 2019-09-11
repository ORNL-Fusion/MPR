%%for rough surfaces, a uniform grid formed as (xmax-xmin)/npoints will
%%lead to different areas for each cell: cellarea(i,j)


%a) if the surface is given analytically -> derivate analytically: dzs
%ns=[-dzsx(xs,Ax,fx),-dzsy(ys,Ay,fy),1];

%%initialize
cell_area(1:npoints,1:npoints)=0.0;
normal_cell_area(1:npoints,1:npoints)=0.0;
sg_theta(1:npoints,1:npoints)=0.0;
sg_phi(1:npoints,1:npoints)=0.0;
zk(1:npoints,1:npoints)=0.0;
nsxg(1:npoints,1:npoints)=0.0;
nsyg(1:npoints,1:npoints)=0.0;
nszg(1:npoints,1:npoints)=1.0;


%%xg, yg, zg = gridded version of the surface; 
dx=xg(2)-xg(1);
dy=yg(2)-yg(1);

run('surfacegrid')

for i = 1:npoints
    for j = 1:npoints
        xi=xg(i)+0.5*dx;
        yj=yg(j)+0.5*dy;
        zk(i,j)=zs(xi,yj,A,SX1,SX2,SY1,SY2,bx,by);
        nsxg(i,j)=-dzsx(xi,yj,A,SX1,SX2,bx,by);
        nsyg(i,j)=-dzsy(xi,yj,A,SY1,SY2,bx,by);
        %nszg(i,j)=1.0;
        sangles = surface_n_angles(nsxg(i,j),nsyg(i,j),nszg(i,j));
        sg_theta(i,j)=sangles(1);
        sg_phi(i,j)=sangles(2);
        normal_cell_area(i,j)=1/(cos(sg_theta(i,j)*sin(sg_phi(i,j)))*cos(sg_theta(i,j)*cos(sg_phi(i,j))));
        if (normal_cell_area(i,j)<=0) 
            Serr=['ERR: cell area =', num2str(cell_area(i,j)), 'giving value dx*dy'];
            disp(Serr)
            normal_cell_area(i,j)=dx*dy;
        end
        cell_area(i,j)=dx*dy/normal_cell_area(i,j);
    end
end

%cell_area(:,:)=dx*dy/(cos(sg_theta(:,:)*sin(sg_phi(:,:)))*cos(sg_theta(:,:)*cos(sg_phi(:,:))));

save(filename, 'nsxg', '-append');
save(filename, 'nsyg', '-append');
save(filename, 'sg_theta', '-append');
save(filename, 'sg_phi', '-append');
save(filename,'cell_area', '-append');
save(filename,'normal_cell_area', '-append');
