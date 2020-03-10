%%xg, yg, zg = gridded version of the surface; 
xg=linspace(surfxmin,surfxmax,npoints);
yg=linspace(surfymin,surfymax,npoints);

zg(1:npoints,1:npoints)=0.0;
for i = 1:npoints
    for j = 1:npoints
        zg(i,j)=zs(xg(i),yg(j),A,SX1,RX1,BX1,CX1,SX2,RX2,BX2,CX2,SY1,RY1,BY1,CY1,SY2,RY2,BY2,CY2,bx,by);
    end
end