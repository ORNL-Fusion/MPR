function y = part_surf_local_angle(p,x0,y0,z0,phi,th,Ax,fx,Ay,fy)

%use disp S7-S10

%%given an analytical surface and particle trajectory, find the intersection
%%point and angle wrt the surface normal 

locang=0.0; %initialize

%normalized velocity components = particle trajectory's angle
%vx=cos(phi)*sin(th); vy=sin(phi)*sin(th); vz=cos(th);

%%zp=z0+vz*t -> t=(z-z0)/cos(th)
%%xp=x0+vx*t = x0+ cos(phi)*sin(th)*((z-z0)/(cos(th)))=x0+(z-z0)*(cos(phi)*tan(th))
%%yp=y0+vy*t = y0+ sin(phi)*sin(th)*((z-z0)/(cos(th)))=y0+(z-z0)*(sin(phi)*tan(th));


%%look for intersecting point: zp=zs

%%a) if the surface is given analytically:
%%as the solution will be along the particle's trajectory, we can write
%%x as a function of y: x=x0+(ys-y0)/tan(phi) 
%%and use it to pass only one variable to zs function

 if (th>pi/2.0 && th<pi) %particle facing downwards
     if (Ax==0) %back to 1D sinusoidal in y
         if (phi==0.) %intersection point = zs(y0)
            q=zs(x0,y0,Ax,fx,Ay,fy); 
         elseif (phi==pi/2.) %tan(phi)-> inf 
            ys=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,Ax,fx,Ay,fy), y0); 
            xs=x0;  
            q=zs(xs,ys,Ax,fx,Ay,fy); 
         else % 0<phi<pi/2
            ys=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0); 
            xs=x0+(ys-y0)/tan(phi);
            q=zs(xs,ys,Ax,fx,Ay,fy);
         end
    
     elseif (Ay==0)
         if (phi==pi/2.) %intersection point = zs(x0)
            q=zs(x0,y0,Ax,fx,Ay,fy);
         elseif (phi==0.) %y=y0, tan(phi)=0
            xs=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,Ax,fx,Ay,fy), x0); 
            ys=y0;
            q=zs(xs,ys,Ax,fx,Ay,fy);
         else % 0<phi<pi/2
            xs=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0+(x-x0)*tan(phi),Ax,fx,Ay,fy), x0); 
            ys=y0+(xs-x0)*tan(phi); 
            q=zs(xs,ys,Ax,fx,Ay,fy); 
         end
     
     else %% Ax,&& Ay > 0
       if (phi==0) %%y=y0, tan(phi)=0
          xs=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,Ax,fx,Ay,fy), x0); 
          ys=y0;
          q=zs(xs,ys,Ax,fx,Ay,fy); 
       elseif (phi==pi/2.0) %%xs=x0, tan(phi)->inf.
          ys=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,Ax,fx,Ay,fy), y0); 
          xs=x0; 
          q=zs(xs,ys,Ax,fx,Ay,fy);
       else
           
        ys=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0);         
        xs=x0+(ys-y0)/tan(phi); 
        q=zs(xs,ys,Ax,fx,Ay,fy); 
    
       end
     end
     
 elseif (th>=pi || th<=pi/2.0)  %this should NOT happen
        q=z0;        
end

%%calculate (x,y) at the intersecting point
xs=x0+(q-z0)*cos(phi)*tan(th); 
ys=y0+(q-z0)*sin(phi)*tan(th); 


%%find normal to the surface at intersection: x-component=0, z-comp=-1

%a) if the surface is given analytically -> derivate analytically: dzs
%ns=[-dzs/dx,-dzs/dy,1]; 
nsx=-dzsx(xs,Ax,fx);
nsy=-dzsy(ys,Ay,fy);
nsz=1.0;
ns=[nsx,nsy,nsz];
%np=[d(zp)/dx,d(zp)/dy,-1] direction opposite to incoming particles'
npx=-cos(phi)*sin(th);
npy=-sin(phi)*sin(th);
npz=-cos(th);
np=[npx,npy,npz];


%%b) if the surface is not given analytically: 
%%see previous versions, commented lines here and above (for plotting)


%%and angle wrt particles trajectory
cosa=dot(np,ns)/(norm(np)*norm(ns));
locang=acos(cosa);
if(locang<0.0)
    locang=-acos(cosa);
    S8=[' For particle ',num2str(p),' impacting in xs = ',num2str(xs), ' ys =',num2str(ys),' ys =',num2str(q)];
    S9=[' angle =', num2str(180/pi*locang), '< 0; -angle will be taken '];    
    disp(S7);
    disp(S8);
end
if (cosa<0.0)
    locang=pi-acos(cosa);
    %{
S4=['For particle ',num2str(p),' impacting in xs =',num2str(xs), ' ys = ',num2str(ys),' zs = ',num2str(q)];
    S5=['angle = ', num2str(180/pi*locang), ' > 90; pi-angle will be taken '];  
    disp(S9);
    disp(S10);
    %}    
end

%{
partlocal(1)=xs;
partlocal(2)=ys;
partlocal(3)=q;
partlocal(4)=locang;
%}

output(1)=xs;
output(2)=ys;
output(3)=q;
output(4)=locang; %angle in rad
y = output;

end



    