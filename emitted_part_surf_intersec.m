function y = emitted_part_surf_intersec(x0,y0,z0,phi,th,A,bx,by,initz0)


%use disp S7-S10

%%given an analytical surface and particle trajectory, find the
%%intersection point

%%zp=z0+vz*t -> t=(z-z0)/vp*cos(th)
%%xp=x0+vx*t = x0+ cos(phi)*sin(th)*((z-z0)/(cos(th)))=x0+(z-z0)*(cos(phi)*tan(th))
%%yp=y0+vy*t = y0+ sin(phi)*sin(th)*((z-z0)/(cos(th)))=y0+(z-z0)*(sin(phi)*tan(th));

%%look for intersecting (redeposition) point (xs, ys, q), given by zp=zs

%%a) if the surface is given analytically:
%%as the solution will be along the particle's trajectory, we can write
%%x as a function of y: x=x0+(ys-y0)*tan(phi)
%%and use it to pass only one variable to zs function

erreps=0.01; %error allowed when comparing the solution to initial position; ~cell size


if (th>0 && th<=pi)
    
    
    %%%%%%%%%%%%%%%%%% A) travelling towards +x, +y %%%%%%%%%%%%
    
    
    if (sin(phi) > 0.0 && cos(phi) >0.0 )
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0+3*pi/2.); %+2pi = as far as the solution can be from y0
        if exist('ys','var')
            disp('ys undefined in case A (+x, +y)')
        end
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,A,bx,by);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        ys_half=(ys+y0)/2.0;
        xs_half=(xs+x0)/2.0;
        q_half=zs(xs_half,ys_half,A,bx,by);
        zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
            elseif (zp_half<q_half)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;                
            elseif ((xs+erreps)<x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
            elseif ((ys+erreps)<y0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=4;
            else
                redep=0;
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            
        end
        
        %%%%%%%%%%%%%%%%%% B) travelling towards -x, +y %%%%%%%%%%%%
        
    elseif (sin(phi) > 0.0 && cos(phi) <0.0 )
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0+3*pi/2.); %+2pi = as far as the solution can be from y0
        if exist('ys','var')
            disp('ys undefined in case B (-x, +y)')
        end
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,A,bx,by);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        ys_half=(ys+y0)/2.0;
        xs_half=(xs+x0)/2.0;
        q_half=zs(xs_half,ys_half,A,bx,by);
        zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
            elseif (zp_half<q_half)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;    
            elseif ((xs+erreps)>x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
            elseif ((ys+erreps)<y0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=4;
            else
                redep=0;
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            
        end
        
        
        %%%%%%%%%%%%%%%%%% C) travelling towards -x, -y %%%%%%%%%%%%
        
        
    elseif (sin(phi) < 0.0 && cos(phi) <0.0 )
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0-3*pi/2.); %+2pi = as far as the solution can be from y0
        if exist('ys','var')
            disp('ys undefined in case C (-x, -y)')
        end
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,A,bx,by);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        ys_half=(ys+y0)/2.0;
        xs_half=(xs+x0)/2.0;
        q_half=zs(xs_half,ys_half,A,bx,by);
        zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
            elseif (zp_half<q_half)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;    
            elseif ((xs+erreps)>x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
            elseif ((ys+erreps)>y0)  
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=4;
            else
                redep=0;
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;         
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;        
        end
        
        %%%%%%%%%%%%%%%%%% D) travelling towards +x, -y %%%%%%%%%%%%
        
    elseif (sin(phi) < 0.0 && cos(phi) >0.0 ) %travelling towards +x, -y
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0-3*pi/2.); %+2pi = as far as the solution can be from y0
        if exist('ys','var')
            disp('ys undefined in case D (+x, -y)')
        end
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,A,bx,by);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        ys_half=(ys+y0)/2.0;
        xs_half=(xs+x0)/2.0;
        q_half=zs(xs_half,ys_half,A,bx,by);
        zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
            elseif (zp_half<q_half)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;    
            elseif ((xs+erreps)<x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
            elseif ((ys+erreps)>y0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=4;
            else
                redep=0;    
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            
        end
    end
        
    
else %if (th<0 OR th>pi)
    q=initz0;
    xs=x0+(q-z0)*(cos(phi)*tan(th));    
    ys=y0+(q-z0)*(sin(phi)*tan(th));
    redep=7;

end

%%%%%%%%%%%%% THERE MIGHT BE MORE EXCEPTION CASES %%%%%%%%%%%%%
%%%%%%%%%%%%%% WITH PHI - THETA - AX- AY  %% TBD %%%%%%%%%%%%%%

%{

if (th>pi/2.0 && th<pi) %particle facing downwards -> directed to surface

     ys=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0);
     xs=x0+(ys-y0)*tan(phi);
     q=zs(xs,ys,Ax,fx,Ay,fy);

 elseif (th>0 && th<=pi/2.0)  %particle emitted upward
%}


%%SKIP CALCULATION OF THE ANGLE AT THE INTERSECTING POINT


output(1)=xs;
output(2)=ys;
output(3)=q;
output(4)=redep;
%%output(4)=locang; %angle in rad
y = output;

end