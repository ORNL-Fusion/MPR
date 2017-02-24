function y = emitted_part_surf_intersec(x0,y0,z0,phi,th,A,bx,by,initz0, dx, dy)


%use disp S7-S10

%%given an analytical surface and particle trajectory, find the
%%intersection point

%%zp=z0+vz*t -> t=(z-z0)/vp*cos(th)
%%xp=x0+vx*t = x0+ cos(phi)*sin(th)*((z-z0)/(cos(th)))=x0+(z-z0)*(cos(phi)*tan(th))
%%yp=y0+vy*t = y0+ sin(phi)*sin(th)*((z-z0)/(cos(th)))=y0+(z-z0)*(sin(phi)*tan(th));

%%look for intersecting (redeposition) point (xs, ys, q), given by zp=zs

%%a) if the surface is given analytically:
%%as the solution will be along the particle's trajectory, we can write
%%x as a function of y: x=x0+(ys-y0)/tan(phi)
%%and use it to pass only one variable to zs function

%redep=0 -> particle redeposited
%redep=1 -> zp=zs only at initial point (within error range of parameter 'erreps')
%redep=6 -> exitflag==-6 -> fzero did not find a solution
%redep=5 -> some other case in fzero

erreps=dx*sqrt(3)/2.; %error allowed when comparing the solution to initial position; ~cell size
initguess=8; %initial guess for fzero = x0+ initiguess*dx (or y0 + ...*dy)

if (tan(th)==0) %th=0 or pi
    %trajectory along z-axis -> NOT REDEPOSITED FOR A WELL DEFINED FUNCTION
    ys=y0;
    xs=x0;
    q=initz0;
    redep=1;
    %disp('redep=1; tan(th)==0') %test redep
    
elseif (th>0 && th<pi) %tan(th)!=0
    
    %%%%%%%%%%%%%%%%%% A) travelling towards +x, +y %%%%%%%%%%%%
    

    if (sin(phi) >= 0.0 && cos(phi) >= 0.0 )
        casename='A';   
        
        if (sin(phi)==0) %tan(phi)=0 -> ys=y0
            [xs,fval,exitflag,outinfo]=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,A,bx,by), x0+initguess*dx);%+initguess*dx (or dy) is an arbitrary distance
            ys=y0;
            
        elseif (cos(phi)==0)
            [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,A,bx,by), y0+initguess*dy); %as guess from initial point
            xs=x0;
            
        else
            [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0+initguess*dy); 
            xs=x0+(ys-y0)/tan(phi);            
        end
        
        if ~exist('ys','var')
            Sysexist=['ys undefined in case' , casename ];
            disp(Sysexist)
        end
        
        q=zs(xs,ys,A,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
       
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; xs=x0, ys=y0') %test redep
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; (xs+erreps)<x0') %test redep
                
            elseif ((xs<x0 || ys<y0)) %check answer is in the right direction:
                
                %if not, set to initial point and throw warning
                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' distance from initial position: ', num2str(dist)];
                disp(SredepWARN)
                Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Sdist)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;
                
                if (dist>10*initguess)
                    Sredepos=['RE-DEPOSITED! in case ', casename, ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos)
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist)
                end
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            disp('redep=6; exitflag==-6') %test redep
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag=!1, -6') %test redep
            disp(exitflag)
        end
        
        %%%%%%%%%%%%%%%%%% B) travelling towards -x, +y %%%%%%%%%%%%
        

    elseif (sin(phi) >= 0.0 && cos(phi) <= 0.0 )
        casename='B';
        
        if (sin(phi)==0)
            [xs,fval,exitflag,outinfo]=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,A,bx,by), x0+initguess*dx);
            ys=y0;
        elseif (cos(phi)==0)
            [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,A,bx,by), y0+initguess*dy);
            xs=x0;
        else
            [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0+initguess*dy);
            xs=x0+(ys-y0)/tan(phi);
        end
        
        if ~exist('ys','var')
            Sysexist=['ys undefined in case' , casename ];
            disp(Sysexist)
        end
        
        q=zs(xs,ys,A,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; xs=x0, ys=y0') %test redep
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; (xs+erreps)<x0') %test redep
                
            elseif (xs>x0 || ys<y0) %check answer is in the right direction:
                
                %if not, set to initial position and show warning:
                
                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' distance from initial position: ', num2str(dist)];
                disp(SredepWARN)
                Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Sdist)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;
                if (dist>10*initguess)
                    Sredepos=['RE-DEPOSITED! in case ', casename, ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos)
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist)
                end
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            disp('redep=6; exitflag==-6') %test redep
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag=!1, -6') %test redep
            disp(exitflag)
        end
        
        
        %%%%%%%%%%%%%%%%%% C) travelling towards -x, -y %%%%%%%%%%%%
        
        
    elseif (sin(phi) < 0.0 && cos(phi) <0.0 ) %sin(phi)==0, cos(phi)==0 already included in B
        casename='C';
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0-initguess*dy);
        xs=x0+(ys-y0)/tan(phi);
        
        if ~exist('ys','var')
            Sysexist=['ys undefined in case' , casename ];
            disp(Sysexist)
        end
        
        q=zs(xs,ys,A,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; xs=x0, ys=y0') %test redep
               
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; (ys+erreps)<y0') %test redep
                
            elseif (xs>x0 || ys>y0) %check answer is in the right direction:
                
                %if not, set to initial position and show warning:
                
                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' distance from initial position: ', num2str(dist)];
                disp(SredepWARN)
                Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Sdist)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;
                if (dist>10*initguess)
                    Sredepos=['RE-DEPOSITED! in case ', casename, ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos)
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist)
                end
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            disp('redep=6; exitflag==-6') %test redep
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag=!1, -6') %test redep
            disp(exitflag)
        end
        
        %%%%%%%%%%%%%%%%%% D) travelling towards +x, -y %%%%%%%%%%%%
        
    elseif (sin(phi) < 0.0 && cos(phi) >0.0 ) %%sin(phi)==0, , cos(phi)==0 already included in A
        casename='D';
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,bx,by), y0-initguess*dy);
        xs=x0+(ys-y0)/tan(phi);
        
        if ~exist('ys','var')
            Sysexist=['ys undefined in case' , casename ];
            disp(Sysexist)
        end
        
        q=zs(xs,ys,A,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; xs=x0, ys=y0') %test redep
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %disp('redep=1; (xs+erreps)<x0') %test redep
                    
            elseif (xs<x0 || ys>y0) %check answer is in the right direction:
                
                %if not, set to initial position and show warning:
                
                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' distance from initial position: ', num2str(dist)];
                disp(SredepWARN)
                Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Sdist)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;
                if (dist>10*initguess)
                    Sredepos=['RE-DEPOSITED! in case ', casename, ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos)
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist)
                end
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            disp('redep=6; exitflag==-6') %test redep
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag=!1, -6') %test redep
            disp(exitflag)
        end
    end
    
    
else %if (th<0 OR th>pi)
    q=initz0;
    xs=x0+(q-z0)*(cos(phi)*tan(th));
    ys=y0+(q-z0)*(sin(phi)*tan(th));
    redep=7;
    
end

%%%%%%%%%%%%% THERE MIGHT BE MORE EXCEPTION CASES %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% DEPENDING ON SURFACE %%%%%%%%%%%%%%%%%%

%{
...
%}


%%SKIP CALCULATION OF THE ANGLE AT THE INTERSECTING POINT


output(1)=xs;
output(2)=ys;
output(3)=q;
output(4)=redep;
%%output(4)=locang; %angle in rad
y = output;

end