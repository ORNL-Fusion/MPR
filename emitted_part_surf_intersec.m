function y = emitted_part_surf_intersec(x0,y0,z0,phi,th,phi_loc,th_loc,phi_s,th_s,A,S,bx,by,initz0, dx, dy,event)


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
%initguess=8; %initial guess for fzero = x0+ initiguess*dx (or y0 + ...*dy)
%initial guess for fzero = x0+ sqrt((1.5*bx)^2+(1.5*by)^2)
%if along x or y, initial guess for fzero = x0 + 1.5*bx 

%use fzero with interval: from x0 to 1.1*bx (1.1 to ensure covering entire trench -- or 1.5)
if (tan(th)==0) %th=0 or pi
    %trajectory along z-axis -> NOT REDEPOSITED FOR A WELL DEFINED FUNCTION
    ys=y0;
    xs=x0;
    q=initz0;
    redep=1;
    disp('redep=1; tan(th)==0') %test redep %%UNCOMMENTED
    
elseif (th>0 && th<pi) %tan(th)!=0
    
    %%%%%%%%%%%%%%%%%% A) travelling towards +x, +y %%%%%%%%%%%%

    if (sin(phi) >= 0.0 && cos(phi) >= 0.0 )
        casename='A';   
        
        if (sin(phi)==0) %tan(phi)=0 -> ys=y0
            %initguess=[x0,1.5*bx];
            %intrv=[x0+dx 2.0*bx]; %1.1;
            range = 2*bx-(x0+dx);
            n = floor(abs(range)/dx); % #steps
            for i=1:n
                l = (x0+dx)+(i-1)*dx;
                r = (x0+dx)+i*dx;
                try
                    [xs,fval,exitflag,outinfo]=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,A,S,bx,by),[l r]);
                catch
                    %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                    %Sintrv=['       for case ', casename, ' sin(phi)==0 ; in interval ', num2str(intrv)];
                    %disp(Sintrv);
                    redep=6;
                    exitflag=-6;
                    xs=x0;
                    ys=y0;

                end
                if (xs ~= x0) %exit loop if we found a solution 
                    %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for x0 ', num2str(x0)];
                    %disp(Ssolve)
                    break
                end
            end
             
            ys=y0;
            
        elseif (cos(phi)==0)
            %initguess=[y0,1.5*by];
            %intrv=[y0+dy 2.0*by]; %1.1;
            range = 2*by-(y0+dy);
            n = floor(abs(range)/dy); % #steps
            for i=1:n
                l = (y0+dy)+(i-1)*dy;
                r = (y0+dy)+i*dy;
                try
                    [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,A,S,bx,by), [l r]);
                catch
                    %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                    %Sintrv=['       for case ', casename, ' cos(phi)==0; in interval', num2str(intrv)];
                    %disp(Sintrv);
                    redep=6;
                    exitflag=-6;
                    xs=x0;
                    ys=y0;
                end
                if (ys ~= y0) %exit loop if we found a solution 
                    %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for y0 ', num2str(y0)];
                    %disp(Ssolve)
                    break
                end
            end
            
            xs=x0;
            
        else
            %initguess=[y0,1.5*by];
            %intrv=[y0+dy 2.0*by];%1.1;
            range = 2*by-(y0+dy);
            n = floor(abs(range)/dy); % #steps
            for i=1:n
                l = (y0+dy)+(i-1)*dy;
                r = (y0+dy)+i*dy;
                try
                    [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,S,bx,by), [l r]);
                catch
                    %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                    %Sintrv=['       for case ', casename, ' in interval ', num2str(intrv)];
                    %disp(Sintrv);
                    redep=6;
                    exitflag=-6;
                    xs=x0;
                    ys=y0;
                end
                if (ys ~= y0) %exit loop if we found a solution 
                    %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for y0 ', num2str(y0)];
                    %disp(Ssolve)
                    break
                end
            end

            xs=x0+(ys-y0)/tan(phi);            
        end
        
        if ~exist('ys','var')
            Syexist=['ys undefined in case' , casename ];
            disp(Syexist)
        end
        
        q=zs(xs,ys,A,S,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
       
        if (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            %disp('redep=6; exitflag==-6') %test redep COMMENTED
            
        elseif (exitflag==1)
            
            if (tfy==1 && tfx==1) %xs=x0 & ys=y0
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                Sredep=['redep=1; xs=x0, ys=y0 ; case ', casename]; %test redep %%UNCOMMENTED
                disp(Sredep)
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                disp('redep=1; (xs+erreps)<x0') %test redep %%UNCOMMENTED
                
            elseif ((xs<x0 || ys<y0)) %check answer is in the right direction:
                
                %if not, set to initial point and throw warning
                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' for event ' , event ];
                disp(SredepWARN)
                Sdist=['         Angles: phi_s=' , num2str(phi_s) , ' th_s=', num2str(th_s), '; phi_loc=', num2str(phi_loc), ' th_loc=', num2str(th_loc), '; phi=', num2str(phi), ' th=', num2str(th) ];
                disp(Sdist)
                Spos=['         Distance from initial position: ', num2str(dist) , ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Spos)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;
                
                %if (dist>abs(sqrt((2*bx)^2+(2*by)^2))) %outside the trench
                if (xs<-1.1*bx || xs>1.1*bx || ys<-1.1*by || ys>1.1*by )
                    Sredepos=['re-deposited in case ', casename, ' event ', event, ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos) %remove output for now %%UNCOMMENTED
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist) %remove output for now %%UNCOMMENTED
                    disp('      discard particle for now') 
                    redep=1;
                else
                    Sredep=['redep=0; SUCCESS! in xs=', num2str(xs), ' ys=', num2str(ys)];
                    %disp(Sredep) %test redep 
                end
            end
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag =! 1, -6') %test redep
            disp(exitflag)
        end
        
        %%%%%%%%%%%%%%%%%% B) travelling towards -x, +y %%%%%%%%%%%%
        

    elseif (sin(phi) >= 0.0 && cos(phi) <= 0.0 )
        casename='B';
        
        if (sin(phi)==0)
            %initguess=[x0,-1.5*bx];
            %intrv=[-2.0*bx x0-dx]; %-1.5;
            range = x0-dx+2*bx;
            n = floor(abs(range)/dx); % #steps
            for i=1:n
                l = (x0-dx)-i*dx;
                r = (x0-dx)-(i-1)*dx;
                try
                    [xs,fval,exitflag,outinfo]=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,A,S,bx,by), [l r]);
                catch
                    %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                    %Sintrv=['case ', casename, ' sin(phi)==0; in interval', num2str(intrv)];
                    %disp(Sintrv);
                    redep=6;
                    exitflag=-6;
                    xs=x0;
                    ys=y0;
                end
                if (xs ~= x0) %exit loop if we found a solution 
                    %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for x0 ', num2str(x0)];
                    %disp(Ssolve)
                    break
                end
            end
            
            ys=y0;
            
        elseif (cos(phi)==0)
            %initguess=[y0,1.5*by];
            %intrv=[y0+dy 2.0*by]; %1.1;
            range = 2*by-(y0+dy);
            n = floor(abs(range)/dy); % #steps
            for i=1:n
                l = (y0+dy)+(i-1)*dy;
                r = (y0+dy)+i*dy;
                try
                    [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,A,S,bx,by),[l r]);
                catch
                    %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                    %Sintrv=['for case ', casename, ' cos(phi)==0; in interval ', num2str(intrv)];
                    %disp(Sintrv);
                    redep=6;
                    exitflag=-6;
                    xs=x0;
                    ys=y0;
                end
                if (ys ~= y0) %exit loop if we found a solution 
                    %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for y0 ', num2str(y0)];
                    %disp(Ssolve)
                    break
                end
            end

            xs=x0;
        else
            %initguess=[y0,1.5*by];
            %intrv=[y0+dy 2.0*by]; %1.1;
            range = 2*by-(y0+dy);
            n = floor(abs(range)/dy); % #steps
            for i=1:n
                l = (y0+dy)+(i-1)*dy;
                r = (y0+dy)+i*dy;
                try
                    [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,S,bx,by),[l r]);
                catch
                    %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                    %Sintrv=['for case ', casename, ' in interval ', num2str(intrv)];
                    %disp(Sintrv);
                    redep=6;
                    exitflag=-6;
                    xs=x0;
                    ys=y0;
                end
                if (ys ~= y0) %exit loop if we found a solution 
                    %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for y0 ', num2str(y0)];
                    %disp(Ssolve)
                    break
                end
            end
            
            xs=x0+(ys-y0)/tan(phi);
        end
        
        if ~exist('ys','var')
            Syexist=['ys undefined in case' , casename ];
            disp(Syexist)
        end
        
        q=zs(xs,ys,A,S,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            %disp('redep=6; exitflag==-6') %test redep COMMENTED
            
        elseif (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                Sredep=['redep=1; xs=x0, ys=y0 ; case ', casename]; %test redep %%UNCOMMENTED
                disp(Sredep)
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                disp('redep=1; (xs+erreps)<x0') %test redep %%UNCOMMENTED
                
            elseif (xs>x0 || ys<y0) %check answer is in the right direction:
                
                %if not, set to initial position and show warning:
                
                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' event ' , event ];
                disp(SredepWARN)
                Sdist=['        Angles: phi_s=' , num2str(phi_s) , ' th_s=', num2str(th_s), '; phi_loc=', num2str(phi_loc), ' th_loc=', num2str(th_loc), '; phi=', num2str(phi), ' th=', num2str(th) ];
                disp(Sdist)
                Spos=['        Distance from initial position: ', num2str(dist) , ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Spos)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;

                %if (dist>abs(sqrt((2*bx)^2+(2*by)^2))) %outside the trench
                if (xs<-1.1*bx || xs>1.1*bx || ys<-1.1*by || ys>1.1*by )
                    Sredepos=['re-deposited in case ', casename, ' event ' , event , ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos) %remove output for now %%UNCOMMENTED
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist)  %remove output for now %%UNCOMMENTED
                    disp('      discard particle for now') 
                    redep=1;
                else
                    Sredep=['redep=0; SUCCESS! in xs=', num2str(xs), ' ys=', num2str(ys)];
                    %disp(Sredep) %test redep
                end
            end

        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag =! 1, -6') %test redep
            disp(exitflag)
        end
        
        
        %%%%%%%%%%%%%%%%%% C) travelling towards -x, -y %%%%%%%%%%%%
        
        
    elseif (sin(phi) < 0.0 && cos(phi) <0.0 ) %sin(phi)==0, cos(phi)==0 already included in B
        casename='C';
        
        %initguess=[y0,-1.5*by];
        %intrv=[-2.0*by y0-dy]; %-1.1;
        range = y0-dy+2*by;
        n = floor(abs(range)/dy); % #steps
        for i=1:n
            l = (y0-dy)-i*dy;
            r = (y0-dy)-(i-1)*dy;
            try
                [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,S,bx,by),[l r]);        
            catch
                %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                %Sintrv=['for case ', casename, ' sin(phi) & cos(phi) <0.0 ; in interval ', num2str(intrv)];
                %disp(Sintrv);
                redep=6;
                exitflag=-6;
                xs=x0;
                ys=y0;
            end
            if (ys ~= y0) %exit loop if we found a solution 
                %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for y0 ', num2str(y0)];
                %disp(Ssolve)
                break
            end
        end

        xs=x0+(ys-y0)/tan(phi);
        
        if ~exist('ys','var')
            Syexist=['ys undefined in case' , casename ];
            disp(Syexist)
        end
        
        q=zs(xs,ys,A,S,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            %disp('redep=6; exitflag==-6') %test redep COMMENTED
            
        elseif (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                Sredep=['redep=1; xs=x0, ys=y0 ; case ', casename]; %test redep %%UNCOMMENTED
                disp(Sredep)
               
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                disp('redep=1; (ys+erreps)<y0') %test redep %UNCOMMENTED
                
            elseif (xs>x0 || ys>y0) %check answer is in the right direction:
                
                %if not, set to initial position and show warning:

                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' event ' , event ];
                disp(SredepWARN)
                Sdist=['        Angles: phi_s=' , num2str(phi_s) , ' th_s=', num2str(th_s), '; phi_loc=', num2str(phi_loc), ' th_loc=', num2str(th_loc), '; phi=', num2str(phi), ' th=', num2str(th) ];
                disp(Sdist)
                Spos=['        Distance from initial position: ', num2str(dist) , ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Spos)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;
                
                %if (dist>abs(sqrt((2*bx)^2+(2*by)^2))) %outside the trench
                if (xs<-1.1*bx || xs>1.1*bx || ys<-1.1*by || ys>1.1*by )
                    Sredepos=['re-deposited in case ', casename, ' event ' , event , ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos) %remove output for now %%UNCOMMENTED
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist) %remove output for now %%UNCOMMENTED
                    disp('      discard particle for now')
                    redep=1;
                else
                    Sredep=['redep=0; SUCCESS! in xs=', num2str(xs), ' ys=', num2str(ys)];
                    %disp(Sredep) %test redep
                end
            end
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag =! 1, -6') %test redep
            disp(exitflag)
        end
        
        %%%%%%%%%%%%%%%%%% D) travelling towards +x, -y %%%%%%%%%%%%
        
    elseif (sin(phi) < 0.0 && cos(phi) >0.0 ) %%sin(phi)==0, , cos(phi)==0 already included in A
        casename='D';
        
        %initguess=[y0-1.5*by];
        %intrv=[-2.0*by y0-dy]; %-1.1;
        range = y0-dy+2*by;
        n = floor(abs(range)/dy); % #steps
        for i=1:n
            l = (y0-dy)-i*dy;
            r = (y0-dy)-(i-1)*dy;
            try
                [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,A,S,bx,by), [l r]);        
            catch
                %disp('no solution found by fzero (err); set redep=6 & exitflag=-6, xs=x0, ys=y0');
                %Sintrv=['for case ', casename, ' sin(phi) <0 & cos(phi) >0.0 ; in interval ', num2str(intrv)];
                %disp(Sintrv);
                redep=6;
                exitflag=-6;
                xs=x0;
                ys=y0;
            end
            if (ys ~= y0) %exit loop if we found a solution 
                %Ssolve=['found a solution after iteration ', num2str(i), ' interval [',num2str(l), ' , ',num2str(r),'] for y0 ', num2str(y0)];
                %disp(Ssolve)
                break
            end
        end

        
        %y0-1.5*by);
        xs=x0+(ys-y0)/tan(phi);
        
        if ~exist('ys','var')
            Syexist=['ys undefined in case' , casename ];
            disp(Syexist)
        end
        
        q=zs(xs,ys,A,S,bx,by);
        
        dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=6;
            %disp('redep=6; exitflag==-6') %test redep COMMENTED
            
        elseif (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                Sredep=['redep=1; xs=x0, ys=y0 ; case ', casename]; %test redep %%UNCOMMENTED
                disp(Sredep)
                
            %in case it found a value almost identical to initial point
            elseif (dist<erreps)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                disp('redep=1; (xs+erreps)<x0') %test redep %UNCOMMENTED
                    
            elseif (xs<x0 || ys>y0) %check answer is in the right direction:
                
                %if not, set to initial position and show warning:
                
                SredepWARN=['WARNING: Found solution in wrong direction, in case ', casename, ' event ' , event ];
                disp(SredepWARN)
                Sdist=['        Angles: phi_s=' , num2str(phi_s) , ' th_s=', num2str(th_s), '; phi_loc=', num2str(phi_loc), ' th_loc=', num2str(th_loc), '; phi=', num2str(phi), ' th=', num2str(th) ];
                disp(Sdist)
                Spos=['        Distance from initial position: ', num2str(dist) , ' Setting from solution: xs=', num2str(xs), ' ys=', num2str(ys) , ' to initial position:  x0=', num2str(x0), ' y0=', num2str(y0) ]; %test redep
                disp(Spos)
                
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                
            else
                redep=0;
                
                %if (dist>abs(sqrt((2*bx)^2+(2*by)^2))) %outside the trench
                if (xs<-1.1*bx || xs>1.1*bx || ys<-1.1*by || ys>1.1*by )
                    Sredepos=['re-deposited in case ', casename, ' event ' , event , ' but at distance from initial position: ', num2str(dist)];
                    disp(Sredepos) %remove output for now %%UNCOMMENTED
                    Sdist=['        For phi=', num2str(phi), ' th=', num2str(th), ' x0=', num2str(x0),' xs=', num2str(xs), ' y0=', num2str(y0),  ' ys=', num2str(ys) ]; %test redep
                    disp(Sdist) %remove output for now %%UNCOMMENTED
                    disp('      discard particle for now')
                    redep=1;
                else
                    Sredep=['redep=0; SUCCESS! in xs=', num2str(xs), ' ys=', num2str(ys)];
                    %disp(Sredep) %test redep 
                end
            end
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            disp('redep=5; exitflag =! 1, -6') %test redep
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