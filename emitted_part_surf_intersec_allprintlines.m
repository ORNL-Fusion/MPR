function y = emitted_part_surf_intersec(x0,y0,z0,phi,th,Ax,fx,Ay,fy,initz0)


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
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0+3*pi/2.); %+2pi = as far as the solution can be from y0
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,Ax,fx,Ay,fy);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;     
                %S7=['     solution found, but only for xs (', num2str(xs), ') = x0(', num2str(x0) ') and ys (', num2str(ys), ') = y0(', num2str(y0) ')  ; assuming particle was NOT re-deposited'];
                %disp(S7)       
            elseif ((xs+erreps)<x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;
                %S7=['     solution found for xs (', num2str(xs), ') < x0(', num2str(x0) ') , but for cos(phi)>0 trajectory -> +x  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            elseif ((ys+erreps)<y0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
                %S7=['     solution found for  ys (', num2str(ys), ') < y0(', num2str(y0) ') , but for sin(phi)>0 trajectory -> +y  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            else
                redep=0;
                %{
                S7=['     particle emitted from x=', num2str(x0), ', y=', num2str(y0), ', z=', num2str(z0)];
                S8=['           RE-DEPOSITED in x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)];
                disp(S7)
                disp(S8)
                %}
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=4;
            %S7=['     exitflag 6: particle was NOT re-deposited: final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)'];
            %disp(S7)
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            %{
            S7=['     invalid value for fzero exitflag, ', num2str(exitflag), 'assuming particle was not re-deposited'];
            S8=['     Final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q), 'fval=',num2str(fval)];
            disp(S7)
            disp(S8)
            %}
            
        end
        
        %%%%%%%%%%%%%%%%%% B) travelling towards -x, +y %%%%%%%%%%%%
        
    elseif (sin(phi) > 0.0 && cos(phi) <0.0 )
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0+3*pi/2.); %+2pi = as far as the solution can be from y0
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,Ax,fx,Ay,fy);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %S7=['     solution found, but only for xs (', num2str(xs), ') = x0(', num2str(x0) ') and ys (', num2str(ys), ') = y0(', num2str(y0) ')  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            elseif ((xs+erreps)>x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;
                %S7=['     solution found for x0 (', num2str(x0), ') < xs(', num2str(xs) ') , but for cos(phi)<0 trajectory -> -x  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            elseif ((ys+erreps)<y0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
                %S7=['     solution found for  ys (', num2str(ys), ') < y0(', num2str(y0) ') , but for sin(phi)>0 trajectory -> +y  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            else
                redep=0;
                %{
                S7=['     particle emitted from x=', num2str(x0), ', y=', num2str(y0), ', z=', num2str(z0)];
                S8=['           RE-DEPOSITED in x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)];
                disp(S7)
                disp(S8)
                %}
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=4;
            %S7=['     exitflag 6: particle was NOT re-deposited: final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)'];
            %disp(S7)
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            %{
            S7=['     invalid value for fzero exitflag, ', num2str(exitflag), 'assuming particle was not re-deposited'];
            S8=['     Final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q), 'fval=',num2str(fval)];
            disp(S7)
            disp(S8)
            %}
            
        end
        
        
        %%%%%%%%%%%%%%%%%% C) travelling towards -x, -y %%%%%%%%%%%%
        
        
    elseif (sin(phi) < 0.0 && cos(phi) <0.0 )
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0-3*pi/2.); %+2pi = as far as the solution can be from y0
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,Ax,fx,Ay,fy);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %S7=['     solution found, but only for xs (', num2str(xs), ') = x0(', num2str(x0) ') and ys (', num2str(ys), ') = y0(', num2str(y0) ')  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            elseif ((xs+erreps)>x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;
                %S7=['     solution found for x0 (', num2str(x0), ') < xs(', num2str(xs) ') , but for cos(phi)<0 trajectory -> -x  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            elseif ((ys+erreps)>y0)  
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
                %S7=['     solution found for  y0 (', num2str(y0), ') < ys(', num2str(ys) ') , but for sin(phi)<0 trajectory -> -y  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            else
                redep=0;
                %{
                S7=['     particle emitted from x=', num2str(x0), ', y=', num2str(y0), ', z=', num2str(z0)];
                S8=['           RE-DEPOSITED in x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)];
                disp(S7)
                disp(S8)
                %}
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=4;
            %S7=['     exitflag 6: particle was NOT re-deposited: final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)'];
            %disp(S7)          
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            %{
            S7=['     invalid value for fzero exitflag, ', num2str(exitflag), 'assuming particle was not re-deposited'];
            S8=['     Final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q), 'fval=',num2str(fval)];
            disp(S7)
            disp(S8)
            %}           
        end
        
        %%%%%%%%%%%%%%%%%% D) travelling towards +x, -y %%%%%%%%%%%%
        
    elseif (sin(phi) < 0.0 && cos(phi) >0.0 ) %travelling towards +x, -y
        
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0-3*pi/2.); %+2pi = as far as the solution can be from y0
        xs=x0+(ys-y0)*tan(phi);
        q=zs(xs,ys,Ax,fx,Ay,fy);
        
        tfy=strcmp(num2str(ys),num2str(y0));
        tfx=strcmp(num2str(xs),num2str(x0));
        
        if (exitflag==1)
            
            if (tfy==1 && tfx==1)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=1;
                %S7=['     solution found, but only for xs (', num2str(xs), ') = x0(', num2str(x0) ') and ys (', num2str(ys), ') = y0(', num2str(y0) ')  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            elseif ((xs+erreps)<x0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=2;
                %S7=['     solution found for xs (', num2str(xs), ') < x0(', num2str(x0) ') , but for cos(phi)>0 trajectory -> +x  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            elseif ((ys+erreps)>y0)
                q=initz0;
                xs=x0+(q-z0)*(cos(phi)*tan(th));
                ys=y0+(q-z0)*(sin(phi)*tan(th));
                redep=3;
                %S7=['     solution found for  y0 (', num2str(y0), ') < ys(', num2str(ys) ') , but for sin(phi)<0 trajectory -> -y  ; assuming particle was NOT re-deposited'];
                %disp(S7)
            else
                redep=0;
                %{
                S7=['     particle emitted from x=', num2str(x0), ', y=', num2str(y0), ', z=', num2str(z0)];
                S8=['           RE-DEPOSITED in x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)];
                disp(S7)
                disp(S8)
                %}     
            end
            
        elseif (exitflag==-6)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=4;
            %S7=['     exitflag 6: particle was NOT re-deposited: final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q)'];
            %disp(S7)
            
        else
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=5;
            %{
            S7=['     invalid value for fzero exitflag, ', num2str(exitflag), 'assuming particle was not re-deposited'];
            S8=['     Final position: x=', num2str(xs), ', y=', num2str(ys), ', z=', num2str(q), 'fval=',num2str(fval)];
            disp(S7)
            disp(S8)
            %}
            
        end
    end
        
    
else %if (th<0 OR th>pi)
    q=initz0;
    xs=x0+(q-z0)*(cos(phi)*tan(th));    
    ys=y0+(q-z0)*(sin(phi)*tan(th));
    redep=6;
    %S7=['     invalid value for emitted particles theta, ', num2str(th*180/pi), 'assuming particle was not re-deposited'];
    %disp(S7)

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