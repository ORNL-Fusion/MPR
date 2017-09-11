function y = emitted_part_surf_intersec(x0,y0,z0,phi,th,Ax,fx,Ay,fy,initz0,dx,dy)


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

erreps=dx*sqrt(3)/2.; %error allowed when comparing the solution to initial position; ~cell size
initguess=pi*(1/fx+1/fy)/2.0; %initial guess for fzero = x0+ initiguess (or y0 + nitiguess)
%+2pi = as far as the solution can be from y0

if (tan(th)==0) %th=0 or pi
    %trajectory along z-axis -> NOT REDEPOSITED FOR A WELL DEFINED FUNCTION
    ys=y0;
    xs=x0;
    q=initz0;
    redep=1;
    %disp('redep=1; tan(th)==0') %test redep
    
elseif (th>0 && th<pi) %tan(th)!=0
    
    if (sin(phi) >= 0.0 && cos(phi) >= 0.0 ) %A) travelling towards +x, +y
        casename='A';
    elseif (sin(phi) >= 0.0 && cos(phi) <=0.0 ) % B) travelling towards -x, +y
        casename='B';
    elseif (sin(phi) < 0.0 && cos(phi) <0.0 ) % C) travelling towards -x, -y
        casename='C';
    elseif (sin(phi) < 0.0 && cos(phi) >0.0 ) %travelling towards +x, -y
        casename='D';
    end
    
    
    if (sin(phi)==0) %along x-axis: tan(phi)=0 -> ys=y0
        [xs,fval,exitflag,outinfo]=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,Ax,fx,Ay,fy), x0+initguess*sign(cos(phi)));
        ys=y0;
        
    elseif (cos(phi)==0) %along y-axis: tan(phi)=inf -> xs=x0
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,Ax,fx,Ay,fy), y0+initguess*sign(sin(phi)));
        xs=x0;
        
    else
        [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0+initguess*sign(sin(phi)));
        xs=x0+(ys-y0)/tan(phi);
    end
    
    if ~exist('ys','var')
        Syexist=['ys undefined in case ' , casename ];
        disp(Syexist)
    end
    
    if ~exist('xs','var')
        Sxexist=['xs undefined in case ' , casename ];
        disp(Sxexist)
    end
    
    q=zs(xs,ys,Ax,fx,Ay,fy);
    dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
    
    ys_half=(ys+y0)/2.0;
    xs_half=(xs+x0)/2.0;
    q_half=zs(xs_half,ys_half,Ax,fx,Ay,fy);
    zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
    
    if (exitflag==1) %found a solution o fzero
        
        if (dist<erreps) %not really a new position
            %disp7=['   NO REDEP: not really a new position, in case ', casename ];
            %disp(disp7)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=1;
        elseif (zp_half<q_half) %half way, particle is under surface
            %disp8=['   NO REDEP: particle travelling under surface, in case ', casename ];
            %disp(disp8)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=2;
        elseif (sign(xs-(x0+erreps))~=sign(cos(phi) ) || sign(ys-(y0+erreps))~=sign(sin(phi))) %wrong direction
            disp9=['   NO REDEP: travelling in the wrong direction, in case ', casename ];
            disp(disp9)
            q=initz0;
            xs=x0+(q-z0)*(cos(phi)*tan(th));
            ys=y0+(q-z0)*(sin(phi)*tan(th));
            redep=3;
        else
            %disp('REDEPOSITED!')
            redep=0;
        end
        
    elseif (exitflag==-6) %not found ~ not redep
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
        
%     elseif (sin(phi) >= 0.0 && cos(phi) <=0.0 )
%         casename='B';
%         
%         if (sin(phi)==0) %along x-axis: tan(phi)=0 -> ys=y0
%             [xs,fval,exitflag,outinfo]=fzero(@(x) zpx(x,x0,z0,th,phi)-zs(x,y0,Ax,fx,Ay,fy), x0-initguess);
%             ys=y0;
%             
%         elseif (cos(phi)==0) %along y-axis: tan(phi)=inf -> xs=x0
%             [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0,y,Ax,fx,Ay,fy), y0+initguess);
%             xs=x0;
%             
%         else
%             [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0+initguess);
%             xs=x0+(ys-y0)/tan(phi);
%         end
%         
%         if ~exist('ys','var')
%             Syexist=['ys undefined in case' , casename ];
%             disp(Syexist)
%         end
%         
%         if ~exist('xs','var')
%             Sxexist=['xs undefined in case' , casename ];
%             disp(Sxexist)
%         end
%         
%         q=zs(xs,ys,Ax,fx,Ay,fy);
%         dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
%               
%         ys_half=(ys+y0)/2.0;
%         xs_half=(xs+x0)/2.0;
%         q_half=zs(xs_half,ys_half,Ax,fx,Ay,fy);
%         zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
%         
%         if (exitflag==1)
%             
%             if (dist<erreps)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=1;
%             elseif (zp_half<q_half)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=2;
%             elseif ((xs+erreps)>x0 || (ys+erreps)<y0)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=3;           
%             else
%                 redep=0;
%             end
%             
%         elseif (exitflag==-6)
%             q=initz0;
%             xs=x0+(q-z0)*(cos(phi)*tan(th));
%             ys=y0+(q-z0)*(sin(phi)*tan(th));
%             redep=6;
%             
%         else
%             q=initz0;
%             xs=x0+(q-z0)*(cos(phi)*tan(th));
%             ys=y0+(q-z0)*(sin(phi)*tan(th));
%             redep=5;
%             
%         end
%         
%         
%         %%%%%%%%%%%%%%%%%% C) travelling towards -x, -y %%%%%%%%%%%%
%         
%         
%     elseif (sin(phi) < 0.0 && cos(phi) <0.0 )
%         casename='C';
%         
%         [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0-initguess); %+2pi = as far as the solution can be from y0
%         xs=x0+(ys-y0)/tan(phi);
%         
%         if ~exist('ys','var')
%             Syexist=['ys undefined in case' , casename ];
%             disp(Syexist)
%         end
%         
%         if ~exist('xs','var')
%             Sxexist=['xs undefined in case' , casename ];
%             disp(Sxexist)
%         end
%         
%         q=zs(xs,ys,Ax,fx,Ay,fy);
%         dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
%         
%         ys_half=(ys+y0)/2.0;
%         xs_half=(xs+x0)/2.0;
%         q_half=zs(xs_half,ys_half,Ax,fx,Ay,fy);
%         zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
%         
%         
%         if (exitflag==1)
%             
%             if (dist<erreps)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=1;
%             elseif (zp_half<q_half)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=2;
%             elseif ((xs+erreps)>x0 || (ys+erreps)>y0)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=3; 
%             else
%                 redep=0;
%             end
%             
%         elseif (exitflag==-6)
%             q=initz0;
%             xs=x0+(q-z0)*(cos(phi)*tan(th));
%             ys=y0+(q-z0)*(sin(phi)*tan(th));
%             redep=6;
%         else
%             q=initz0;
%             xs=x0+(q-z0)*(cos(phi)*tan(th));
%             ys=y0+(q-z0)*(sin(phi)*tan(th));
%             redep=5;
%         end
%         
%         %%%%%%%%%%%%%%%%%% D) travelling towards +x, -y %%%%%%%%%%%%
%         
%     elseif (sin(phi) < 0.0 && cos(phi) >0.0 ) %travelling towards +x, -y
%         casename='D';
%         
%         [ys,fval,exitflag,outinfo]=fzero(@(y) zpy(y,y0,z0,th,phi)-zs(x0+(y-y0)/tan(phi),y,Ax,fx,Ay,fy), y0-initguess); %+2pi = as far as the solution can be from y0
%         xs=x0+(ys-y0)/tan(phi);
%         
%         if ~exist('ys','var')
%             Syexist=['ys undefined in case' , casename ];
%             disp(Syexist)
%         end
%         
%         if ~exist('xs','var')
%             Sxexist=['xs undefined in case' , casename ];
%             disp(Sxexist)
%         end
%         
%         q=zs(xs,ys,Ax,fx,Ay,fy);
%         dist=sqrt((xs-x0)^2+(ys-y0)^2+(q-z0)^2);
%         
%         ys_half=(ys+y0)/2.0;
%         xs_half=(xs+x0)/2.0;
%         q_half=zs(xs_half,ys_half,Ax,fx,Ay,fy);
%         zp_half=z0+(ys_half-y0)/(tan(th)*sin(phi));
%         
%         if (exitflag==1)
%             
%             if (dist<erreps)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=1;
%             elseif (zp_half<q_half)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=2;
%             elseif ((xs+erreps)<x0 || (ys+erreps)>y0)
%                 q=initz0;
%                 xs=x0+(q-z0)*(cos(phi)*tan(th));
%                 ys=y0+(q-z0)*(sin(phi)*tan(th));
%                 redep=3;
%             else
%                 redep=0;
%             end
%             
%         elseif (exitflag==-6)
%             q=initz0;
%             xs=x0+(q-z0)*(cos(phi)*tan(th));
%             ys=y0+(q-z0)*(sin(phi)*tan(th));
%             redep=6;
%         else
%             q=initz0;
%             xs=x0+(q-z0)*(cos(phi)*tan(th));
%             ys=y0+(q-z0)*(sin(phi)*tan(th));
%             redep=5;
%             
%         end
%    end
    
    
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
     xs=x0+(ys-y0)/tan(phi);
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