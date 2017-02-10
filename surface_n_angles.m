function y = surface_n_angles(nsx,nsy,nsz)
%angle of surface normal wrt global coordinate system
%ns=[-dzsx(xs,Ax,fx),-dzsy(ys,Ay,fy),1];
%tan(th)=sqrt(nsx2+nsy2)/nsz; tan(phi)=nsy/nsx

%nsz=1.0, so cannot diverge (i.e., well defined surface)
surfacenth=atan(sqrt(nsx^2+nsy^2)/nsz);

if (surfacenth>pi/2.0 || surfacenth<0.0)
    dispErrth=['ERR: surface normal theta = ', num2str(surfacenth), 'is <0 or >pi/2; set to 0'];
    disp(dispErrth)
    surfacenth=0.0;
end

%phi can take any value: to avoid numerical issue, deal separately with 
%tan(phi)->inf: nsx=0=-(Ax/fx)*cos(x/fx) <-> x/fx=(n+1)pi/2 <-> perp to x-axis 

if (nsx==0.0 && nsy==0.0) %surfacenth = 0.0, phi not defined -> say 0.0
        surfacenphi=0.0;
    
elseif (nsy==0.0)
    
    if  (nsy>0) %cos(y/fy)<0 <-> pi/2<y/fy<3pi/2; point at +y
        surfacenphi=pi/2.0;
    elseif (nsy<0) %%cos(y/fy)>0 <-> -pi/2<y/fy<pi/2; point at -y
        surfacenphi=-pi/2.0;
    end
    
elseif (nsy==0.0) 
        
    if (nsx>0) %cos(x/fx)<0 <-> pi/2<x/fx<3pi/2; point at +x
        surfacenphi=0.0;
    elseif (nsx<0) %%cos(x/fx)>0 <-> -pi/2<x/fx<pi/2; point at -x
        surfacenphi=-pi;
    end
    
else
    tmp_phi=atan(nsy/nsx); %-pi/2<tmp_phi<pi/2
    if (nsx>0.0)
        surfacenphi=tmp_phi; %-pi/2<surfacenphi<pi/2
    elseif (nsx<0.0 && nsy>0.0) 
        surfacenphi=tmp_phi+pi; %-pi/2<tmp_phi<0 -> pi/2<surfacenphi<pi
    elseif (nsx<0.0 && nsy<0.0) 
        surfacenphi=tmp_phi-pi; %0<tmp_phi<pi/2 -> -pi<surfacenphi<-pi/2
    end
end


output(1)=surfacenth;
output(2)=surfacenphi;
y = output;

end