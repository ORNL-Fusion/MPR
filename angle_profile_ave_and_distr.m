%angle profile: search for points along y-axis, around x=xprofile 
xprofile=pi/4.0;
ip=0;
yangleprof(1)=0;
angleprofy(1)=0;
for p = 1:NP
    if (xa(p)<(xprofile+pi/4.0) && xa(p)>(xprofile-pi/4.0))
        ip=ip+1;
        yangleprof(ip)=ya(p);
        angleprofy(ip)=pangle(p); %in degrees
    end
end

yprofile=pi/4.0;
jp=0;
xangleprof(1)=0;
angleprofx(1)=0;
for p = 1:NP
    if (ya(p)<(yprofile+pi/4.0) && ya(p)>(yprofile-pi/4.0))
        jp=jp+1;
        xangleprof(jp)=xa(p);
        angleprofx(jp)=pangle(p); %in degrees;
    end
end

save(filename,'xprofile','yprofile', '-append');
save(filename,'yangleprof','angleprofy', '-append');
save(filename,'xangleprof','angleprofx', '-append');

%%cumulative distributions:

%angle
adistr(1:90)=0.0; %180?
ai(1:90)=0.0;
for aint = 1:90 %180?
    ai(aint) = aint;
end
for k = 1:NP
    aint=floor((pangle(k)*180/pi))+1; %ai in [1:pi] or [1:pi/2]
    adistr(aint)=adistr(aint)+1;
end

save(filename,'ai','aint','adistr', '-append');

%average angle of each cell

cumulangle(1:npoints,1:npoints)=0.0;
avecellangle(1:npoints,1:npoints)=0.0;

for p=1:NP
    i=floor((xa(p)-surfxmin)/di)+1;
    j=floor((ya(p)-surfymin)/dj)+1;
    cumulangle(i,j)=cumulangle(i,j)+pangle(p); %in rad;
end

for i=1:npoints
    for j=1:npoints 
        if (Ncounts(i,j)>0)
            avecellangle(i,j)=cumulangle(i,j)/Ncounts(i,j);
        else
            avecellangle(i,j)=0.0;
        end
    end
end

save(filename,'avecellangle', '-append');