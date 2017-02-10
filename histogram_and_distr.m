%histrogram
Ncounts(1:npoints,1:npoints)=0;

di=(surfxmax-surfxmin)/(npoints-1);
dj=(surfymax-surfymin)/(npoints-1);

S11='binning impact points...';
disp(S11)


for p = 1:NP
    i=floor((xa(p)-surfxmin)/di)+1;
    j=floor((ya(p)-surfymin)/dj)+1;
    
    if (i<=0)
        erri=['WARNING, impact out of surface area; xs=', num2str(xa(p)), ' i=', num2str(i) ' ; setting i=1 and xs=surfxmin=', num2str(surfxmin) ];
        disp(erri)
        i=1;
        xa(p)=surfxmin;
    elseif (i>npoints)
        erri=['WARNING, impact out of surface area; xs=', num2str(xa(p)), ' i=', num2str(i) ' ; setting i=npoints and xs=surfxmax=' , num2str(surfxmax) ];
        i=npoints;
        xa(p)=surfxmax;
    end
    if (j<=0)
        errj=['WARNING, impact out of surface area; ys=', num2str(ya(p)), ' j=', num2str(j) ' ; setting j=1 and ys=surfymin=' , num2str(surfymin) ];
        disp(errj)
        j=1;
        ya(p)=surfymin;
    elseif (j>npoints)
        errj=['WARNING, impact out of surface area; ys=', num2str(ya(p)), ' j=', num2str(j) ' ; setting j=npoints and ys=surfymax=' , num2str(surfymax) ];
        disp(errj)
        j=npoints;
        ya(p)=surfymax;
    end
    
    Ncounts(i,j)=Ncounts(i,j)+1;
    
    if (mod(10*p/NP,1)==0)
        S4=['   ...',num2str(100*p/NP),'% done'];
        disp(S4)
    end
end

save(filename,'Ncounts', '-append');

%%cumulative distributions:

%histogram
ncdistr(1:max(max(Ncounts(:,:))))=0.0;
inc(1:max(max(Ncounts(:,:))))=0.0;
ncint=0.0;
for i = 1:npoints
    for j = 1:npoints
        if Ncounts(i,j)>0
            ncint=Ncounts(i,j);
            inc(ncint)=ncint;
            ncdistr(ncint)=ncdistr(ncint)+1;
        end
    end
end

save(filename,'inc','ncdistr', '-append');