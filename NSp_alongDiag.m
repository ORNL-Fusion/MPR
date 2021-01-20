%2 ways to run this file:
    % - at the end of a simulation, before clearing variables
        %comment out definition of xg, yg, dx, dy, nfig=0
    % - else, if all variables are not loaded, run 'main' until 'run case',
    %   load the simulation output mat-file & then run this file
    %keep the following lines uncommented:
nfig=0;
xg=linspace(surfxmin,surfxmax,npoints);
yg=linspace(surfymin,surfymax,npoints);
dx=xg(2)-xg(1);
dy=yg(2)-yg(1);

nfig=nfig+1;
ip=0;
jp=0;
Spint(1:npoints)=0;
Nint=npoints*Lint/(initxmax-initxmin);
for i = 1:npoints
    jmin=max(1,i-Nint);
    jmax=min(i+Nint,npoints);
    for j = jmin:jmax
        Spint(i)=Spint(i)+nSppart(i,j);
    end
end

figure(nfig+1)

x_int=1:1:npoints;
xrange=x_int*dx+initxmin;
%scatter(xrange,Spint)
plot(xrange,Spint)
axis ([surfxmin surfxmax 0 max(Spint)]); %180?
sTitle=['profile of number of sputtered particles along x=y, +/- ', num2str(Lint), 'um'];
title(sTitle)
xlabel('x')
ylabel('number of Sp particles')
hold off;

Sprint=['NSpParticles_alongDiagonal_',num2str(Lint),'um'];
print(Sprint,'-dpng')
%print('NSpParticles_alongDiagonal','-dpng')