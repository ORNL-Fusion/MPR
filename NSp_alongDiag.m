ip=0;
jp=0;
Spint(1:npoints)=0;
Lint=5.0;
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