%Use Von Neumann rejection method to generate values for delta
%it's a method easy to implement, but not necessarily very efficient. 
%subject to change
r = 0; %counter
while (r~=1)
    dlttemp = 45*rand+45; %uniformly dist. random number for theta values between 45 and 90 degrees.
    %Cutoff value is pi/4 or 45
    %degrees, to simplify search method
    ftemp = rand*0.08; %uniformly dist. random number between 0 and 0.08, the upper limit for distribution
    distrFit = ga1*exp(-((dlttemp-gb1)/gc1)^2)+ga2*exp(-((dlttemp-gb2)/gc2)^2);
    if (ftemp<=distrFit)
        partglobal(p,4)=dlttemp*pi/180;      %convert from degrees to radian, and set new value for delta for this code
        partglobal(p,5)=pi-(dlttemp*pi/180); %convert from degrees to radian, and set new value for theta for this code
        r=1; %end the Von Neumann rejection process
    else
    end
end