%use disp S4 - S6
S4='Calculate intersection of surface and trajectory ';
disp(S4)

%initialize; increase range of dim2 if adding output to 'local_angle' 
partglobal(1:NP,1:5)=0.0;
partlocal(1:NP,1:4)=0.0;

%%------------- HERE FOR-LOOP ON PARTICLES:
for p = 1:NP
        
    %S=['part  ',num2str(p)];
    %disp(S)
    
    %%different options for (x0,y0)
    %%a) random initial position (x0,y0), between (initmin,initmax)
    %rng(0,'twister');
    x0=0.5*(initxmax-initxmin)*(1-2*rand); %test print
    y0=0.5*(initymax-initymin)*(1-2*rand); %test print
    %%b) uniform grid of (x0,y0)
    %%NOT IMPLEMENTED YET
    
    %%so global (input) values of the particle are:
    partglobal(p,1)=x0;
    partglobal(p,2)=y0;
    partglobal(p,3)=z0;
    
    
    %initiate angular distribution of th for each particle, based on a given analytical distribution
    if (distr=='Curr85')
        %Curreli data set is in degrees, from theta 0 - 90
        %theta in Curreli data corresponds to delta in this code
        ga1=0.02985; %Double gaussian fit values for data
        gb1=75.01;
        gc1=6.697;
        ga2=0.03242;
        gb2=67.75;
        gc2=10.41;
        %Use Von Neumann rejection method to generate values for delta
        r = 0; %counter
        while (r~=1)
            dlttemp = 45*rand+45; %uniformly dist. random number for theta values between 45 and 90 degrees. 
                                        %Cutoff value is pi/4 or 45
                                        %degrees, to simplify search method
            ftemp = rand*0.08; %uniformly dist. random number between 0 and 0.08, the upper limit for distribution
            Curreli85dis = ga1*exp(-((dlttemp-gb1)/gc1)^2)+ga2*exp(-((dlttemp-gb2)/gc2)^2);
            if (ftemp<=Curreli85dis)
                partglobal(p,4)=dlttemp*pi/180;      %convert from degrees to radian, and set new value for delta for this code
                partglobal(p,5)=pi-(dlttemp*pi/180); %convert from degrees to radian, and set new value for theta for this code
                r=1; %end the Von Neumann rejection process
            else 
            end
        end
    elseif (distr=='Boro85')
        %Borodkina data set is in degrees, from theta 0 - 90
        %theta in Borodkina data corresponds to delta in this code
        ga1=0.02713; %Double gaussian fit values for data
        gb1=60.47;
        gc1=9.064;
        ga2=0.03234;
        gb2=71.59;
        gc2=10.61;
        %Use Von Neumann rejection method to generate values for delta
        r = 0; %counter
        while (r~=1)
            dlttemp = 45*rand+45; %uniformly dist. random number for theta values between 45 and 90 degrees. 
                                        %Cutoff value is pi/4 or 45
                                        %degrees, to simplify search method
            ftemp = rand*0.08; %uniformly dist. random number between 0 and 0.08, the upper limit for distribution
            Borod85dis = ga1*exp(-((dlttemp-gb1)/gc1)^2)+ga2*exp(-((dlttemp-gb2)/gc2)^2);
            if (ftemp<=Borod85dis)
                partglobal(p,4)=dlttemp*pi/180;      %convert from degrees to radian, and set new value for delta for this code
                partglobal(p,5)=pi-(dlttemp*pi/180); %convert from degrees to radian, and set new value for theta for this code
                r=1; %end the Von Neumann rejection process
            else 
            end
        end
    elseif (distr=='Boro88')
        %Borodkina data set is in degrees, from theta 0 - 90
        %theta in Borodkina data corresponds to delta in this code
        ga1=-0.05643; %Double gaussian fit values for data
        gb1=69.05;
        gc1=8.433;
        ga2=0.1111;
        gb2=70.2;
        gc2=9.454;
        %Use Von Neumann rejection method to generate values for delta
        r = 0; %counter
        while (r~=1)
            dlttemp = 45*rand+45; %uniformly dist. random number for theta values between 45 and 90 degrees. 
                                        %Cutoff value is pi/4 or 45
                                        %degrees, to simplify search method
            ftemp = rand*0.08; %uniformly dist. random number between 0 and 0.08, the upper limit for distribution
            Borod88dis = ga1*exp(-((dlttemp-gb1)/gc1)^2)+ga2*exp(-((dlttemp-gb2)/gc2)^2);
            if (ftemp<=Borod88dis)
                partglobal(p,4)=dlttemp*pi/180;      %convert from degrees to radian, and set new value for delta for this code
                partglobal(p,5)=pi-(dlttemp*pi/180); %convert from degrees to radian, and set new value for theta for this code
                r=1; %end the Von Neumann rejection process
            else 
            end
        end
    elseif (distr=='Boro89')
        %Borodkina data set is in degrees, from theta 0 - 90
        %theta in Borodkina data corresponds to delta in this code
        ga1=0.0145; %Double gaussian fit values for data
        gb1=75.64;
        gc1=5.011;
        ga2=0.05906;
        gb2=74.59;
        gc2=8.594;
        %Use Von Neumann rejection method to generate values for delta
        r = 0; %counter
        while (r~=1)
            dlttemp = 45*rand+45; %uniformly dist. random number for theta values between 45 and 90 degrees. 
                                        %Cutoff value is pi/4 or 45
                                        %degrees, to simplify search method
            ftemp = rand*0.08; %uniformly dist. random number between 0 and 0.08, the upper limit for distribution
            Borod89dis = ga1*exp(-((dlttemp-gb1)/gc1)^2)+ga2*exp(-((dlttemp-gb2)/gc2)^2);
            if (ftemp<=Borod89dis)
                partglobal(p,4)=dlttemp*pi/180;      %convert from degrees to radian, and set new value for delta for this code
                partglobal(p,5)=pi-(dlttemp*pi/180); %convert from degrees to radian, and set new value for theta for this code
                r=1; %end the Von Neumann rejection process
            else 
            end
        end
    elseif (distr=='Chrobk')
        %Chrobak data set is in degrees, from theta 0 - 90
        %C. Chrobak et. al., Nucl. Fusion 58 (2018) 106019
        %polar angle in Chrobak's paper data corresponds to delta in this code
        ga1=0.0135; %quatruple gaussian fit values for data
        gb1=82.87;  %to capture the multiple peaks correctly
        gc1=0.8441;
        ga2=0.07816;
        gb2=80.57;
        gc2=3.323;
        ga3=0.06401; 
        gb3=85.37;  
        gc3=2.346;
        ga4=0.04371;
        gb4=75.2;
        gc4=3.404;
        %Use Von Neumann rejection method to generate values for delta
        r = 0; %counter
        while (r~=1)
            dlttemp = 45*rand+45; %uniformly dist. random number for theta values between 45 and 90 degrees. 
                                        %Cutoff value is pi/4 or 45
                                        %degrees, to simplify search method
            ftemp = rand*0.09; %uniformly dist. random number (0, 0.09), the upper limit for Chrobak's distribution
            ChrobakDis = ga1*exp(-((dlttemp-gb1)/gc1)^2)+ga2*exp(-((dlttemp-gb2)/gc2)^2)+ga3*exp(-((dlttemp-gb3)/gc3)^2)+ga4*exp(-((dlttemp-gb4)/gc4)^2);
            if (ftemp<=ChrobakDis)
                partglobal(p,4)=dlttemp*pi/180;      %convert from degrees to radian, and set new value for delta for this code
                partglobal(p,5)=pi-(dlttemp*pi/180); %convert from degrees to radian, and set new value for theta for this code
                r=1; %end the Von Neumann rejection process
            else 
            end
        end
    elseif (distr=='ChrobC')
        %Chrobak data set is in degrees, from theta 0 - 90
        %C. Chrobak et. al., Nucl. Fusion 58 (2018) 106019
        %ChrobC follows the Chrobak distribution for D, -5deg, based on statement
        %on page 8 "IADs for C were similar to thos of D within 5 deg"
        %polar angle in Chrobak's paper data corresponds to delta in this code
        ga1=-0.216 ; %quatruple gaussian fit values for data
        gb1=77.91;  %to capture the multiple peaks correctly
        gc1=1.889;
        ga2=0.1267;
        gb2=77.87 ;
        gc2=1.495;
        ga3=0.1684; 
        gb3=78.04;  
        gc3=3.465;
        ga4=0.05222;
        gb4=71.25;
        gc4=4.174;
        %Use Von Neumann rejection method to generate values for delta
        r = 0; %counter
        while (r~=1)
            dlttemp = 45*rand+45; %uniformly dist. random number for theta values between 45 and 90 degrees.
                                        %Cutoff value is pi/4 or 45
                                        %degrees, to simplify search method
            ftemp = rand*0.085; %uniformly dist. random number (0, 0.085), the upper limit for Chrobak's distribution
            ChrobakDisC = ga1*exp(-((dlttemp-gb1)/gc1)^2)+ga2*exp(-((dlttemp-gb2)/gc2)^2)+ga3*exp(-((dlttemp-gb3)/gc3)^2)+ga4*exp(-((dlttemp-gb4)/gc4)^2);
            if (ftemp<=ChrobakDisC)
                partglobal(p,4)=dlttemp*pi/180;      %convert from degrees to radian, and set new value for delta for this code
                partglobal(p,5)=pi-(dlttemp*pi/180); %convert from degrees to radian, and set new value for theta for this code
                r=1; %end the Von Neumann rejection process
            else 
            end
            
        end
    else
        partglobal(p,4)=dlt;
        partglobal(p,5)=th;
    end
            
    %%given an analytical surface and particle trajectory, find the intersection
    %%point and angle wrt the surface normal
    %%save all output (local particle's values) as:
    %first index (p) = particle index
    %%component 1:3 = impact point; 4 = angle wrt surface normal
    partlocal(p,:) = part_surf_local_angle(p,x0,y0,z0,phi,partglobal(p,5),A,SX1,RX1,BX1,CX1,SX2,RX2,BX2,CX2,SY1,RY1,BY1,CY1,SY2,RY2,BY2,CY2,bx,by);

    %run('local_angle')
    
    
    if (mod(10*p/NP,1)==0)
        S5=['   ...', num2str(100*p/NP),'% done'];
        disp(S5)
    end
    
end
%%------------- LOOP DONE

save(filename,'partlocal','-append');
