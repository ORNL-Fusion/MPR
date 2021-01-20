function y = thompson_distr(Ein, Esb, Mtg, Mpr, npart, npoints)

u(1,npart)=0.0; %array of npart random number, uniformly distributed (0,1)

output(1:npart)=0.0;
func(1:npoints)=0.0; %array of discretized function
fsum(1:(npoints))=0.0; %array of cumulative values of the function 

gamma=4*Mtg*Mpr/(Mtg+Mpr)^2;
Eth=Esb/(gamma*(1-gamma));
Emax=Ein*gamma*(1-gamma)-Esb;

Exi=linspace(0,Emax, npoints);

%initial values:
func(1)=2*Esb*Exi(1)*(Ein/Eth)*(Ein/Eth)/(((Exi(1)+Esb)^3)*((Ein/Eth-1)^2));
fsum(1)=0.0; %or func(1)?

%not sure what to do for self-sputtering: gamma=1 -> Eth=inf...
%and I cannot find the source for this formulation of the Th distrib
tf=strcmp(num2str(Mtg),num2str(Mpr)); %compage target and projectile components (masses)

if (tf~=1) %NOT self-bombardment (things diverge!)
    disp('components of different masses; proceed with calculating the Thompson distributio')
    disp('for the energy of sputtered particles (although not used for anything here)')
    for i=2:npoints
        func(i)=2*Esb*Exi(i)*(Ein/Eth)*(Ein/Eth)/(((Exi(i)+Esb)^3)*((Ein/Eth-1)^2));
        fsum(i)=fsum(i-1)+func(i); %fsum(0)=0.0!
    end

    %normalize fsum, just in case:
    for i=1:npoints
        fsum(i)=fsum(i)/fsum(npoints);
    end


    for j=1:npart
        u(j)=rand;
    %binary search
        kmin=1;
        kmax=npoints;
        kmid=floor(npoints/2.0);
        err=0;

        while kmin<kmax-1 %while there's a difference of 2 or more indexes
            %{
            Stest=[ 'energy binary search for particle ' , num2str(j), ' u(j)=', num2str(u(j))];
            Stest2=['kmin=', num2str(kmin), ' kmid=' , num2str(kmid), ' kmax=' , num2str(kmax), ' fsum(kmin)=', num2str(fsum(kmin)), ' fsum(kmid)=' , num2str(fsum(kmin)), ' fsum(kmax)=' , num2str(fsum(kmax)) ];
            disp(Stest)
            disp(Stest2)
            %}
            if (u(j)>fsum(kmid) && u(j)<=fsum(kmax))
                kmin=kmid;
                kmid=floor((kmax+kmin)/2.0);
            elseif (u(j)<fsum(kmid) && u(j)>=fsum(kmin))
                kmax=kmid;
                kmid=floor((kmax+kmin)/2.0);
            else
                Serror=['error in binary search; u(j)=', num2str(u(j)) , ' fsum(kmin)=', num2str(fsum(kmin)), ' kmin=', num2str(kmin), ' fsum(kmid)=' , num2str(fsum(kmin)), ' kmid=' , num2str(kmin),' fsum(kmax)=' , num2str(fsum(kmax)),' kmax=' , num2str(kmax) ];
                disp(Serror)  
                err=1;
            end
        end

        %{
        if (err==0)
            Sbin=['found k | F(k)<u(j)<F(k+1) ; where u(j)=', num2str(u(j)) , ' k=', num2str(kmid) , ' fsum(kmin)=' , num2str(fsum(kmin)), ' fsum(k)=', num2str(fsum(kmid)), ' fsum(kmax)=', num2str(fsum(kmax)),  ' kmin=' , num2str(kmin),' kmax=' , num2str(kmax)];
            Sxi=['therefore the energy is: x(k)=', num2str(Exi(kmid))];
            disp(Sbin)           
            disp(Sxi)        
        end
        %}

        output(j)=Exi(kmid);
    end

else %if self sputtering (tf=1)
    %energy of sputtered particles = surface binding energy (made up!!)
    %this is currently not used anyway
    output(1:npart)=Esb;
    %warning moved to 'particle_emission', to avoid printing for each cell
	%disp('WARNING: you are running a case of self-sputtering/reflection (M1=M2)')
    %disp('         As implemented, the Thompson distrib would diverge -> fixed output E = Esb ')
    %disp('         This is completely made up, but it isnt used for anything, anyway')
    
end

y=output;
end