function y = thompson_distr(Ein, Esb, Mtg, Mpr, npart, npoints)

u(1,npart)=0.0; %array of npart random number, uniformly distributed (0,1)
uN(1,npart)=0.0; %as u(), but in range (0,npoint-1)
k(1:npart)=0.0;
func(1:npoints)=0.0; %array of discretized function
fsum(1:(npoints))=0.0; %array of cumulative values of the function 

gamma=4*Mtg*Mpr/(Mtg+Mpr)^2;
Eth=Esb/(gamma*(1-gamma));
Emax=Ein*gamma*(1-gamma)-Esb;

xi=linspace(0,Emax, npoints);

%initial values:
func(1)=2*Esb*xi(1)*(Ein/Eth)*(Ein/Eth)/(((xi(1)+Esb)^3)*((Ein/Eth-1)^2));
fsum(1)=0.0; %or func(1)?


for i=2:npoints
    func(i)=2*Esb*xi(i)*(Ein/Eth)*(Ein/Eth)/(((xi(i)+Esb)^3)*((Ein/Eth-1)^2));
    fsum(i)=fsum(i-1)+func(i); %fsum(0)=0.0!
end

%normalize fsum, just in case:
for i=1:npoints
    fsum(i)=fsum(i)/fsum(npoints);
end


for j=1:npart
    u(j)=rand;
    uN(j)=u(j)*(npoint-1);
    k(j)=floor(uN(j))+1;
end

y=k;


end