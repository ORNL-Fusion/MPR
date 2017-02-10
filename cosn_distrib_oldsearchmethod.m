%generate random numbers with an analytically given distribution:
function y = cosn_distrib(a1, n1, a2, n2, npart, npoints)

u(1,npart)=0.0; %array of npart random number, uniformly distributed (0,1)
uN(1,npart)=0.0; %as u(), but in range (0,npoints-1)
k(1:npart)=0.0;

func(1:npoints)=0.0; %array of discretized function
fsum(1:npoints)=0.0; %array of cumulative values of the function 

xi=linspace(0,pi/2.0,npoints);

%initial values:
func(1)=a1*(cos(xi(1)))^n1+a2*(cos(xi(1)))^n2;
fsum(1)=0.0; %or func(1)?

for i=2:npoints
    func(i)=a1*(cos(xi(i)))^n1+a2*(cos(xi(i)))^n2;
    fsum(i)=fsum(i-1)+func(i); 
end
%normalize fsum, just in case:

for i=1:npoints
    fsum(i)=fsum(i)/fsum(npoints);
end


for j=1:npart
    u(j)=rand;
    uN(j)=u(j)*(npoints-1);
    k(j)=floor(uN(j))+1;
end

y=k;

end