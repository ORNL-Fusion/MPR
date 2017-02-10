%generate random numbers with an analytically given distribution:
function y = cosn_distrib(a1, n1, a2, n2, npart, npoints)

u(1,npart)=0.0; %array of npart random number, uniformly distributed (0,1)

output(1:npart)=0.0;

func(1:npoints)=0.0; %array of discretized function
fsum(1:npoints)=0.0; %array of cumulative values of the function 

cosxi=linspace(0,pi/2.0,npoints);


%initial values:
func(1)=a1*(cos(cosxi(1)))^n1+a2*(cos(cosxi(1)))^n2;
fsum(1)=0.0; %or func(1)?

for i=2:npoints
    func(i)=a1*(cos(cosxi(i)))^n1+a2*(cos(cosxi(i)))^n2;
    fsum(i)=fsum(i-1)+func(i); 
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
        Stest=[ 'angle binary search for particle ' , num2str(j), ' u(j)=', num2str(u(j))];
        Stest2=['kmin=', num2str(kmin), ' kmid=' , num2str(kmin), ' kmax=' , num2str(kmax), ' fsum(kmin)=', num2str(fsum(kmin)), ' fsum(kmid)=' , num2str(fsum(kmin)), ' fsum(kmax)=' , num2str(fsum(kmax)) ];
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
        Sxi=['therefore the angle is: x(k)=', num2str(cosxi(kmid))];
        disp(Sbin)           
        disp(Sxi)
    end
    %}
    
    output(j)=cosxi(kmid);
    
end

y=output;

end