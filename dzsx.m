function t = dzsx(x,y,A,SX1,RX1,BX1,CX1,SX2,RX2,BX2,CX2,bx,by)

if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif ((x>0)&&(x<bx))
    t = (A/2)*SX1*sech(SX1*(x-bx))*sech(SX1*(x-bx));
elseif (x>bx)
    t = (A/2)*SX1*sech(SX1*(x-bx))*sech(SX1*(x-bx)) + (BX1/2)*RX1*sech(CX1+RX1*(bx-x))*sech(CX1+RX1*(bx-x));
elseif ((x<0)&&(x>-bx))
    t = -(A/2)*SX2*sech(SX2*(bx+x))*sech(SX2*(bx+x));
elseif (x<-bx)
    t = -(A/2)*SX2*sech(SX2*(bx+x))*sech(SX2*(bx+x)) - (BX2/2)*RX2*sech(CX2+RX2*(bx+x))*sech(CX2+RX2*(bx+x));
end 

% y = A*tanh(S*x) --> y = A*tanh(S*x) + B*tanh(R*x - C)  = f(x)+g(x)
% d(tanh(x))/dx=(sech(x))^2

% f(x) = A*tanh(x+C)
% df/dt= A*sech^2(C+x)

% g(x) = (B*2/A)*tanh(R*(x-bx)-C)
%dg/dx= (2*B*R/A)*sech(C + R(bx-x)*sech(C + R*(bx-x))

%a)
%t(x)= A/2*(tanh(S*(x-bx))+(B*2/A)*tanh(R*(x-bx)-C)
%dt/dx=A/2*(S*sech^2(S*(x-bx))+(B*2*R/A)*sech^2(C+R*(bx-x)))

%b)
%t=A/2*(tanh(-S*(x+bx))-(B*2/A)*tanh(-R*(x+bx)+C)
%dt/dx= -A/2*(S*sech^2(S*(bx + x)) - (B*2/A)*R*sech^2(C-R*(bx + x)))