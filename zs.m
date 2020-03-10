function t = zs(x,y,A,SX1,RX1,BX1,CX1,SX2,RX2,BX2,CX2,SY1,RY1,BY1,CY1,SY2,RY2,BY2,CY2,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;  %Non-trenched, smooth surface

%the main trench, (-bx,bx)&(-by,by)
elseif ((x>0)&&(y>0)&&(x<bx)&&(y<by))
    t = (A/2)*(tanh(SX1*(x-bx))+tanh(SY1*(y-by)));
elseif ((x<0)&&(y>0)&&(x>-bx)&&(y<by))
    t = (A/2)*(tanh(-SX2*(x+bx))+tanh(SY1*(y-by))); 
elseif ((x>0)&&(y<0)&&(x<bx)&&(y>-by))
    t = (A/2)*(tanh(SX1*(x-bx))+tanh(-SY2*(y+by)));
elseif ((x<0)&&(y<0)&&(x>-bx)&&(y>-by))
    t = (A/2)*(tanh(-SX2*(x+bx))+tanh(-SY2*(y+by)));    
    
%now add 2nd slope
%not sure which of these geometries works best (derivatives are not affected)
% (a) when both up/downstreams go together
elseif ((x>bx)||(y>by))
    t = (A/2)*(tanh(SX1*(x-bx))+tanh(SY1*(y-by)))+(BX1/2)*tanh(RX1*(x-bx)-CX1)+(BY1/2)*tanh(RY1*(y-by)-CY1);
elseif ((x<-bx)||(y<-by))
    t = (A/2)*(tanh(-SX2*(x+bx))+tanh(-SY2*(y+by)))+(BX2/2)*tanh(-RX2*(x+bx)-CX2)+(BY2/2)*tanh(-RY2*(y+by)-CY2);

%(b) when each side is changed independently
%{
elseif ((x>bx)&&(y>by))
    t = (A/2)*(tanh(SX1*(x-bx))+tanh(SY1*(y-by)))+(BX1/2)*tanh(RX1*(x-bx)-CX1)+(BY1/2)*tanh(RY1*(y-by)-CY1);
elseif ((x<-bx)&&(y<-by))
    t = (A/2)*(tanh(-SX2*(x+bx))+tanh(-SY2*(y+by)))+(BX2/2)*tanh(-RX2*(x+bx)-CX2)+(BY2/2)*tanh(-RY2*(y+by)-CY2);
elseif (x>bx)
    t = (A/2)*(tanh(SX1*(x-bx))+tanh(SY1*(y-by)))+(BX1/2)*tanh(RX1*(x-bx)-CX1);
elseif (y>by)
    t = (A/2)*(tanh(SX1*(x-bx))+tanh(SY1*(y-by)))+(BY1/2)*tanh(RY1*(y-by)-CY1);
elseif (x<-bx)
    t = (A/2)*(tanh(-SX2*(x+bx))+tanh(-SY2*(y+by)))+(BY2/2)*tanh(-RY2*(y+by)-CY2);
elseif (y<-by)
    t = (A/2)*(tanh(-SX2*(x+bx))+tanh(-SY2*(y+by)))+(BY2/2)*tanh(-RY2*(y+by)-CY2);
  %}
else
    t=0;
end 

%add 2nd slope:
%surface: y = A*tanh(S*x) --> y = A*tanh(S*x) + B*tanh(R*x - C)     
%params: SX1 --> SX1,RX1,BX1,CX1 
%        SX2 --> SX2,RX2,BX2,CX2
%        SY1 --> SY1,RY1,BY1,CY1
%        SY2 --> SY2,RY2,BY2,CY2