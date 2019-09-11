function t = zs(x,y,A,SX1,SX2,SY1,SY2,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;  %Non-trenched, smooth surface
elseif ((x>0)&&(y>0))
    t = A/2*(tanh(SX1*(x-bx))+tanh(SY1*(y-by)));
elseif ((x<0)&&(y>0))
    t = A/2*(tanh(-SX2*(x+bx))+tanh(SY1*(y-by)));
elseif ((x>0)&&(y<0))
    t = A/2*(tanh(SX1*(x-bx))+tanh(-SY2*(y+by)));
else % ((x<0)&&(y<0))
    t = A/2*(tanh(-SX2*(x+bx))+tanh(-SY2*(y+by)));

end 
