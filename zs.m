function t = zs(x,y,A,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif ((x>0)&&(y>0))
    t = A/2*(tanh(x-bx)+tanh(y-by));
elseif ((x<0)&&(y>0))
    t = A/2*(tanh(-x-bx)+tanh(y-by));
elseif ((x>0)&&(y<0))
    t = A/2*(tanh(x-bx)+tanh(-y-by));
else % ((x<0)&&(y<0))
    t = A/2*(tanh(-x-bx)+tanh(-y-by));

end 
