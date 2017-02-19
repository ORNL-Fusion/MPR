function t = dzsy(x,y,A,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif (y>0)
    t = A/2*(sech(y-by)*sech(y-by));
else
    t = -A/2*(sech(-y-by)*sech(-y-by));
end 
