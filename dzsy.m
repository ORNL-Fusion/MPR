function t = dzsy(x,y,A,SY1,SY2,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif (y>0)
    t = A/2*SY1*(sech(SY1*(y-by))*sech(SY1*(y-by)));
else
    t = -A/2*SY2*(sech(-SY2*(y+by))*sech(-SY2*(y+by)));
end 
