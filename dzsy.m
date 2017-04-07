function t = dzsy(x,y,A,S,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif (y>0)
    t = A/2*S*(sech(S*(y-by))*sech(S*(y-by)));
else
    t = -A/2*S*(sech(-S*(y+by))*sech(-S*(y+by)));
end 
