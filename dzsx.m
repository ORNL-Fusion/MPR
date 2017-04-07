function t = dzsx(x,y,A,S,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif (x>0)
    t = A/2*S*(sech(S*(x-bx))*sech(S*(x-bx)));
else
    t = -A/2*S*(sech(-S*(x+bx))*sech(-S*(x+bx)));
end 