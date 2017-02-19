function t = dzsx(x,y,A,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif (x>0)
    t = A/2*(sech(x-bx)*sech(x-bx));
else
    t = -A/2*(sech(-x-bx)*sech(-x-bx));
end 