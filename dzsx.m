function t = dzsx(x,y,A,SX1,SX2,bx,by)
if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif (x>0)
    t = A/2*SX1*(sech(SX1*(x-bx))*sech(SX1*(x-bx)));
else
    t = -A/2*SX2*(sech(-SX2*(x+bx))*sech(-SX2*(x+bx)));
end 