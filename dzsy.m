function t = dzsy(x,y,A,SY1,RY1,BY1,CY1,SY2,RY2,BY2,CY2,bx,by)

if (y>(-x+bx+by)||y>(x+bx+by)||y<(-x-bx-by)||y<(x-bx-by))
    t = 0;
elseif ((y>0)&&(y<by))
    t = (A/2)*SY1*sech(SY1*(y-by))*sech(SY1*(y-by));
elseif (y>by)
    t = (A/2)*SY1*sech(SY1*(y-by))*sech(SY1*(y-by)) + (BY1/2)*RY1*sech(CY1+RY1*(by-y))*sech(CY1+RY1*(by-y));
elseif ((y<0)&&(y>-by))
    t = -(A/2)*SY2*sech(SY2*(by+y))*sech(SY2*(by+y));
elseif (y<-by)
    t = -(A/2)*SY2*sech(SY2*(by+y))*sech(SY2*(by+y)) - (BY2/2)*RY2*sech(CY2+RY2*(by+y))*sech(CY2+RY2*(by+y));
end 

%detailes in derivatives in 'dzsx.m'