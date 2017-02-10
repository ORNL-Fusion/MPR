function t = zpx(x,x0,z0,th,phi)
t = z0+(x-x0)/(tan(th)*cos(phi));
end