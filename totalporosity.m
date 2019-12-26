function tp=totalporosity(r,t,n)
tt=0;
tv=0;
for i=1:n
tt=tt+4/3*pi*r(i)^3;
tv=tv+4/3*pi*(r(i)-t(i))^3;
end
tp=tv/tt;
end