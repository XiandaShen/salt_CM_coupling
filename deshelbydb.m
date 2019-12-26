% eshelby tensor
function dSdb=deshelbydb(a,b,c)
v=0.3;
theta=asin(sqrt(1-c^2/a^2));
k=sqrt(a^2-b^2)/sqrt(a^2-c^2);
F=ellipticF(theta,k^2);
E=ellipticE(theta,k^2);

dkdb=-((a^2-b^2)*(a^2-c^2))^(-0.5)*b;
fundFdb = @(w) (1-k.^2.*(sin(w)).^2).^(-1.5).*(sin(w)).^2;
dFdb=k*dkdb*integral(fundFdb,0,theta);
fundEdb = @(w) (1-k.^2.*(sin(w)).^2).^(-0.5).*(sin(w)).^2;
dEdb=-k*dkdb*integral(fundEdb,0,theta);

Ia=4*pi*a*b*c/sqrt(a^2-c^2)/(a^2-b^2)*(F-E);
Ic=4*pi*a*b*c/sqrt(a^2-c^2)/(b^2-c^2)*(b*sqrt(a^2-c^2)/a/c-E);
Ib=4*pi-Ia-Ic;

dIadb=(dFdb-dEdb)*4*pi*a*b*c/(a^2-b^2)/sqrt(a^2-c^2)+(a^2-c^2)^(-0.5)*(4*pi*a*c/(a^2-b^2)+8*pi*a*b^2*c/(a^2-b^2)^2)*(F-E);
dIcdb=((a^2-c^2)^0.5/a/c-dEdb)*4*pi*a*b*c/(b^2-c^2)/sqrt(a^2-c^2)+(a^2-c^2)^(-0.5)*(4*pi*a*c/(b^2-c^2)-8*pi*a*b^2*c/(b^2-c^2)^2)*(b*(a^2-c^2)^0.5/a/c-E);
dIbdb=-dIadb-dIcdb;

Iab=(Ib-Ia)/(3*a^2-3*b^2);
Iac=(Ic-Ia)/(3*a^2-3*c^2);
Iaa=4*pi/3/a^2-Iab-Iac;

dIabdb=2*b/3/(a^2-b^2)^2*(Ib-Ia)+(dIbdb-dIadb)/3/(a^2-b^2);
dIacdb=(dIcdb-dIadb)/(3*a^2-3*c^2);
dIaadb=-dIabdb-dIacdb;

Ibc=(Ic-Ib)/(3*b^2-3*c^2);
Ibb=4*pi/3/b^2-Iab-Ibc;
Icc=4*pi/3/c^2-Iac-Ibc;

dIbcdb=-2*b/3/(b^2-c^2)^2*(Ic-Ib)+(dIcdb-dIbdb)/3/(b^2-c^2);
dIbbdb=-8*pi/3/b^3-dIabdb-dIbcdb;
dIccdb=-dIacdb-dIbcdb;

dSdb(1:3,1:3,1:3,1:3)=0;
Q=3/(8*pi*(1-v));
R=(1-2*v)/(8*pi*(1-v));

dSdb(1,1,1,1)=Q*a^2*dIaadb+R*dIadb;
dSdb(1,1,1,2)=0;
dSdb(1,1,1,3)=0;
dSdb(1,1,2,1)=0;
dSdb(1,1,2,2)=Q*b^2*dIabdb+2*Q*b*Iab-R*dIadb;
dSdb(1,1,2,3)=0;
dSdb(1,1,3,1)=0;
dSdb(1,1,3,2)=0;
dSdb(1,1,3,3)=Q*c^2*dIacdb-R*dIadb;
dSdb(1,2,1,1)=0;
dSdb(1,2,1,2)=Q*(a^2+b^2)/2*dIabdb+Q*b*Iab+R*(Ia+dIbdb)/2;
dSdb(1,2,1,3)=0;
dSdb(1,2,2,1)=Q*(a^2+b^2)/2*dIabdb+Q*b*Iab+R*(Ia+dIbdb)/2;
dSdb(1,2,2,2)=0;
dSdb(1,2,2,3)=0;
dSdb(1,2,3,1)=0;
dSdb(1,2,3,2)=0;
dSdb(1,2,3,3)=0;
dSdb(1,3,1,1)=0;
dSdb(1,3,1,2)=0;
dSdb(1,3,1,3)=Q/2*(a^2+c^2)*dIacdb+R/2*(dIadb+dIcdb);
dSdb(1,3,2,1)=0;
dSdb(1,3,2,2)=0;
dSdb(1,3,2,3)=0;
dSdb(1,3,3,1)=Q/2*(a^2+c^2)*dIacdb+R/2*(dIadb+dIcdb);
dSdb(1,3,3,2)=0;
dSdb(1,3,3,3)=0;
dSdb(2,1,1,1)=0;
dSdb(2,1,1,2)=Q*(a^2+b^2)/2*dIabdb+Q*b*Iab+R*(Ia+dIbdb)/2;
dSdb(2,1,1,3)=0;
dSdb(2,1,3,1)=0;
dSdb(2,1,3,2)=0;
dSdb(2,1,3,3)=0;
dSdb(2,1,2,1)=Q*(a^2+b^2)/2*dIabdb+Q*b*Iab+R*(Ia+dIbdb)/2;
dSdb(2,1,2,2)=0;
dSdb(2,1,2,3)=0;
dSdb(2,1,3,1)=0;
dSdb(2,1,3,2)=0;
dSdb(2,1,3,3)=0;
dSdb(2,2,1,1)=Q*a^2*dIabdb-R*dIbdb;
dSdb(2,2,1,2)=0;
dSdb(2,2,1,3)=0;
dSdb(2,2,2,1)=0;
dSdb(2,2,2,2)=Q*b^2*dIbbdb+Q*2*b*Ibb+R*dIbdb;
dSdb(2,2,2,3)=0;
dSdb(2,2,3,1)=0;
dSdb(2,2,3,2)=0;
dSdb(2,2,3,3)=Q*c^2*dIbcdb-R*dIbdb;
dSdb(2,3,1,1)=0;
dSdb(2,3,1,2)=0;
dSdb(2,3,1,3)=0;
dSdb(2,3,2,1)=0;
dSdb(2,3,2,2)=0;
dSdb(2,3,2,3)=Q/2*(b^2+c^2)*dIbcdb+Q*b*Ibc+R/2*(dIbdb+dIcdb);
dSdb(2,3,3,1)=0;
dSdb(2,3,3,2)=Q/2*(b^2+c^2)*dIbcdb+Q*b*Ibc+R/2*(dIbdb+dIcdb);
dSdb(2,3,3,3)=0;
dSdb(3,1,1,1)=0;
dSdb(3,1,1,2)=0;
dSdb(3,1,1,3)=Q/2*(a^2+c^2)*dIacdb+R/2*(dIadb+dIcdb);
dSdb(3,1,2,1)=0;
dSdb(3,1,2,2)=0;
dSdb(3,1,2,3)=0;
dSdb(3,1,3,1)=Q/2*(a^2+c^2)*dIacdb+R/2*(dIadb+dIcdb);
dSdb(3,1,3,2)=0;
dSdb(3,1,3,3)=0;
dSdb(3,2,1,1)=0;
dSdb(3,2,1,2)=0;
dSdb(3,2,1,3)=0;
dSdb(3,2,2,1)=0;
dSdb(3,2,2,2)=0;
dSdb(3,2,2,3)=Q/2*(b^2+c^2)*dIbcdb+Q*b*Ibc+R/2*(dIbdb+dIcdb);
dSdb(3,2,3,1)=0;
dSdb(3,2,3,2)=Q/2*(b^2+c^2)*dIbcdb+Q*b*Ibc+R/2*(dIbdb+dIcdb);
dSdb(3,2,3,3)=0;
dSdb(3,3,1,1)=Q*a^2*dIacdb-R*dIcdb;
dSdb(3,3,1,2)=0;
dSdb(3,3,1,3)=0;
dSdb(3,3,2,1)=0;
dSdb(3,3,2,2)=Q*b^2*dIbcdb+Q*b*2*Ibc-R*dIcdb;
dSdb(3,3,3,1)=0;
dSdb(3,3,3,2)=0;
dSdb(3,3,3,3)=Q*c^2*dIccdb+R*dIcdb;
end




