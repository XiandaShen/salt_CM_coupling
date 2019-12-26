% eshelby tensor
function dSda=deshelbyda(a,b,c)
v=0.3;
theta=asin(sqrt(1-c^2/a^2));
k=sqrt(a^2-b^2)/sqrt(a^2-c^2);
F=ellipticF(theta,k^2);
E=ellipticE(theta,k^2);

dkda=-a*(a^2-b^2)^(-0.5)*(a^2-c^2)^(-1.5)*(c^2-b^2);
dthetada=1/a^2*c*(1-c^2/a^2)^(-0.5);
fundFda = @(w) (1-k.^2.*(sin(w)).^2).^(-1.5).*(sin(w)).^2;
dFda=k*dkda*integral(fundFda,0,theta)+dthetada/(1-k^2*(sin(theta))^2)^0.5;
fundEda = @(w) (1-k.^2.*(sin(w)).^2).^(-0.5).*(sin(w)).^2;
dEda=-k*dkda*integral(fundEda,0,theta)+dthetada*(1-k^2*(sin(theta))^2)^0.5;

Ia=4*pi*a*b*c/sqrt(a^2-c^2)/(a^2-b^2)*(F-E);
Ic=4*pi*a*b*c/sqrt(a^2-c^2)/(b^2-c^2)*(b*sqrt(a^2-c^2)/a/c-E);
Ib=4*pi-Ia-Ic;

dIada=(dFda-dEda)*4*pi*a*b*c/(a^2-b^2)/sqrt(a^2-c^2)+(4*pi*c*b/(a^2-b^2)/sqrt(a^2-c^2)-4*pi*c*b*a^2/(a^2-b^2)/(a^2-c^2)^1.5-8*pi*c*b*a^2/(a^2-b^2)^2/(a^2-c^2)^0.5)*(F-E);
dIcda=-dEda*4*pi*a*b*c/(b^2-c^2)/sqrt(a^2-c^2)-E/(b^2-c^2)*(4*pi*b*c/(a^2-c^2)^0.5-4*pi*a^2*b*c/(a^2-c^2)^1.5);
dIbda=-dIada-dIcda;

Iab=(Ib-Ia)/(3*a^2-3*b^2);
Iac=(Ic-Ia)/(3*a^2-3*c^2);
Iaa=4*pi/3/a^2-Iab-Iac;

dIabda=(dIbda-dIada)/3/(a^2-b^2)-2*(Ib-Ia)*a/3/(a^2-b^2)^2;
dIacda=(dIcda-dIada)/3/(a^2-c^2)-2*(Ic-Ia)*a/3/(a^2-c^2)^2;
dIaada=-8*pi/3/a^3-dIabda-dIacda;

Ibc=(Ic-Ib)/(3*b^2-3*c^2);
Ibb=4*pi/3/b^2-Iab-Ibc;
Icc=4*pi/3/c^2-Iac-Ibc;

dIbcda=(dIcda-dIbda)/(3*b^2-3*c^2);
dIbbda=-dIabda-dIbcda;
dIccda=-dIacda-dIbcda;


dSda(1:3,1:3,1:3,1:3)=0;
Q=3/(8*pi*(1-v));
R=(1-2*v)/(8*pi*(1-v));

dSda(1,1,1,1)=Q*a^2*dIaada+Q*2*a*Iaa+R*Ia;
dSda(1,1,1,2)=0;
dSda(1,1,1,3)=0;
dSda(1,1,2,1)=0;
dSda(1,1,2,2)=Q*b^2*dIabda-R*dIada;
dSda(1,1,2,3)=0;
dSda(1,1,3,1)=0;
dSda(1,1,3,2)=0;
dSda(1,1,3,3)=Q*c^2*dIacda-R*dIada;
dSda(1,2,1,1)=0;
dSda(1,2,1,2)=Q*a*Iab+Q*(a^2+b^2)/2*dIabda+R*(dIada+dIbda)/2;
dSda(1,2,1,3)=0;
dSda(1,2,2,1)=Q*a*Iab+Q*(a^2+b^2)/2*dIabda+R*(dIada+dIbda)/2;
dSda(1,2,2,2)=0;
dSda(1,2,2,3)=0;
dSda(1,2,3,1)=0;
dSda(1,2,3,2)=0;
dSda(1,2,3,3)=0;
dSda(1,3,1,1)=0;
dSda(1,3,1,2)=0;
dSda(1,3,1,3)=Q*a*Iac+Q/2*(a^2+c^2)*dIacda+R/2*(dIada+dIcda);
dSda(1,3,2,1)=0;
dSda(1,3,2,2)=0;
dSda(1,3,2,3)=0;
dSda(1,3,3,1)=Q*a*Iac+Q/2*(a^2+c^2)*dIacda+R/2*(dIada+dIcda);
dSda(1,3,3,2)=0;
dSda(1,3,3,3)=0;
dSda(2,1,1,1)=0;
dSda(2,1,1,2)=Q*a*Iab+Q*(a^2+b^2)/2*dIabda+R*(dIada+dIbda)/2;
dSda(2,1,1,3)=0;
dSda(2,1,3,1)=0;
dSda(2,1,3,2)=0;
dSda(2,1,3,3)=0;
dSda(2,1,2,1)=Q*a*Iab+Q*(a^2+b^2)/2*dIabda+R*(dIada+dIbda)/2;
dSda(2,1,2,2)=0;
dSda(2,1,2,3)=0;
dSda(2,1,3,1)=0;
dSda(2,1,3,2)=0;
dSda(2,1,3,3)=0;
dSda(2,2,1,1)=Q*2*a*Iab+Q*a^2*dIabda-R*dIbda;
dSda(2,2,1,2)=0;
dSda(2,2,1,3)=0;
dSda(2,2,2,1)=0;
dSda(2,2,2,2)=Q*b^2*dIbbda+R*dIbda;
dSda(2,2,2,3)=0;
dSda(2,2,3,1)=0;
dSda(2,2,3,2)=0;
dSda(2,2,3,3)=Q*c^2*dIbcda-R*dIbda;
dSda(2,3,1,1)=0;
dSda(2,3,1,2)=0;
dSda(2,3,1,3)=0;
dSda(2,3,2,1)=0;
dSda(2,3,2,2)=0;
dSda(2,3,2,3)=Q/2*(b^2+c^2)*dIbcda+R/2*(dIbda+dIcda);
dSda(2,3,3,1)=0;
dSda(2,3,3,2)=Q/2*(b^2+c^2)*dIbcda+R/2*(dIbda+dIcda);
dSda(2,3,3,3)=0;
dSda(3,1,1,1)=0;
dSda(3,1,1,2)=0;
dSda(3,1,1,3)=Q*a*Iac+Q/2*(a^2+c^2)*dIacda+R/2*(dIada+dIcda);
dSda(3,1,2,1)=0;
dSda(3,1,2,2)=0;
dSda(3,1,2,3)=0;
dSda(3,1,3,1)=Q*a*Iac+Q/2*(a^2+c^2)*dIacda+R/2*(dIada+dIcda);
dSda(3,1,3,2)=0;
dSda(3,1,3,3)=0;
dSda(3,2,1,1)=0;
dSda(3,2,1,2)=0;
dSda(3,2,1,3)=0;
dSda(3,2,2,1)=0;
dSda(3,2,2,2)=0;
dSda(3,2,2,3)=Q/2*(b^2+c^2)*dIbcda+R/2*(dIbda+dIcda);
dSda(3,2,3,1)=0;
dSda(3,2,3,2)=Q/2*(b^2+c^2)*dIbcda+R/2*(dIbda+dIcda);
dSda(3,2,3,3)=0;
dSda(3,3,1,1)=Q*a^2*dIacda+Q*2*a*Iac-R*dIcda;
dSda(3,3,1,2)=0;
dSda(3,3,1,3)=0;
dSda(3,3,2,1)=0;
dSda(3,3,2,2)=Q*b^2*dIbcda-R*dIcda;
dSda(3,3,3,1)=0;
dSda(3,3,3,2)=0;
dSda(3,3,3,3)=Q*c^2*dIccda+R*dIcda;
end




