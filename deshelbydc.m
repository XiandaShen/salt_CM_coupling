% eshelby tensor
function dSdc=deshelbydc(a,b,c)
v=0.3;
theta=asin(sqrt(1-c^2/a^2));
k=sqrt(a^2-b^2)/sqrt(a^2-c^2);
F=ellipticF(theta,k^2);
E=ellipticE(theta,k^2);

funFc=@(w) (1-k.^2.*(sin(w)).^2).^(-0.5);
Fc=integral(funFc,0,theta);
funEc=@(w) (1-k.^2.*(sin(w)).^2).^(0.5);
Ec=integral(funEc,0,theta);

dkdc=c*(a^2-b^2)^0.5/(a^2-c^2)^1.5;
dthetadc=-1/a/(1-c^2/a^2)^0.5;
fundFdc = @(w) (1-k.^2.*(sin(w)).^2).^(-1.5).*(sin(w)).^2;
dFdc=k*dkdc*integral(fundFdc,0,theta)+dthetadc/(1-k^2*(sin(theta))^2)^0.5;
fundEdc = @(w) (1-k.^2.*(sin(w)).^2).^(-0.5).*(sin(w)).^2;
dEdc=-k*dkdc*integral(fundEdc,0,theta)+dthetadc*(1-k^2*(sin(theta))^2)^0.5;

Ia=4*pi*a*b*c/sqrt(a^2-c^2)/(a^2-b^2)*(F-E);
Ic=4*pi*a*b*c/sqrt(a^2-c^2)/(b^2-c^2)*(b*sqrt(a^2-c^2)/a/c-E);
Ib=4*pi-Ia-Ic;

dIadc=(dFdc-dEdc)*4*pi*a*b*c/(a^2-b^2)/sqrt(a^2-c^2)+(a^2-b^2)^(-1)*(4*pi*a*b/sqrt(a^2-c^2)+4*pi*a*b*c^2/(a^2-c^2)^1.5)*(F-E);
dIcdc=8*pi*b^2*c/(b^2-c^2)^2-dEdc*4*pi*a*b*c/(b^2-c^2)/sqrt(a^2-c^2)-E*(4*pi*a*b/(b^2-c^2)/sqrt(a^2-c^2)+8*pi*a*b*c^2/(b^2-c^2)^2/sqrt(a^2-c^2)+4*pi*a*b*c^2/(b^2-c^2)/(a^2-c^2)^1.5);
dIbdc=-dIadc-dIcdc;

Iab=(Ib-Ia)/(3*a^2-3*b^2);
Iac=(Ic-Ia)/(3*a^2-3*c^2);
Iaa=4*pi/3/a^2-Iab-Iac;

dIabdc=(dIbdc-dIadc)/(3*a^2-3*b^2);
dIacdc=(dIcdc-dIadc)/(3*a^2-3*c^2)+(Ic-Ia)/3/(a^2-c^2)^2*2*c;
dIaadc=-dIabdc-dIacdc;

Ibc=(Ic-Ib)/(3*b^2-3*c^2);
Ibb=4*pi/3/b^2-Iab-Ibc;
Icc=4*pi/3/c^2-Iac-Ibc;

dIbcdc=(dIcdc-dIbdc)/(3*b^2-3*c^2)+(Ic-Ib)/3/(b^2-c^2)^2*2*c;
dIbbdc=-dIabdc-dIbcdc;
dIccdc=-8*pi/3/c^3-dIacdc-dIbcdc;

dSdc(1:3,1:3,1:3,1:3)=0;
Q=3/(8*pi*(1-v));
R=(1-2*v)/(8*pi*(1-v));

dSdc(1,1,1,1)=Q*a^2*dIaadc+R*dIadc;
dSdc(1,1,1,2)=0;
dSdc(1,1,1,3)=0;
dSdc(1,1,2,1)=0;
dSdc(1,1,2,2)=Q*b^2*dIabdc-R*dIadc;
dSdc(1,1,2,3)=0;
dSdc(1,1,3,1)=0;
dSdc(1,1,3,2)=0;
dSdc(1,1,3,3)=Q*c^2*dIacdc+Q*2*c*Iac-R*dIadc;
dSdc(1,2,1,1)=0;
dSdc(1,2,1,2)=Q*(a^2+b^2)/2*dIabdc+R*(dIadc+dIbdc)/2;
dSdc(1,2,1,3)=0;
dSdc(1,2,2,1)=Q*(a^2+b^2)/2*dIabdc+R*(dIadc+dIbdc)/2;
dSdc(1,2,2,2)=0;
dSdc(1,2,2,3)=0;
dSdc(1,2,3,1)=0;
dSdc(1,2,3,2)=0;
dSdc(1,2,3,3)=0;
dSdc(1,3,1,1)=0;
dSdc(1,3,1,2)=0;
dSdc(1,3,1,3)=Q/2*(a^2+c^2)*dIacdc+Q*c*Iac+R/2*(dIadc+dIcdc);
dSdc(1,3,2,1)=0;
dSdc(1,3,2,2)=0;
dSdc(1,3,2,3)=0;
dSdc(1,3,3,1)=Q/2*(a^2+c^2)*dIacdc+Q*c*Iac+R/2*(dIadc+dIcdc);
dSdc(1,3,3,2)=0;
dSdc(1,3,3,3)=0;
dSdc(2,1,1,1)=0;
dSdc(2,1,1,2)=Q*(a^2+b^2)/2*dIabdc+R*(dIadc+dIbdc)/2;
dSdc(2,1,1,3)=0;
dSdc(2,1,3,1)=0;
dSdc(2,1,3,2)=0;
dSdc(2,1,3,3)=0;
dSdc(2,1,2,1)=Q*(a^2+b^2)/2*dIabdc+R*(dIadc+dIbdc)/2;
dSdc(2,1,2,2)=0;
dSdc(2,1,2,3)=0;
dSdc(2,1,3,1)=0;
dSdc(2,1,3,2)=0;
dSdc(2,1,3,3)=0;
dSdc(2,2,1,1)=Q*a^2*dIabdc-R*dIbdc;
dSdc(2,2,1,2)=0;
dSdc(2,2,1,3)=0;
dSdc(2,2,2,1)=0;
dSdc(2,2,2,2)=Q*b^2*dIbbdc+R*dIbdc;
dSdc(2,2,2,3)=0;
dSdc(2,2,3,1)=0;
dSdc(2,2,3,2)=0;
dSdc(2,2,3,3)=Q*c^2*dIbcdc+Q*2*c*Ibc-R*dIbdc;
dSdc(2,3,1,1)=0;
dSdc(2,3,1,2)=0;
dSdc(2,3,1,3)=0;
dSdc(2,3,2,1)=0;
dSdc(2,3,2,2)=0;
dSdc(2,3,2,3)=Q/2*(b^2+c^2)*dIbcdc+Q*c*Ibc+R/2*(dIbdc+dIcdc);
dSdc(2,3,3,1)=0;
dSdc(2,3,3,2)=Q/2*(b^2+c^2)*dIbcdc+Q*c*Ibc+R/2*(dIbdc+dIcdc);
dSdc(2,3,3,3)=0;
dSdc(3,1,1,1)=0;
dSdc(3,1,1,2)=0;
dSdc(3,1,1,3)=Q/2*(a^2+c^2)*dIacdc+Q*c*Iac+R/2*(dIadc+dIcdc);
dSdc(3,1,2,1)=0;
dSdc(3,1,2,2)=0;
dSdc(3,1,2,3)=0;
dSdc(3,1,3,1)=Q/2*(a^2+c^2)*dIacdc+Q*c*Iac+R/2*(dIadc+dIcdc);
dSdc(3,1,3,2)=0;
dSdc(3,1,3,3)=0;
dSdc(3,2,1,1)=0;
dSdc(3,2,1,2)=0;
dSdc(3,2,1,3)=0;
dSdc(3,2,2,1)=0;
dSdc(3,2,2,2)=0;
dSdc(3,2,2,3)=Q/2*(b^2+c^2)*dIbcdc+Q*c*Ibc+R/2*(dIbdc+dIcdc);
dSdc(3,2,3,1)=0;
dSdc(3,2,3,2)=Q/2*(b^2+c^2)*dIbcdc+Q*c*Ibc+R/2*(dIbdc+dIcdc);
dSdc(3,2,3,3)=0;
dSdc(3,3,1,1)=Q*a^2*dIacdc-R*dIcdc;
dSdc(3,3,1,2)=0;
dSdc(3,3,1,3)=0;
dSdc(3,3,2,1)=0;
dSdc(3,3,2,2)=Q*b^2*dIbcdc-R*dIcdc;
dSdc(3,3,3,1)=0;
dSdc(3,3,3,2)=0;
dSdc(3,3,3,3)=Q*c^2*dIccdc+Q*2*c*Icc+R*dIcdc;
end



