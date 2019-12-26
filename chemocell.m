function [nrr,gvce,dfi]=chemocell(gsigma,aa,Q,cthick,t,k)
% function [nrr,vce,dfi,ina,tVa,tVb,tVc]=chemocell(sigma,aa,thick,vs,t)      % test chemo
% time unit s, length unit mm
lsigma=Q*gsigma*Q';
ratio=1;
% modified
porepressure=0.101; % water pressure is 1 atm, 0.101MPa
sigma1=lsigma(1,1)*aa^2/(aa^2-(aa-cthick)^2)*ratio-porepressure;  % axial effective stress of cracks  (MPa) 
sigma2=lsigma(2,2)*aa^2/(aa^2-(aa-cthick)^2)*ratio-porepressure;
sigma3=lsigma(3,3)*aa^2/(aa^2-(aa-cthick)^2)*ratio-porepressure;
% modified


% make DS to be e-19 
D=0.0005*k;   % mm2/s
S=2*10^(-6);  % mm
OMEGA=2.7*10^4; % mm^3/mol
a=zeros(t+1,1);
ina=zeros(t+1,1);
tVc=zeros(t+1,1);
tVa=zeros(t+1,1);
tVb=zeros(t+1,1);
dr=zeros(t,1);
a(1)=aa;  % mm aa is diameter, a is radius

ina(1)=a(1)-cthick;

vctest=0;
% porosity(1)=0;
% a=0.1;  % mm
% b=0.05; % mm
% c=0.025; % mm
% ina=0.01;
% inb=0.005;
% inc=0.0025; 

Cpore=6.48*10^(-6);  % mol/mm^3
R=8.37*10^3;   % mJ/mol/K
T=293;    % K
for i=1:t

Vi=-2*(a(i)^2-ina(i)^2)*OMEGA^2*S*D/R/T/4*Cpore/(a(i)^4*(log(ina(i))-log(a(i)))-a(i)^2*(ina(i)^2-a(i)^2)+(ina(i)^4-a(i)^4)/4);
VVi=(a(i)^2*pi-ina(i)^2*pi);

if sigma3>0
Vc=sigma3*Vi; % dissolution speed perpendicular to crack
Vvc=VVi*Vc;  % dissolution volume along crack
tVc(i+1)=Vc+tVc(i);
else
   Vc=0;
   Vvc=VVi*Vc;  % dissolution volume along crack
   tVc(i+1)=Vc+tVc(i);
end


if sigma1>0
Va=sigma1*Vi;
Vva=VVi*Va;
tVa(i+1)=Va+tVa(i);
else
   Va=0;
Vva=VVi*Va;
tVa(i+1)=Va+tVa(i);
end

if sigma2>0
Vb=sigma2*Vi;
Vvb=VVi*Vb;
tVb(i+1)=Vb+tVb(i);
else
Vb=0;
Vvb=VVi*Vb;
tVb(i+1)=Vb+tVb(i);
end
    
    
Vt=Vva+Vvb+Vvc;
Isurface=4*pi*ina(i)^2;
dr(i)=Vt/Isurface;

% porosity(i)=1-(ina(i)*inb(i)*inc(i))/(a(i)*b(i)*c(i));



if ina(i)-dr(i)<0 
   vctest=1;
   ina(i)=0.00001;
    break
end

a(i+1)=a(i);
ina(i+1)=ina(i)-dr(i);

end

if vctest==0

nrr=a(i+1)-ina(i+1);
lvce=[tVa(i+1)/aa 0 0;0 tVb(i+1)/aa 0;0 0 tVc(i+1)/aa];


end

if vctest==1

nrr=a(i)-ina(i);
lvce=[tVa(i)/aa 0 0;0 tVb(i)/aa 0;0 0 tVc(i)/aa];

end
gvce=Q'*lvce*Q;
dfi=((aa-nrr)^3-ina(1)^3)/aa^3;


end
%     figure('Name','axis length','NumberTitle','off')
%     l1=plot(a,'r')
%     hold on 
%     l2=plot(b,'g')
%     hold on 
%     l3=plot(c,'b')
%     hold on 
%     l4=plot(ina,'k')
%     hold on 
%     l5=plot(inb,'m')
%     hold on 
%     l6=plot(inc,'c')
%     xlabel('Time(s)','fontsize',16)
%     ylabel('Axis length(mm)','fontsize',16)
%      legend([l1,l2,l3,l4,l5,l6],'a','b','c','ina','inb','inc','Location','northeast')
%     
%      figure('Name','porosity','NumberTitle','off')
%     l1=plot(porosity,'r')
%     xlabel('Time(s)','fontsize',16)
%     ylabel('Porosity','fontsize',16)
% r=(a(i)-ina(i))/2;
% integralc=-1/24/pi^3*(6*(2*b(i)+a(i)*(-2+pi))^2*(2*a(i)+b(i)*(-2*pi))*pi*r ...
% +6*(2*b(i)+a(i)*(-2+pi))^2*pi^2*r^2-8*pi^3*(2*a(i)*(-1+pi)+b(i)*(2+pi))*r^3 ...
% +12*r^4*pi^4+3*(2*a(i)^2*(-2+pi)+2*b(i)^2*(-2+pi)+a(i)*b(i)*(8-4*pi+pi^2))^2*log(-2*a(i)+2*b(i)-b(i)*pi+2*pi*r) ...
% -3*(2*a(i)^2*(-2+pi)+2*b(i)^2*(-2+pi)+a(i)*b(i)*(8-4*pi+pi^2))^2*log(-2*a(i)+2*b(i)-b(i)*pi));
% Vc=4*sigma1*pi*(a(i)*b(i)-ina(i)*inb(i))*OMEGA^2*S*D/integralc/R/T*Cpore;


