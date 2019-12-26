% time unit s, length unit mm
clear;
sigma1=2;  % axial stress (MPa)
sigma2=2;
sigma3=2;
D=0.0013;   % mm2/s
S=20*10^(-6);  % mm
OMEGA=2.7*10^4; % mm^3/mol
a(1)=0.1;  % mm
b(1)=0.09; % mm
c(1)=0.08; % mm
ina(1)=0.05;
inb(1)=0.04;
inc(1)=0.03; 
porosity(1)=0;
% a=0.1;  % mm
% b=0.05; % mm
% c=0.025; % mm
% ina=0.01;
% inb=0.005;
% inc=0.0025; 

Cpore=6.48*10^(-6);  % mol/mm^3
R=8.37*10^3;   % mJ/mol/K
T=293;    % K
for i=1:600000;
funVc = @(r) ((a(i)+b(i))*pi.*r-2*pi.*r.^2).^2./(-2*pi.*r+pi*b(i)+2*(a(i)-b(i)));
Vc=sigma3*pi*(a(i)*b(i)-ina(i)*inb(i))*OMEGA^2*S*D/integral(funVc,0,(a(i)-ina(i))/2)/R/T*Cpore;
Vvc=Vc*(a(i)*b(i)*pi-ina(i)*inb(i)*pi)/4;

funVa = @(r) ((b(i)+c(i))*pi.*r-2*pi.*r.^2).^2./(-2*pi.*r+pi*c(i)+2*(b(i)-c(i)));
Va=sigma1*pi*(b(i)*c(i)-inb(i)*inc(i))*OMEGA^2*S*D/integral(funVa,0,(a(i)-ina(i))/2)/R/T*Cpore;
Vva=Va*(c(i)*b(i)*pi-inc(i)*inb(i)*pi)/4;

funVb = @(r) ((a(i)+c(i))*pi.*r-2*pi.*r.^2).^2./(-2*pi.*r+pi*c(i)+2*(a(i)-c(i)));
Vb=sigma2*pi*(a(i)*c(i)-ina(i)*inc(i))*OMEGA^2*S*D/integral(funVb,0,(a(i)-ina(i))/2)/R/T*Cpore;
Vvb=Vb*(a(i)*c(i)*pi-ina(i)*inc(i)*pi)/4;

Vt=Vva+Vvb+Vvc;
Isurface=4*pi*(((ina(i)*inb(i))^(1.6075)+(inc(i)*inb(i))^(1.6075)+(ina(i)*inc(i))^(1.6075))/3)^(1/1.6075);
dr(i)=Vt/Isurface;

porosity(i)=1-(ina(i)*inb(i)*inc(i))/(a(i)*b(i)*c(i));

a(i+1)=a(i)-Va;
b(i+1)=b(i)-Vb;
c(i+1)=c(i)-Vc;
ina(i+1)=ina(i)-Va-2*dr(i);
inb(i+1)=inb(i)-Vb-2*dr(i);
inc(i+1)=inc(i)-Vc-2*dr(i);

if inc(i+1)<0 || ina(i+1)<0 || inb(i+1)<0
   break
end

end

    figure('Name','axis length','NumberTitle','off')
    l1=plot(a,'r')
    hold on 
    l2=plot(b,'g')
    hold on 
    l3=plot(c,'b')
    hold on 
    l4=plot(ina,'k')
    hold on 
    l5=plot(inb,'m')
    hold on 
    l6=plot(inc,'c')
    xlabel('Time(s)','fontsize',16)
    ylabel('Axis length(mm)','fontsize',16)
     legend([l1,l2,l3,l4,l5,l6],'a','b','c','ina','inb','inc','Location','northeast')
    
     figure('Name','porosity','NumberTitle','off')
    l1=plot(porosity,'r')
    xlabel('Time(s)','fontsize',16)
    ylabel('Porosity','fontsize',16)
% r=(a(i)-ina(i))/2;
% integralc=-1/24/pi^3*(6*(2*b(i)+a(i)*(-2+pi))^2*(2*a(i)+b(i)*(-2*pi))*pi*r ...
% +6*(2*b(i)+a(i)*(-2+pi))^2*pi^2*r^2-8*pi^3*(2*a(i)*(-1+pi)+b(i)*(2+pi))*r^3 ...
% +12*r^4*pi^4+3*(2*a(i)^2*(-2+pi)+2*b(i)^2*(-2+pi)+a(i)*b(i)*(8-4*pi+pi^2))^2*log(-2*a(i)+2*b(i)-b(i)*pi+2*pi*r) ...
% -3*(2*a(i)^2*(-2+pi)+2*b(i)^2*(-2+pi)+a(i)*b(i)*(8-4*pi+pi^2))^2*log(-2*a(i)+2*b(i)-b(i)*pi));
% Vc=4*sigma1*pi*(a(i)*b(i)-ina(i)*inb(i))*OMEGA^2*S*D/integralc/R/T*Cpore;


