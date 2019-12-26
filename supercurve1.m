function []=supercurvecase1(creepinterval,cellnumber,j,cellr,thick,epsilon,sigma,totalp)
% cellnumber=50;
ctime=j-1;
ftime=1:1:ctime;
ttime=creepinterval.*ftime;
voidr=zeros(ctime,cellnumber);
celle=zeros(ctime,cellnumber);
celles=zeros(ctime,cellnumber);
cells=zeros(ctime,cellnumber);
cellss=zeros(ctime,cellnumber);
globalp=zeros(1,ctime);
globalve=zeros(1,ctime);
globalh1e=zeros(1,ctime);
globalh2e=zeros(1,ctime);
dglobalve=zeros(1,ctime);
dglobalh1e=zeros(1,ctime);
dglobalh2e=zeros(1,ctime);
globale=cell(1,ctime);
for i=1:ctime
globale(ctime)={zeros(3,3)};
end

for t=1:ctime
    for cn=1:cellnumber
voidr(t,cn)=cellr(cn)-thick(t,cn);
celle(t,cn)=max(eig(epsilon{t,cn}));
celles(t,cn)=min(eig(epsilon{t,cn}));
cells(t,cn)=max(eig(sigma{t,cn}));
cellss(t,cn)=min(eig(sigma{t,cn}));
if cn==1
    globale(t)={epsilon{t,cn}./cellnumber};
else
globale(t)={globale{t}+epsilon{t,cn}./cellnumber};
end
    end
    globalp(t)=totalp(t);
    globalve(t)=globale{t}(3,3);
    globalh1e(t)=globale{t}(1,1);
    globalh2e(t)=globale{t}(2,2);
%     if t==2
%     globalve(1)=globalve(2);
%     globalh1e(1)=globalh1e(2);
%     globalh2e(1)=globalh2e(2);
%     end
    dglobalve(t)=globalve(t)-globalve(1);
    dglobalh1e(t)=globalh1e(t)-globalh1e(1);
    dglobalh2e(t)=globalh2e(t)-globalh2e(1);
end

% curves for each cells' void radius
figure('Name','Void radius','NumberTitle','off')    
voidr1=plot(ttime,voidr(:,cellnumber/10),'r-');
hold on 
voidr2=plot(ttime,voidr(:,cellnumber/10*2),'g-');
hold on 
voidr3=plot(ttime,voidr(:,cellnumber/10*3),'k-');
hold on 
voidr4=plot(ttime,voidr(:,cellnumber/10*4),'b-');
hold on 
voidr5=plot(ttime,voidr(:,cellnumber/10*5),'m-');
hold on 
voidr6=plot(ttime,voidr(:,cellnumber/10*6),'r--');
hold on 
voidr7=plot(ttime,voidr(:,cellnumber/10*7),'g--');
hold on 
voidr8=plot(ttime,voidr(:,cellnumber/10*8),'k--');
hold on 
voidr9=plot(ttime,voidr(:,cellnumber/10*9),'b--');
hold on 
voidr10=plot(ttime,voidr(:,cellnumber),'m--');
hold on 
xlabel('Time(s)','fontsize',16)
ylabel('Void radious(mm)','fontsize',16)
legend([voidr1,voidr2,voidr3,voidr4,voidr5,voidr6,voidr7,voidr8,voidr9,voidr10],'cell5','cell10','cell15','cell20','cell25','cell30','cell35','cell40','cell45','cell50','Location','southeast')
saveas(gcf,'../case1voidr.png');



% curves for each cells' principle strain 
figure('Name','Principal strain','NumberTitle','off')    
% figure('Name','Vertical strain','NumberTitle','off') 
celle1=plot(ttime,celle(:,cellnumber/10),'r-');
hold on 
celle2=plot(ttime,celle(:,cellnumber/10*2),'g-');
hold on 
celle3=plot(ttime,celle(:,cellnumber/10*3),'k-');
hold on 
celle4=plot(ttime,celle(:,cellnumber/10*4),'b-');
hold on 
celle5=plot(ttime,celle(:,cellnumber/10*5),'m-');
hold on 
celle6=plot(ttime,celle(:,cellnumber/10*6),'r--');
hold on 
celle7=plot(ttime,celle(:,cellnumber/10*7),'g--');
hold on 
celle8=plot(ttime,celle(:,cellnumber/10*8),'k--');
hold on 
celle9=plot(ttime,celle(:,cellnumber/10*9),'b--');
hold on 
celle10=plot(ttime,celle(:,cellnumber),'m--');
hold on 
xlabel('Time(s)','fontsize',16)
 ylabel('Principal strain','fontsize',16)
% ylabel('Vertical strain','fontsize',16)
legend([celle1,celle2,celle3,celle4,celle5,celle6,celle7,celle8,celle9,celle10],'cell5','cell10','cell15','cell20','cell25','cell30','cell35','cell40','cell45','cell50','Location','southeast')
saveas(gcf,'../case1celle.png');




% curves for each cells' principle stress 
 figure('Name','Principal stress','NumberTitle','off')    
% figure('Name','Vertical stress','NumberTitle','off')  
cells1=plot(ttime,cells(:,cellnumber/10),'r-');
hold on 
cells2=plot(ttime,cells(:,cellnumber/10*2),'g-');
hold on 
cells3=plot(ttime,cells(:,cellnumber/10*3),'k-');
hold on 
cells4=plot(ttime,cells(:,cellnumber/10*4),'b-');
hold on 
cells5=plot(ttime,cells(:,cellnumber/10*5),'m-');
hold on 
cells6=plot(ttime,cells(:,cellnumber/10*6),'r--');
hold on 
cells7=plot(ttime,cells(:,cellnumber/10*7),'g--');
hold on 
cells8=plot(ttime,cells(:,cellnumber/10*8),'k--');
hold on 
cells9=plot(ttime,cells(:,cellnumber/10*9),'b--');
hold on 
cells10=plot(ttime,cells(:,cellnumber),'m--');
hold on 
xlabel('Time(s)','fontsize',16)
 ylabel('Principal stress','fontsize',16)
% ylabel('Vertical stress(MPa)','fontsize',16)
legend([cells1,cells2,cells3,cells4,cells5,cells6,cells7,cells8,cells9,cells10],'cell5','cell10','cell15','cell20','cell25','cell30','cell35','cell40','cell45','cell50','Location','northeast')
saveas(gcf,'../case1cells.png');



% curves for Global porosity
figure('Name','Global porosity','NumberTitle','off')    
globalp1=plot(ttime,globalp,'r-');
xlabel('Time(s)','fontsize',16)
ylabel('Global porosity','fontsize',16)
saveas(gcf,'../case1globalp.png');



% curves for Global principal strain
% figure('Name','Global principal strain','NumberTitle','off') 
figure('Name','Global vertical strain','NumberTitle','off') 
globalpe1=plot(ttime,globalve,'r-');
xlabel('Time(s)','fontsize',16)
% ylabel('Global principal strain','fontsize',16)
ylabel('Global vertical strain','fontsize',16)
saveas(gcf,'../case1globalve.png');




figure('Name','Global horizontal1 strain','NumberTitle','off') 
globalh1e1=plot(ttime,globalh1e,'r-');
xlabel('Time(s)','fontsize',16)
% ylabel('Global principal strain','fontsize',16)
ylabel('Global horizontal1 strain','fontsize',16)
saveas(gcf,'../case1globalh1e.png');



figure('Name','Global horizontal2  strain','NumberTitle','off') 
globalh2e1=plot(ttime,globalh2e,'r-');
xlabel('Time(s)','fontsize',16)
% ylabel('Global principal strain','fontsize',16)
ylabel('Global horizontal2 strain','fontsize',16)
saveas(gcf,'../case1globalh2e.png');


% cheng's figure
% figure('Name','Principal strain','NumberTitle','off') 
% rl=celle(ctime,:);
% xl=sin(theta).*rl;
% yl=cos(theta).*rl;
% rs=celles(ctime,:);
% xs=sin(theta).*rs;
% ys=cos(theta).*rs;
% re=plot(xl,yl,'r',xs,ys,'b');
% re(1).Marker='*';
% re(2).Marker='*';
% re(1).LineStyle='none';
% re(2).LineStyle='none';
% xlabel('Major principal strain','fontsize',16)
% ylabel('Major principal strain','fontsize',16)
% legend([re(1),re(2)],'Major principal strain','Minor principal strain','Location','southeast')
% 
% figure('Name','Principal stress','NumberTitle','off') 
% rl=cells(ctime,:);
% xl=sin(theta).*rl;
% yl=cos(theta).*rl;
% rs=cellss(ctime,:);
% xs=sin(theta).*rs;
% ys=cos(theta).*rs;
% rrs=plot(xl,yl,'r',xs,ys,'b');
% rrs(1).Marker='*';
% rrs(2).Marker='*';
% rrs(1).LineStyle='none';
% rrs(2).LineStyle='none';
% xlabel('Major principal stress','fontsize',16)
% ylabel('Major principal stress','fontsize',16)
% legend([rrs(1),rrs(2)],'Major principal stress','Minor principal stress','Location','southeast')

%% change
figure('Name','Global vertical creep strain','NumberTitle','off') 
dglobalpe1=plot(ttime,dglobalve,'r-');
xlabel('Time(s)','fontsize',16)
% ylabel('Global principal strain','fontsize',16)
ylabel('Global vertical creep strain','fontsize',16)
saveas(gcf,'../case1globalcve.png');




figure('Name','Global horizontal1 creep strain','NumberTitle','off') 
dglobalh1e1=plot(ttime,dglobalh1e,'r-');
xlabel('Time(s)','fontsize',16)
% ylabel('Global principal strain','fontsize',16)
ylabel('Global horizontal1 creep strain','fontsize',16)
saveas(gcf,'../case1globalch1e.png');



figure('Name','Global horizontal2 creep strain','NumberTitle','off') 
dglobalh2e1=plot(ttime,dglobalh2e,'r-');
xlabel('Time(s)','fontsize',16)
% ylabel('Global principal strain','fontsize',16)
ylabel('Global horizontal2 creep strain','fontsize',16)
saveas(gcf,'../case1globalch2e.png');



