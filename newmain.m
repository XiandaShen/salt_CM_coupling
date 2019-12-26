%% parameter definition
% stress unit is MPa, length unit is mm
clear;
tic
cellnumber=10; % number of cells in REV
timesteps=35000; % creep time steps
creepinterval=1; % interval for each creep time step
totaliter=1000;  % iteration numer
muo=8846; % initial mu

psi=zeros(1,cellnumber);
phi=zeros(1,cellnumber);
theta=zeros(1,cellnumber);
% orientaion of each cell
% for i=1:cellnumber
% % phi(i)=mod(i,5);
% theta(i)=fix((i-1)/5);
% 
% 
% end



Q=cell(1,cellnumber); % transformation matrix
for i=1:cellnumber
    Q(i)={transmatrix(psi(i),theta(i),phi(i))};
  %    Q(i)={eye(3)};  % no orientation
%    Q{i}*Q{i}'
end    



dinter=cell(totaliter,cellnumber);  %  delta stress for each cell, column means cell orders, row means iteraion times
dsigma=cell(timesteps,cellnumber); % delta stress for each cell, column means cell orders
beforeisigma=cell(1,cellnumber); % start sigma for each timestep iteration
sigma=cell(timesteps,cellnumber); % stress for each cell
    for i=1:cellnumber
    sigma(1,i)={zeros(3,3)};
    end
depsilon=cell(timesteps,cellnumber); % delta strain for each cell
epsilon=cell(timesteps,cellnumber); % strain for each cell
    for i=1:cellnumber
    epsilon(1,i)={zeros(3,3)};
    end
dvcepsilon=cell(timesteps,cellnumber); % delta viscostrain for each cell
depsilone=cell(timesteps,cellnumber); % delta elastic strain for each cell
vcepsilon=cell(timesteps,cellnumber); % viscostrain for each cell
    for i=1:cellnumber
    vcepsilon(1,i)={zeros(3,3)};
    end
epsilone=cell(timesteps,cellnumber); % elastic strain for each cell
    for i=1:cellnumber
    epsilone(1,i)={zeros(3,3)};
    end
% dtsigma=[0.1 0 0; 0 0.1 0; 0 0 0.1]; % total delta sigma



F0=zeros(6*cellnumber,1); % new equations vector error
x0=zeros(6*cellnumber,1); % new dsigma iteration vector error
F1=zeros(6*cellnumber,1); % current equations vector error 
x1=zeros(6*cellnumber,1); % current dsigma iteration vector error
xnew=zeros(6*cellnumber,1); % current dsigma iteration vector error
invBFGS=eye(6*cellnumber);
% dtsigma=[0.05 0 0; 0 0.02 0; 0 0 0.01]; % total delta sigma
dtsigma=[0 0 0; 0 0 0; 0 0 0]; % total delta sigma
dtsigmai=[0 0 0; 0 0 0; 0 0 0.1]; % mechanical delta sigma
tsigmai=[0 0 0; 0 0 0; 0 0 9.9]; % total mechanical delta sigma
% dtepsilon=zeros(3,3); % total delta strain
tepsilon=zeros(3,3); % total strain 
tvolume=zeros(1,timesteps); % total volume



% modified
% intial size of cell


cellr=0.1375*ones(1,cellnumber); 

gooddata=0;
while gooddata==0
gooddata=1;
realm = 0.103125;
realv = 0.0001;
logm = log((realm^2)/sqrt(realv+realm^2));
logsigma = sqrt(log(realv/(realm^2)+1));
iniX = lognrnd(logm,logsigma,1,cellnumber);

for i=1:cellnumber
    if iniX(i)>0.13
       gooddata=0 ;
    end
end
if abs(mean(iniX)-realm)>0.1 || abs(var(iniX)-realv)>0.00002
          gooddata=0 ;
     
    
end
end

thick=cellr-iniX;


% modified end 

% cellr=0.3*ones(1,cellnumber);
% thick=0.20:0.01:(0.01*(cellnumber-1)+0.20);  
% thick=0.25:0.005:(0.005*(cellnumber-1)+0.25);%  thick=0.3-0.05/cellnumber:-0.05/cellnumber:0.25;
% thick=0.3-0.1/cellnumber:-0.1/cellnumber:0.2;
% thick=0.2*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size
% for i=1:cellnumber
% vsolid(i)=4/3*(cellr(i)^3-(cellr(i)-thick(i))^3);
% end
% cella=0.3:0.01:(0.01*(cellnumber-1)+0.3); 
% cellb=0.28:0.01:(0.01*(cellnumber-1)+0.28);
% cellc=0.27:-0.005:(-0.005*(cellnumber-1)+0.27);
% thick=0.01*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

%  cella=0.1:0.01:(0.01*(cellnumber-1)+0.1); 
%  cellb=0.0999:0.01:(0.01*(cellnumber-1)+0.0999);
%  cellc=0.0998:0.01:(0.01*(cellnumber-1)+0.0998);
%  thick=0.005*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

%  cella=0.3*ones(1,cellnumber); 
%  cellb=0.2999*ones(1,cellnumber);
%  cellc=0.2998*ones(1,cellnumber);
%  thick=0.01*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

% cellS=cell(1,cellnumber); % Eshelby tensor for each cell
cellL=cell(1,cellnumber); % Stiffness tensor for each cell
cellLF=cell(1,cellnumber); % Stiffness tensor for each cell Last step
cellH=cell(1,cellnumber); % Hill tensor for each cell
cellP=zeros(1,cellnumber); % cell porosity
totalp=zeros(1,timesteps+1);  % total porosity

dHdL=doubledotef(majorsymidendityeighth,(inversegeneral(sphereS)-symidendityf));
%% First increment 
totalp(1)=totalporosity(cellr(1,:),thick(1,:),cellnumber);
totalL=stiffness(muo,totalp(1));
dtepsilon=100*doubledotft(inversegeneral(totalL),dtsigmai);
tepsilon=tepsilon+dtepsilon; 

for cn=1:cellnumber
    cellP(cn)=porosity(cellr(cn),thick(1,cn));  
    cellL(cn)={stiffness(muo,cellP(cn))};
    cellH(cn)={sphereH(muo,totalp(1))};
    sigmaright=symidendityf+doubledotff(cellH{cn},inversegeneral(totalL));
    sigmaleft=inversegeneral(symidendityf+doubledotff(cellH{cn},inversegeneral(cellL{cn})));
    dsigma(1,cn)={doubledotft(doubledotff(sigmaleft,sigmaright),dtsigmai)};
    
    % strian and stress accumulation, notice it is 100 times 
    depsilon(1,cn)={doubledotft(inversegeneral(cellL{cn}),dsigma{1,cn})};
    depsilone(1,cn)=depsilon(1,cn);
    epsilon(1,cn)={100*depsilon{1,cn}};
    epsilone(1,cn)=epsilon(1,cn);
    sigma(1,cn)={100*dsigma{1,cn}};
    dsigma(1,cn)={zeros(3,3)};
    
end


%% Creep
for j=1:timesteps-1
beforeisigma=sigma(j,:);
storetep=tepsilon; % avoid iteration changes tepsilon
dinter(1,:)=dsigma(j,:);
storex1=zeros(6*cellnumber,totaliter);
storeF1=zeros(6*cellnumber,totaliter);
arfak=1;
iteralpha=0;
% update property for each cell

checkdts=cell(totaliter,1); % check the change of total stress
for i=1:totaliter
    checkdts(i)={zeros(3,3)};
end

% iteration 
for itnum=1:totaliter

    tepsilon=storetep; 
 



if itnum==1
for cn=1:cellnumber    
    xcn=dinter{itnum,cn};
    x1(6*(cn-1)+1)=xcn(1,1);
    x1(6*(cn-1)+2)=xcn(1,2);
    x1(6*(cn-1)+3)=xcn(1,3);
    x1(6*(cn-1)+4)=xcn(2,2);
    x1(6*(cn-1)+5)=xcn(2,3);
    x1(6*(cn-1)+6)=xcn(3,3);
%    load('/Users/caihejushi/Documents/academic/PHD/chemo_self_consistent/general_ellipsoid/cmx1.mat')
end   

% pure mechanics
% [F1,totalp,tepsilon,dinter,dsigma,depsilon,epsilon,cella,cellb,cellc,thick]=puremechanics(j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter,x1);
% chemomechanics
 [F1,totalp,tepsilon,dinter,dsigma,sigma,depsilon,epsilon,epsilone,cellr,thick]=chemomechanics(creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,x1,muo);


end

pk=-invBFGS*F1;
% BFGS=inv(invBFGS);
% alpha iteration
% if itnum>50
% [arfak, iteralpha] = back_tracting(iteralpha, pk, F1, x1, j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter);
% end
%  arfak
% arfa=1;
sk=arfak*pk;
xnew=sk+x1;

% pure mechanics
%[Fnew,totalp,tepsilon,dinter,dsigma,depsilon,epsilon,cella,cellb,cellc,thick]=puremechanics(j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter,xnew);
 % chemomechanics
 [Fnew,totalp,tepsilon,dinter,dsigma,sigma,depsilon,epsilon,epsilone,cellr,thick]=chemomechanics(creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,xnew,muo);



yk=Fnew-F1;
invBFGS=invBFGS+(sk'*yk+yk'*invBFGS*yk)*(sk*sk')/(sk'*yk)^2-(invBFGS*yk*sk'+sk*yk'*invBFGS)/(sk'*yk);
storex1(:,itnum)=x1(:,1);
storeF1(:,itnum)=F1(:,1);
F1=Fnew;
x1=xnew;
if norm(F1,'inf')<1e-10
     itnum
j
%     iteralpha
    break
end
end


storetep=tepsilon; % avoid iteration changes tepsilon
dsigma(j+1,:)=dinter(itnum,:);
for cn=1:cellnumber 
sigma(j+1,cn)={beforeisigma{1,cn}+dinter{itnum,cn}}; 
end

if totalp(j)<1e-9
    break;
end
end
toc
save('C:\Users\xshen60\Downloads\oriented sphere\50totalplo.mat','totalp');
save('C:\Users\xshen60\Downloads\oriented sphere\50sigmalo.mat','sigma');
save('C:\Users\xshen60\Downloads\oriented sphere\50epsilonlo.mat','epsilon');
save('C:\Users\xshen60\Downloads\oriented sphere\50thicklo.mat','thick');