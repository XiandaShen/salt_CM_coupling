%% en means for consolidation 6n dsigma+ global strain (3, 3)
%% parameter definition
% stress unit is MPa, length unit is mm
clear;
tic
cellnumber=50; % number of cells in REV
timesteps=2500; % creep time steps
creepinterval=100; % interval for each creep time step
totaliter=5000;  % iteration numer
muo=8846; % initial mu
whatever=1;
psi=zeros(1,cellnumber);
phi=zeros(1,cellnumber);
theta=zeros(1,cellnumber);
failsign=0;
% orientaion of each cell
%  for i=1:cellnumber
% % % phi(i)=mod(i,5);
%  theta(i)=pi/2/cellnumber*fix(i-1);
% % 
% % 
%  end
phi=pi.*rand(1,cellnumber);
theta=pi/2.*rand(1,cellnumber);



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



F0=zeros(6*cellnumber+1,1); % new equations vector error
x0=zeros(6*cellnumber+1,1); % new dsigma iteration vector error
F1=zeros(6*cellnumber+1,1); % current equations vector error 
x1=zeros(6*cellnumber+1,1); % current dsigma iteration vector error
xnew=zeros(6*cellnumber+1,1); % current dsigma iteration vector error
invBFGS=eye(6*cellnumber+1);
% dtsigma=[0.05 0 0; 0 0.02 0; 0 0 0.01]; % total delta sigma
dtsigma=[0 0 0; 0 0 0; 0 0 0]; % total delta sigma
dtsigmai=[0 0 0; 0 0 0; 0 0 3.1]; % mechanical delta sigma
% tsigmai=[0 0 0; 0 0 0; 0 0 9.9]; % total mechanical delta sigma
dtepsilon=zeros(3,3); % total delta strain
tepsilon=zeros(3,3); % total strain 
tvolume=zeros(1,timesteps); % total volume



% modified
% intial size of cell
% 
% 
cellr=0.1375*ones(1,cellnumber); 

gooddata=0;
while gooddata==0
gooddata=1;
realm = 0.103125;
% realm = 0.103125;
realv = 0.0005;
logm = log((realm^2)/sqrt(realv+realm^2));
logsigma = sqrt(log(realv/(realm^2)+1));
iniX = lognrnd(logm,logsigma,1,cellnumber);

for i=1:cellnumber
    if iniX(i)>0.13
       gooddata=0 ;
    end
end
if abs(mean(iniX)-realm)>0.1 || abs(var(iniX)-realv)>0.0001
          gooddata=0 ;
     
    
end
end
figure 
plot (iniX)
mean (iniX)
var(iniX)
thick=cellr-iniX;


% modified end 
% 
%  cellr=0.1375*ones(1,cellnumber); 
%   thick=0.034*ones(1,cellnumber);
% thick=0.20:0.01:(0.01*(cellnumber-1)+0.20);  
% thick=0.25:0.005:(0.005*(cellnumber-1)+0.25);
% thick=0.3-0.05/cellnumber:-0.05/cellnumber:0.25;
% thick=0.3-0.1/cellnumber:-0.1/cellnumber:0.2;
% thick=0.2*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

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


% modified
dtepsilon(3,3)=dtsigmai(3,3)/totalL(3,3,3,3);
dtsigmai(1,1)=totalL(1,1,3,3)*dtepsilon(3,3); % get the expression for confining global stress
dtsigmai(2,2)=totalL(2,2,3,3)*dtepsilon(3,3);
% modified

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
    epsilon(1,cn)=depsilon(1,cn);
    epsilone(1,cn)=epsilon(1,cn);
    sigma(1,cn)=dsigma(1,cn);
    dsigma(1,cn)={zeros(3,3)};
    dtsigma=dtsigma+dsigma{1,cn}./cellnumber;
end


%% Creep
for j=1:timesteps-1
    testdis=0;
    disa=1;
    itnum=0;
    redotime=0;
beforeisigma=sigma(j,:);
storetep=tepsilon; % avoid iteration changes tepsilon
dinter(1,:)=dsigma(j,:);
storex1=zeros(6*cellnumber+1,totaliter);
storeF1=zeros(6*cellnumber+1,totaliter);
arfak=1;
iteralpha=0;
% update property for each cell

checkdts=cell(totaliter,1); % check the change of total stress
for i=1:totaliter
    checkdts(i)={zeros(3,3)};
end

% iteration 
% for itnum=1:totaliter
while whatever==1
    itnum=itnum+1;
    
    
    
    
    
    tepsilon=storetep; 
 



if itnum==1
% for cn=1:cellnumber    
%     xcn=dinter{itnum,cn};
%     x1(6*(cn-1)+1)=xcn(1,1);
%     x1(6*(cn-1)+2)=xcn(1,2);
%     x1(6*(cn-1)+3)=xcn(1,3);
%     x1(6*(cn-1)+4)=xcn(2,2);
%     x1(6*(cn-1)+5)=xcn(2,3);
%     x1(6*(cn-1)+6)=xcn(3,3);
% %     x1(6*cellnumber+1)=x1(6*cellnumber+1)+x1(6*(cn-1)+1);
% %     x1(6*cellnumber+2)=x1(6*cellnumber+2)+x1(6*(cn-1)+2);
% %     x1(6*cellnumber+3)=x1(6*cellnumber+3)+x1(6*(cn-1)+3);
% %     x1(6*cellnumber+4)=x1(6*cellnumber+4)+x1(6*(cn-1)+4);
% %     x1(6*cellnumber+5)=x1(6*cellnumber+5)+x1(6*(cn-1)+5);
% %    load('/Users/caihejushi/Documents/academic/PHD/chemo_self_consistent/general_ellipsoid/cmx1.mat')
% end   





% pure mechanics
% [F1,totalp,tepsilon,dinter,dsigma,depsilon,epsilon,cella,cellb,cellc,thick]=puremechanics(j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter,x1);
% chemomechanics
 [F1,totalp,tepsilon,dinter,dsigma,sigma,depsilon,epsilon,epsilone,cellr,thick,dvcepsilon]=encalibrationchemomechanics(creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,x1,muo);


end

pk=-invBFGS*F1;
% BFGS=inv(invBFGS);
% alpha iteration
% if itnum>50
% [arfak, iteralpha] = back_tracting(iteralpha, pk, F1, x1, j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter);
% end
%  arfak
% arfa=1;


% 
% if j>1
% [arfak, iteralpha] = back_tracting( iteralpha, pk, F1, x1, creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,muo); 
% itnum
% arfak
% iteralpha
% end
% 
if norm(x1,'inf')>500

disa=disa+0.05;
    x1=orgx1./disa;

    iteralpha=0;
    invBFGS=eye(6*cellnumber+1);
     [F1]=endiscalibrationchemomechanics(creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,x1,muo);
testdis=1;
pk=-invBFGS*F1;
end

sk=arfak*pk;
xnew=sk+x1;

% pure mechanics
%[Fnew,totalp,tepsilon,dinter,dsigma,depsilon,epsilon,cella,cellb,cellc,thick]=puremechanics(j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter,xnew);
 % chemomechanics
 [Fnew,totalp,tepsilon,dinter,dsigma,sigma,depsilon,epsilon,epsilone,cellr,thick,dvcepsilon]=encalibrationchemomechanics(creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,xnew,muo);



yk=Fnew-F1;
invBFGS=invBFGS+(sk'*yk+yk'*invBFGS*yk)*(sk*sk')/(sk'*yk)^2-(invBFGS*yk*sk'+sk*yk'*invBFGS)/(sk'*yk);

F1=Fnew;
x1=xnew;
storex1(:,itnum)=x1(:,1);
storeF1(:,itnum)=F1(:,1);

% Fstress=zeros(cellnumber,1);
% Fstrain=zeros(5*cellnumber,1);
% 
% for erri=1:cellnumber
% Fstress(erri)=F1((erri-1)*6+1);
% end
% 
% for erri=1:cellnumber*5
% Fstrain(erri)=F1(erri+fix((erri-1)/5)+1);
% end




if norm(F1,'inf')<1e-14
      itnum;
%1  1    
      j;
    break
end

if itnum>totaliter
failsign=1;
break
end


% if itnum==totaliter
%     itnum
%     itnum=1;
%     redotime=redotime+1;
% end
% 
% if redotime==10
%     break
% end



end
testdis;
orgx1=x1;
orginvBFGS=invBFGS;
% storetep=tepsilon; % avoid iteration changes tepsilon
dsigma(j+1,:)=dinter(itnum,:);
for cn=1:cellnumber 
sigma(j+1,cn)={beforeisigma{1,cn}+dinter{itnum,cn}}; 
end

if totalp(j)<1e-9 || failsign==1
    break;
end



end

time=toc;
save('../case1time.mat','time');
save('../case1totalp.mat','totalp');
save('../case1sigma.mat','sigma');
save('../case1epsilon.mat','epsilon');
save('../case1thick.mat','thick');


supercurvecase1(creepinterval,cellnumber,j,cellr,thick,epsilon,sigma,totalp)





