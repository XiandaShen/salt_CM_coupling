%% parameter definition
% stress unit is MPa, length unit is mm
clear;
eir=0.5;  % explicit implicit ratio
cellnumber=10; % number of cells in REV
timesteps=2; % creep time steps
creepinterval=10; % interval for each creep time step

dsigma=cell(timesteps,cellnumber); % delta stress for each cell, column means cell orders
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
dtsigma=[0 0 0; 0 0 0; 0 0 0]; % total delta sigma
dtsigmai=[0.1 0 0; 0 0.1 0; 0 0 0.1]; % mechanical delta sigma
tsigmai=[9.9 0 0; 0 9.9 0; 0 0 9.9]; % total mechanical delta sigma
% dtepsilon=zeros(3,3); % total delta strain
tepsilon=zeros(3,3); % total strain 
tvolume=zeros(1,timesteps); % total volume
% intial size of cell
cella=0.3:0.01:(0.01*(cellnumber-1)+0.3); 
cellb=0.28:0.01:(0.01*(cellnumber-1)+0.28);
cellc=0.27:-0.005:(-0.005*(cellnumber-1)+0.27);
thick=0.01*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

%  cella=0.1:0.01:(0.01*(cellnumber-1)+0.1); 
%  cellb=0.0999:0.01:(0.01*(cellnumber-1)+0.0999);
%  cellc=0.0998:0.01:(0.01*(cellnumber-1)+0.0998);
%  thick=0.005*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

%  cella=0.3*ones(1,cellnumber); 
%  cellb=0.2999*ones(1,cellnumber);
%  cellc=0.2998*ones(1,cellnumber);
%  thick=0.01*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

cellS=cell(1,cellnumber); % Eshelby tensor for each cell
cellL=cell(1,cellnumber); % Stiffness tensor for each cell
cellLF=cell(1,cellnumber); % Stiffness tensor for each cell Last step
cellH=cell(1,cellnumber); % Hill tensor for each cell
cellP=zeros(1,cellnumber); % cell porosity
totalp=zeros(1,timesteps+1);  % total porosity

%% First increment 
totalp(1)=totalporosity(cella(1,:),cellb(1,:),cellc(1,:),thick(1,:),cellnumber);
totalL=stiffness(totalp(1));
dtepsilon=doubledotft(inversegeneral(totalL),dtsigmai);
tepsilon=tepsilon+dtepsilon; 

for cn=1:cellnumber
    cellS(cn)={eshelby(cella(1,cn),cellb(1,cn),cellc(1,cn))}; % for cell data, {} is the content of cell, () is the location of cell
    cellP(cn)=porosity(cella(1,cn),cellb(1,cn),cellc(1,cn),thick(1,cn));  
    cellL(cn)={stiffness(cellP(cn))};
    cellH(cn)={doubledotff(totalL,inversegeneral(cellS{cn}))-totalL};
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
storetep=tepsilon; % avoid iteration changes tepsilon
aaa=zeros(1,cellnumber); % define a b c da db dc dfi for each time step, avoid reading large matrix over times
bbb=zeros(1,cellnumber);
ccc=zeros(1,cellnumber);
rrr=zeros(1,cellnumber);
da=zeros(1,cellnumber); 
db=zeros(1,cellnumber);
dc=zeros(1,cellnumber);
dfi=zeros(1,cellnumber);
% update property for each cell
totaliter=50;
checkdts=cell(totaliter,1); % check the change of total stress
for i=1:totaliter
    checkdts(i)={zeros(3,3)};
end

record=zeros(totaliter,1);
dinter=cell(totaliter,cellnumber);  
siter=cell(totaliter,cellnumber);
for itnum=1:totaliter
dtemptepsilon=zeros(3,3); % temper dtepsilon for summation
tepsilon=storetep; 

for cn=1:cellnumber

% update a, b, c, r, vcstrain 
[cella(j+1,cn),cellb(j+1,cn),cellc(j+1,cn),thick(j+1,cn),vce,da(cn),db(cn),dc(cn),dfi(cn)]=chemocell(sigma{j,cn},cella(j,cn),cellb(j,cn),cellc(j,cn),thick(j,cn),creepinterval);
aaa(cn)=cella(j+1,cn);
bbb(cn)=cellb(j+1,cn);
ccc(cn)=cellc(j+1,cn);
rrr(cn)=thick(j+1,cn);
dvcepsilon(j+1,cn)={vce};
% update porosity, stiffness, eshelby tensor
cellP(cn)=porosity(aaa(cn),bbb(cn),ccc(cn),rrr(cn));  
cellL(cn)={stiffness(cellP(cn))};
cellLF(cn)=cellL(cn);
cellS(cn)={eshelby(aaa(cn),bbb(cn),ccc(cn))}; 
% update strain, cell volume
depsilone(j+1,cn)={doubledotft(inversegeneral(2*cellL{cn}-cellLF{cn}),(dsigma{j,cn}-doubledotft((cellL{cn}-cellLF{cn}),epsilone{j,cn})))};
depsilon(j+1,cn)={depsilone{j+1,cn}+dvcepsilon{j+1,cn}};
epsilone(j+1,cn)={epsilone{j,cn}+depsilone{j+1,cn}};
cellvolume=4/3*pi*aaa(cn)*bbb(cn)*ccc(cn);
% sum strain and volume 
dtemptepsilon=dtemptepsilon+depsilon{j+1,cn};
% dtemptepsilon=dtemptepsilon+cellvolume*depsilon{j+1,cn};
tvolume(j+1)=tvolume(j+1)+cellvolume;
epsilon(j+1,cn)={epsilon{j,cn}+depsilon{j+1,cn}};
end

% update global delta strain, global strain
dtepsilon=dtemptepsilon/cellnumber;
% dtepsilon=dtemptepsilon/tvolume(j+1);
tepsilon=tepsilon+dtepsilon;
% update global porosity, global stiffness
totalp(j+1)=totalporosity(aaa,bbb,ccc,rrr,cellnumber);
totalL=stiffness(totalp(j+1));

% update cellS, cellP,cellL,cellH,dcellH,dsigma
for cn=1:cellnumber
    cellH(cn)={doubledotff(totalL,inversegeneral(cellS{cn}))-totalL};
    dSa=deshelbyda(aaa(cn),bbb(cn),ccc(cn));
    dSb=deshelbydb(aaa(cn),bbb(cn),ccc(cn));
    dSc=deshelbydc(aaa(cn),bbb(cn),ccc(cn));
    dinSS=dineshelbydS(inversegeneral(cellS{cn}));
    dHdS=doubledotfe(totalL,dinSS);
    dHdL=doubledotef(majorsymidendityeighth,(inversegeneral(cellS{cn})-symidendityf));
    dcellH=fourdotef(dHdS,deshelbyda(aaa(cn),bbb(cn),ccc(cn)))*da(cn)+fourdotef(dHdS,deshelbydb(aaa(cn),bbb(cn),ccc(cn)))*db(cn)+fourdotef(dHdS,deshelbydc(aaa(cn),bbb(cn),ccc(cn)))*dc(cn) ...
        +fourdotef(dHdL,dstiffnessdfi)*(totalp(j+1)-totalp(j));
     dinter(itnum,cn)={dtsigma-doubledotft(cellH{cn},(depsilon{j+1,cn}-dtepsilon))-doubledotft(dcellH,(epsilon{j+1,cn}-tepsilon))};
     
     % 
     stablecheck=dsigma{j,cn};
     sigma(j,cn)={sigma{j,cn}+eir*dinter{itnum,cn}}; 
     siter(itnum,cn)=sigma(j,cn);
    dsigma(j,cn)={eir*dinter{itnum,cn}};
    checkdts(itnum)={checkdts{itnum}+dsigma{j,cn}};
    
%      dsigma(j+1,cn)={dtsigma-doubledotft(dcellH,(epsilon{j+1,cn}-tepsilon))-doubledotft(cellH{cn},(depsilon{j+1,cn}-dtepsilon))};
%      sigma(j+1,cn)={sigma{j,cn}+dsigma{j+1,cn}};
     
end

if abs((stablecheck(1,1)-dsigma{j,cn}(1,1))/stablecheck(1,1))<0.1
    eir=0.6;
end
if abs((stablecheck(1,1)-dsigma{j,cn}(1,1))/stablecheck(1,1))<0.05
    eir=0.75;
end
if abs((stablecheck(1,1)-dsigma{j,cn}(1,1))/stablecheck(1,1))<0.02
    eir=0.9;
end
if abs((stablecheck(1,1)-dsigma{j,cn}(1,1))/stablecheck(1,1))<0.01
    eir=1;
end
record(itnum)=eir;
% update global vcepsilon

% uddate cell stressp;
end
 save(['/Users/caihejushi/Documents/academic/PHD/chemo_self_consistent/dinter',num2str(j),'.mat'],'dinter');


end















% scell=cell(2,2)
% scell(1,1)={[2 2; 2 2]}
