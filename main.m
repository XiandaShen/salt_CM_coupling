%% parameter definition
% stress unit is MPa, length unit is mm
clear;
cellnumber=10; % number of cells in REV
timesteps=22; % creep time steps
creepinterval=10; % interval for each creep time step
dsigma=cell(1,cellnumber); % delta stress for each cell, column means cell orders
sigma=cell(1,cellnumber); % stress for each cell
    for i=1:cellnumber
    sigma(i)={zeros(3,3)};
    end
depsilon=cell(1,cellnumber); % delta strain for each cell
epsilon=cell(1,cellnumber); % strain for each cell
    for i=1:cellnumber
    epsilon(i)={zeros(3,3)};
    end
dvcepsilon=cell(1,cellnumber); % delta viscostrain for each cell
vcepsilon=cell(1,cellnumber); % viscostrain for each cell
    for i=1:cellnumber
    vcepsilon(i)={zeros(3,3)};
    end
% dtsigma=[0.1 0 0; 0 0.1 0; 0 0 0.1]; % total delta sigma
dtsigma=[0 0 0; 0 0 0; 0 0 0]; % total delta sigma
dtsigmai=[10 0 0; 0 10 0; 0 0 10]; % total delta sigma
dtepsilon=zeros(3,3); % total delta strain
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

%  cella=1.1*ones(1,cellnumber); 
%  cellb=1.09999*ones(1,cellnumber);
%  cellc=1.09998*ones(1,cellnumber);
%  thick=0.005*ones(1,cellnumber); % thickness of each cell, can be used to calculate the iner abc, it should not be larger than cell size

cellS=cell(1,cellnumber); % Eshelby tensor for each cell
cellL=cell(1,cellnumber); % Stiffness tensor for each cell
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
    dsigma(cn)={doubledotft(doubledotff(sigmaleft,sigmaright),dtsigmai)};
    
    % strian and stress accumulation
    depsilon(cn)={doubledotft(inversegeneral(cellL{cn}),dsigma{cn})};
    epsilon(cn)={epsilon{cn}+depsilon{cn}};
    sigma(cn)={sigma{cn}+dsigma{cn}};
end


%% Creep
for j=1:timesteps
dtemptepsilon=zeros(3,3); % temper dtepsilon for summation
aaa=zeros(1,cellnumber); % define a b c da db dc dfi for each time step, avoid reading large matrix over times
bbb=zeros(1,cellnumber);
ccc=zeros(1,cellnumber);
rrr=zeros(1,cellnumber);
da=zeros(1,cellnumber); 
db=zeros(1,cellnumber);
dc=zeros(1,cellnumber);
dfi=zeros(1,cellnumber);
% update property for each cell
for cn=1:cellnumber
% update a, b, c, r, vcstrain 
[cella(j+1,cn),cellb(j+1,cn),cellc(j+1,cn),thick(j+1,cn),vce,da(cn),db(cn),dc(cn),dfi(cn)]=chemocell(sigma{cn},cella(j,cn),cellb(j,cn),cellc(j,cn),thick(j,cn),creepinterval);
aaa(cn)=cella(j+1,cn);
bbb(cn)=cellb(j+1,cn);
ccc(cn)=cellc(j+1,cn);
rrr(cn)=thick(j+1,cn);
dvcepsilon(cn)={vce};
% update porosity, stiffness, eshelby tensor
cellP(cn)=porosity(aaa(cn),bbb(cn),ccc(cn),rrr(cn));  
cellL(cn)={stiffness(cellP(cn))};
cellS(cn)={eshelby(aaa(cn),bbb(cn),ccc(cn))}; 
% update strain, cell volume
depsilon(cn)={doubledotft(inversegeneral(cellL{cn}),dsigma{cn})+dvcepsilon{cn}};
epsilon(cn)={epsilon{cn}+depsilon{cn}};
cellvolume=4/3*pi*aaa(cn)*bbb(cn)*ccc(cn);
% sum strain and volume 
dtemptepsilon=dtemptepsilon+cellvolume*depsilon{cn};
tvolume(j)=tvolume(j)+cellvolume;
end

% update global delta strain, global strain
dtepsilon=dtemptepsilon/tvolume(j);
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
        +fourdotef(dHdL,dstiffnessdfi)*dfi(cn);
     dsigma(cn)={dtsigma-doubledotft(dcellH,(epsilon{cn}-tepsilon))-doubledotft(cellH{cn},(depsilon{cn}-dtepsilon))};
     sigma(cn)={sigma{cn}+dsigma{cn}};
     % dsigma(cn)={dtsigma-doubledotft(dcellH,(epsilon{cn}-tepsilon))-doubledotft(cellH{cn},(depsilon{cn}-dtepsilon))};
end

% update global vcepsilon

% uddate cell stressp;

end















% scell=cell(2,2)
% scell(1,1)={[2 2; 2 2]}
