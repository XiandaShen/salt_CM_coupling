function [FF1,totalp,tepsilon,dinter,dsigma,depsilon,epsilon,cella,cellb,cellc,thick]=puremechanics(j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter,xx1)

aaa=zeros(1,cellnumber); % define a b c da db dc dfi for each time step, avoid reading large matrix over times
bbb=zeros(1,cellnumber);
ccc=zeros(1,cellnumber);
rrr=zeros(1,cellnumber);
dtemptepsilon=zeros(3,3); % temper dtepsilon for summation
cellS=cell(1,cellnumber); % Eshelby tensor for each cell
cellL=cell(1,cellnumber); % Stiffness tensor for each cell
cellLF=cell(1,cellnumber); % Stiffness tensor for each cell Last step
cellH=cell(1,cellnumber); % Hill tensor for each cell
cellP=zeros(1,cellnumber); % cell porosity
FF1=zeros(3*cellnumber,1); % current equations vector error 


for putback=1:cellnumber   
    dinter(itnum,putback)={[xx1(3*(putback-1)+1) 0 0;0 xx1(3*(putback-1)+2) 0;0 0 xx1(3*(putback-1)+3)]};
   dsigma(j,putback)=dinter(itnum,putback); 
 end  

for cn=1:cellnumber
cella(j+1,cn)=cella(j,cn);
cellb(j+1,cn)=cellb(j,cn);
cellc(j+1,cn)=cellc(j,cn);
thick(j+1,cn)=thick(j,cn);
aaa(cn)=cella(j+1,cn);
bbb(cn)=cellb(j+1,cn);
ccc(cn)=cellc(j+1,cn);
rrr(cn)=thick(j+1,cn);

% update porosity, stiffness, eshelby tensor
cellP(cn)=porosity(aaa(cn),bbb(cn),ccc(cn),rrr(cn));  
cellL(cn)={stiffness(cellP(cn))};
cellLF(cn)=cellL(cn);
cellS(cn)={eshelby(aaa(cn),bbb(cn),ccc(cn))}; 
% update strain
depsilon(j+1,cn)={doubledotft(inversegeneral(cellL{cn}),dsigma{j,cn})};
% sum strain and volume 
dtemptepsilon=dtemptepsilon+depsilon{j+1,cn};
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

    % construct a vector with 0 value equations
 %   dinter(itnum,cn)={dtsigma-doubledotft(cellH{cn},(depsilon{j+1,cn}-dtepsilon))-doubledotft(dcellH,(epsilon{j+1,cn}-tepsilon))}; 
     Fcn=dinter{itnum,cn}-dtsigma-doubledotft(cellH{cn},(depsilon{j+1,cn}-dtepsilon));
     FF1(3*(cn-1)+1)=Fcn(1,1);
     FF1(3*(cn-1)+2)=Fcn(2,2);
     FF1(3*(cn-1)+3)=Fcn(3,3);

     
end
end