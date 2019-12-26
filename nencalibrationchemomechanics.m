function [FF1,totalp,tepsilon,dinter,dsigma,sigma,depsilon,epsilon,epsilone,cellr,thick,dvcepsilon]=nencalibrationchemomechanics(creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,xx1,muo,DSs)

 % define r dr t dfi for each time step, avoid reading large matrix over times

 dtemptepsilon=zeros(3,3); % temper dtepsilon for summation
% cellS=cell(1,cellnumber); % Eshelby tensor for each cell
cellL=cell(1,cellnumber); % Stiffness tensor for each cell
cellLF=cell(1,cellnumber); % Stiffness tensor for each cell Last step

cellH=cell(1,cellnumber); % Hill tensor for each cell
cellP=zeros(1,cellnumber); % cell porosity
FF1=zeros(6*cellnumber,1); % current equations vector error 
% dr=zeros(1,cellnumber); 
dfi=zeros(1,cellnumber);
dtepsilon=zeros(3,3);

for putback=1:cellnumber   
    dinter(1,putback)={[xx1(6*(putback-1)+1) xx1(6*(putback-1)+2) xx1(6*(putback-1)+3);xx1(6*(putback-1)+2) xx1(6*(putback-1)+4) xx1(6*(putback-1)+5);xx1(6*(putback-1)+3) xx1(6*(putback-1)+5) xx1(6*(putback-1)+6)]};
   dsigma(j,putback)=dinter(1,putback); 
       sigma(j+1,putback)={sigma{j,putback}+dsigma{j,putback}};
       dtsigma=dtsigma+dsigma{j,putback}./cellnumber;
end  

dtsigma(3,3)=0;
% dtepsilon(3,3)=xx1(cellnumber*6+1)/1000;







for cn=1:cellnumber


 [thick(j+1,cn),vce,dfi(cn)]=chemocell(sigma{j+1,cn},cellr(cn),Q{cn},thick(j,cn),creepinterval,DSs);
% thick(j+1,cn)=thick(j,cn);
% dfi(cn)=0;

%!!!!!!!!!!!!    
% vce=zeros(3,3);



% S does not change
% cella(j+1,cn)=cella(j,cn);
% cellb(j+1,cn)=cellb(j,cn);
% cellc(j+1,cn)=cellc(j,cn);
% da(cn)=0;
% db(cn)=0;
% dc(cn)=0;
% 
% aaa(cn)=cella(j+1,cn);
% bbb(cn)=cellb(j+1,cn);
% ccc(cn)=cellc(j+1,cn);
% rrr(cn)=thick(j+1,cn);

dvcepsilon(j+1,cn)={vce};
% update porosity, stiffness, eshelby tensor
cellLF(cn)={stiffness(muo,porosity(cellr(cn),thick(j,cn)))};
cellP(cn)=porosity(cellr(cn),thick(j+1,cn));


%!!!!!!!!!!!!    
 cellL(cn)={stiffness(muo,cellP(cn))};
% cellL(cn)=cellLF(cn);



% update strain, cell volume
% depsilone(j+1,cn)={doubledotft(inversegeneral(2*cellL{cn}-cellLF{cn}),(dsigma{j,cn}-doubledotft((cellL{cn}-cellLF{cn}),epsilone{j,cn})))};


%!!!!!!!!!!!! 
 depsilone(j+1,cn)={doubledotft(inversegeneral(2*cellL{cn}-cellLF{cn}),(dsigma{j,cn}-doubledotft((cellL{cn}-cellLF{cn}),epsilone{j,cn})))};
% depsilone(j+1,cn)={doubledotft(inversegeneral(cellL{cn}),dsigma{j,cn})};


depsilon(j+1,cn)={depsilone{j+1,cn}+dvcepsilon{j+1,cn}};
epsilone(j+1,cn)={epsilone{j,cn}+depsilone{j+1,cn}};

% sum strain and volume 
 dtemptepsilon=dtemptepsilon+depsilon{j+1,cn};
epsilon(j+1,cn)={epsilon{j,cn}+depsilon{j+1,cn}};
end

% update global delta strain, global strain
 dtepsilon(3,3)=dtemptepsilon(3,3)/cellnumber;
tepsilon=tepsilon+dtepsilon;

% update global porosity, global stiffness





 totalp(j+1)=totalporosity(cellr(1,:),thick(j+1,:),cellnumber);
% totalp(j+1)=totalp(j);
% totalL=stiffness(totalp(j+1));

% update cellS, cellP,cellL,cellH,dcellH,dsigma
for cn=1:cellnumber

    
    cellH(cn)={sphereH(muo,totalp(j+1))};

    dcellH=fourdotef(dHdL,dstiffnessdfi)*(totalp(j+1)-totalp(j));
     
    
    Fcn=-dinter{1,cn}+dtsigma-doubledotft(dcellH,(epsilon{j+1,cn}-tepsilon))-doubledotft(cellH{cn},(depsilon{j+1,cn}-dtepsilon));
 
    
    
    
    % construct a vector with 0 value equations
 %   dinter(itnum,cn)={dtsigma-doubledotft(cellH{cn},(depsilon{j+1,cn}-dtepsilon))-doubledotft(dcellH,(epsilon{j+1,cn}-tepsilon))}; 

 % modified
    FF1(6*(cn-1)+1)=Fcn(1,1);
    FF1(6*(cn-1)+2)=Fcn(1,2);
    FF1(6*(cn-1)+3)=Fcn(1,3);
    FF1(6*(cn-1)+4)=Fcn(2,2);
    FF1(6*(cn-1)+5)=Fcn(2,3);
    FF1(6*(cn-1)+6)=Fcn(3,3);
% modified
     
end
%    FF1(6*cellnumber+1)=dtsigma(3,3);


end