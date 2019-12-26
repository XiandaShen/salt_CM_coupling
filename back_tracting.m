% function [alpha, iteralpha] = back_tracting( iteralpha, d, f, x, j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,dinter)
function [alpha, iteralpha] = back_tracting( iteralpha, d, f, x, creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,muo) 
alpha = 1;
rho = 0.8;
c = 1e-3;
 insidec=0;
condition_hold = 0;
while ~condition_hold && insidec<10;
    xalpha = x + alpha * d;
%    [falpha,totalp,tepsilon,dinter,dsigma,depsilon,epsilon,cella,cellb,cellc,thick]=puremechanics(j,itnum,cellnumber,cella,cellb,cellc,thick,dsigma,dtsigma,depsilon,epsilon,tepsilon,totalp,ditner,xalpha);
    [falpha,totalp,tepsilon,dinter,dsigma,sigma,depsilon,epsilon,epsilone,cellr,thick]=encalibrationchemomechanics(creepinterval,j,itnum,cellnumber,cellr,Q,thick,dHdL,dsigma,dtsigma,sigma,depsilon,epsilon,tepsilon,epsilone,depsilone,dvcepsilon,totalp,dinter,xalpha,muo);


    rhs = f' * f*(1 - c * alpha);

    
    
    
       insidec=insidec+1;
    iteralpha=iteralpha+1;
   
    
%     iteralpha
%     falpha' * falpha 
    
    
    
    if (falpha' * falpha <= rhs)
        condition_hold = 1;

    else
        alpha = rho * alpha;
    end
    
%     if insidec==20
%         alpha=-1;
%     end
end