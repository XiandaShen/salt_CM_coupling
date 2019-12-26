function [QQ]=transmatrix(psiQ, thetaQ, phiQ)
QQ=[cos(psiQ)*cos(thetaQ)*cos(phiQ)-sin(psiQ)*sin(phiQ)      sin(psiQ)*cos(thetaQ)*cos(phiQ)+cos(psiQ)*sin(phiQ)   -sin(thetaQ)*cos(phiQ)
    -cos(psiQ)*cos(thetaQ)*sin(phiQ)-sin(psiQ)*cos(phiQ)     -sin(psiQ)*cos(thetaQ)*sin(phiQ)+cos(psiQ)*cos(phiQ)  sin(thetaQ)*sin(phiQ)           
    cos(psiQ)*sin(thetaQ)                                  sin(psiQ)*sin(thetaQ)                  cos(thetaQ)];
end
