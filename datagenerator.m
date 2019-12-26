cellnumber=10;
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

plot (iniX)
mean (iniX)
var(iniX)
thick=cellr-iniX;