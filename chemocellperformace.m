clear;
cellr=0.3;
thick=0.1;
cellnumber=10; % number of cells in REV
timesteps=400000; % creep time steps
vsolid=4/3*(cellr^3-(cellr-thick)^3);
time=1:timesteps;
sigma=[10 0 0; 0 5 0; 0 0 5]; 


[nrr,vce,dfi,inai,tvai,tvbi,tvci]=chemocell(sigma,cellr,thick,vsolid,timesteps);
for i=1:timesteps
    ina(i)=inai(i);
    tva(i)=tvai(i);
    tvb(i)=tvbi(i);
    tvc(i)=tvci(i);

end
fi=(ina.^3)./vsolid;

figure('Name','Cell strain','NumberTitle','off')
    zz=plot(time,tva./cellr,'r');
    hold on 
    xx=plot(time,tvb./cellr,'g');
    hold on 
    yy=plot(time,tvc./cellr,'b');
    
        xlabel('Time(s)','fontsize',16)
        ylabel('Cell strain','fontsize',16)
      legend([zz,xx,yy],'z direction','x direction','y direction','Location','northwest')

      figure('Name','Porosity','NumberTitle','off')      
       pp=plot(time,fi,'r');   
           xlabel('Time(s)','fontsize',16)
        ylabel('Porosity','fontsize',16)
      
      figure('Name','Void radius','NumberTitle','off')      
       pp=plot(time,ina,'r');   
           xlabel('Time(s)','fontsize',16)
        ylabel('Void radius','fontsize',16)