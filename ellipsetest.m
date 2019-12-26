gsigma=[3.1 0 0; 0 3.1 0; 0 0 3.1];
casenumber=9;
%% orientation
psi=[0; pi/6; pi/3];
theta=[0; pi/6; pi/3];
phi=[0; pi/6; pi/3];


% orientaion of each cell
%  for i=1:cellnumber
% % % phi(i)=mod(i,5);
%  theta(i)=pi/2/cellnumber*fix(i-1);
% % 
% % 
%  end
Q=cell(1,27); % transformation matrix
countq=0;
for i=1:3
for    j=1:3
for    k=1:3
    countq=countq+1;
    Q(countq)={transmatrix(psi(i),theta(j),phi(k))};
  %    Q(i)={eye(3)};  % no orientation
%    Q{i}*Q{i}'
end
end
end   



totoalt=60000;
interval=100;


for cnum=1:27

ga=zeros(1,totoalt);
gb=zeros(1,totoalt);
gc=zeros(1,totoalt);
gr=zeros(1,totoalt);
lsigma=Q{cnum}*gsigma*Q{cnum}';
ga(1)=0.3;
gb(1)=0.25;
gc(1)=0.2;
gr(1)=0.05;
vce=zeros(1,totoalt);
gvce=zeros(1,totoalt);
fi=zeros(1,totoalt);
fi(1)=(ga(1)-2*gr(1))*(b(1)-2*gr(1))*(gc(1)-2*gr(1))/ga(1)/gb(1)/gc(1);
for time=1:totoalt

[ga(time+1),gb(time+1),gc(time+1),gr(time+1),vce(time+1),fi(time+1)]=ellipcell(lsigma,ga(time),gb(time),gc(time),gr(time),interval);

gvce(time+1)=Q{cnum}'*vce(time+1)*Q{cnum};


end
save(['C:\Users\xshen60\Downloads\ellipse cell\oi',num2str(cnum),'a.mat'],'ga');
save(['C:\Users\xshen60\Downloads\ellipse cell\oi',num2str(cnum),'b.mat'],'gb');
save(['C:\Users\xshen60\Downloads\ellipse cell\oi',num2str(cnum),'c.mat'],'gc');
save(['C:\Users\xshen60\Downloads\ellipse cell\oi',num2str(cnum),'r.mat'],'gr');
save(['C:\Users\xshen60\Downloads\ellipse cell\oi',num2str(cnum),'vce.mat'],'vce');
save(['C:\Users\xshen60\Downloads\ellipse cell\oi',num2str(cnum),'fi.mat'],'fi');
end

