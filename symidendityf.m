function R=symidendityf
R(1:3,1:3,1:3,1:3) = 0;
R(1,1,1,1)=1;
R(2,2,2,2)=1;
R(3,3,3,3)=1;
R(2,3,2,3)=0.5;
R(2,3,3,2)=0.5;
R(3,2,3,2)=0.5;
R(3,2,2,3)=0.5;
R(2,1,2,1)=0.5;
R(2,1,1,2)=0.5;
R(1,2,2,1)=0.5;
R(1,2,1,2)=0.5;
R(1,3,1,3)=0.5;
R(1,3,3,1)=0.5;
R(3,1,1,3)=0.5;
R(3,1,3,1)=0.5;
end

