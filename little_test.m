clear;
sss=eshelby(0.3,0.2999,0.2998);
H=doubledotff(stiffness(0),inversegeneral(sss))-stiffness(0);

HH(1:3,1:3,1:3,1:3) = 0;
muo=8846;
v=0.3;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                dij=0;
                dkl=0;
                dik=0;
                djl=0;
                dil=0;
                djk=0;
               
              if i==j
                  dij=1;
              end
              
              if k==l
                  dkl=1;
              end
              
              if i==k
                  dik=1;
              end
              
              if j==l
                  djl=1;
              end
              
              if i==l
                  dil=1;
              end
              
              if j==k
                  djk=1;
              end
              
              HH(i,j,k,l)=muo/(4-5*v)*((3-5*v)*dij*dkl+(7-5*v)/2*(dik*djl+dil*djk));
                
            end
        end
    end
end
