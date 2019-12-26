function SS=sphereS 
v=0.3;
SS(1:3,1:3,1:3,1:3)=0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                dij=0; dkl=0;dik=0;djl=0;dil=0;djk=0;
                
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
                SS(i,j,k,l)=(5*v-1)/15/(1-v)*dij*dkl+(4-5*v)/15/(1-v)*(dik*djl+dil*djk);
                
                
            end
        end
    end
end
end