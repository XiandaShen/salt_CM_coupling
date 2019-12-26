
function SH=sphereH (muo,fi)
mu=muo*fi;
v=0.3;
SH(1:3,1:3,1:3,1:3)=0;
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
                SH(i,j,k,l)=mu/(4-5*v)*((3-5*v)*dij*dkl+(7-5*v)/2*(dik*djl+dil*djk));
                
                
            end
        end
    end
end
end