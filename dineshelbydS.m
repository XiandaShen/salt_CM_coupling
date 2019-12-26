function dinSdS=dineshelbydS(inS)

dinSdS(1:3,1:3,1:3,1:3,1:3,1:3,1:3,1:3)=0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    for n=1:3
                        for o=1:3
                            for p=1:3                           
                              dinSdS(i,j,k,l,m,n,o,p)=-inS(i,j,k,l)*inS(m,n,o,p);                               
                            end
                        end
                    end
                end
            end
        end
    end
end
end