function   ddfe=doubledotfe(A, B)
ddfe(1:3,1:3,1:3,1:3,1:3,1:3,1:3,1:3)=0;
for i=1:3
    for j=1:3
       for k=1:3
           for l=1:3
                for m=1:3
                    for n=1:3
                        for o=1:3
                            for p=1:3
                                for a=1:3
                                    for b=1:3
                              ddfe(i,j,k,l,m,n,o,p)=ddfe(i,j,k,l,m,n,o,p)+A(i,j,b,a)*B(a,b,k,l,m,n,o,p);
                                    end
                                end
                            end
                        end
                    end
                end
           end
       end
    end
end
end