function   fdef=fourdotef(A, B)
fdef(1:3,1:3,1:3,1:3)=0;
for i=1:3
    for j=1:3
       for k=1:3
           for l=1:3
                for m=1:3
                    for n=1:3
                        for o=1:3
                            for p=1:3
                              fdef(i,j,k,l)=fdef(i,j,k,l)+A(i,j,k,l,m,n,o,p)*B(p,o,n,m);
                            end
                        end
                    end
                end
           end
       end
    end
end
end