function I=majorsymidendityeighth




I(1:3,1:3,1:3,1:3,1:3,1:3,1:3,1:3)=0;
for i=1:3
    for j=1:3
       for k=1:3
           for l=1:3
                for m=1:3
                    for n=1:3
                        for o=1:3
                            for p=1:3
dim=0;
djn=0;
din=0;
djm=0;
dko=0;
dlp=0;
dkp=0;
dlo=0;
dio=0;
djp=0;
dip=0;
djo=0;
dkm=0;
dln=0;
dkn=0;
dlm=0;

                              if i==m
                                  dim=1;
                              end
                              if j==n
                                  djn=1;
                              end
                              if i==n
                                  din=1;
                              end                              
                              if j==m
                                  djm=1;
                              end                                
                              if k==o
                                  dko=1;
                              end                                
                              if l==p
                                  dlp=1;
                              end
                              if k==p
                                  dkp=1;
                              end                              
                              if l==o
                                  dlo=1;
                              end 
                              if i==o
                                  dio=1;
                              end
                              if i==p
                                  dip=1;
                              end
                              if j==o
                                  djo=1;
                              end
                              if j==p
                                  djp=1;
                              end
                              if k==m
                                  dkm=1;
                              end
                              if k==n
                                  dkn=1;
                              end
                              if l==m
                                  dlm=1;
                              end
                              if l==n
                                  dln=1;
                              end
                              I(i,j,k,l,m,n,o,p)=1/4*(dim*djn+din*djm)*(dko*dlp+dkp*dlo)+1/4*(dio*djp+dip*djo)*(dkm*dln+dkn*dlm);
                            end
                        end
                    end
                end
           end
       end
    end
end
end

        

