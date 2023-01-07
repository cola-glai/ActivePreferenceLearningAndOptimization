% for DTLZ4
function z=DTLZ4(x,nObj)

    alpha=100;
    nVar=numel(x);
    i=nObj:nVar;
    g=sum(power(x(i)-0.5,2));
    
    z=zeros(nObj,1);
    
    for k=1:nObj
        z(k)=1+g;
    end
    
    for k=1:nObj
        for m=1:nObj-k
            z(k)=z(k)*cos(0.5*pi*power(x(m),alpha));
        end
        if k~=1
            z(k)=z(k)*sin(0.5*pi*power(x(nObj-k+1),alpha));
        end
    end

end