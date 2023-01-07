% for DTLZ1
function z=DTLZ1(x,nObj)
    nVar=numel(x);
    i=nObj:nVar;
    g=sum(power(x(i)-0.5,2)-cos(20*pi*(x(i)-0.5)));
    g=100*(numel(i)+g);
    
    z=zeros(nObj,1);
    
    for k=1:nObj
        z(k)=(1+g)*0.5;
    end
    
    for k=1:nObj
        for m=1:nObj-k
            z(k)=z(k)*x(m);
        end
        if k~=1
            z(k)=z(k)*(1-x(nObj-k+1));
        end
    end
    
end