% for DTLZ3
function z=DTLZ3(x,nObj)

    nVar=numel(x);
    i=nObj:nVar;
    g=sum(power(x(i)-0.5,2)-cos(20*pi*(x(i)-0.5)));
    g=100*(numel(i)+g);
    
    z=zeros(nObj,1);
    
    for k=1:nObj
        z(k)=1+g;
    end
    
    for k=1:nObj
        for m=1:nObj-k
            z(k)=z(k)*cos(0.5*pi*x(m));
        end
        if k~=1
            z(k)=z(k)*sin(0.5*pi*x(nObj-k+1));
        end
    end

end