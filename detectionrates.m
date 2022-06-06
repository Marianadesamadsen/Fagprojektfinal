function [truepositive,falsepositive,falsenegative] = detectionrates(stride,D,D_detected,Ts)

falsenegative   = 0;
falsepositive   = 0;
truepositive    = 0;
truemealdetec   = zeros(1,length(D(1,:)));
tempmealdetec   = zeros(1,length(D_detected));


for i = 1:length(D(1,:))
    
    if D(1,i) >= 50/Ts
        D(1,i) = 1;
        
    else
        D(1,i) = 0;
    end 
    
end 


for i = 1 : length(D(1,:))
    if D(1,i) == 1
        truemealdetec(i) = i;
    end 
end

for i = 1 : length(D(1,:))
    if D_detected(i) == 1
        tempmealdetec(i) = i;
    end 
end

idxdetecmeals = find(tempmealdetec);
idxtruemeals = find(truemealdetec);


for i = 1:length(idxtruemeals)
   
    k = idxtruemeals(i);
    
    if sum(D_detected(k:k+stride)) == 0 
        falsenegative = falsenegative + 1;
    end  
end

for i = 1:length(idxdetecmeals)
   
    k = idxdetecmeals(i);
    
    if sum(D(1,k-stride:k)) == 0 
        falsepositive = falsepositive + 1;
    end  
end

for i = 1:length(idxdetecmeals)
   
    k = idxdetecmeals(i);
    
    if sum(D(1,k-stride:k)) == 1 
        truepositive = truepositive + 1;
    end  
end


end