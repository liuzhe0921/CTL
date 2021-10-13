%------------------------------------------------------------------------------------------------------------------------
% Contributed by Zhe Liu
% Ref:
% Credal Transfer Learning with Multi-estimation for Missing Data (IEEE ACCESS 2020).
%------------------------------------------------------------------------------------------------------------------------
function [I, D]= Find_KNN(data_com, data_inc, KK)

[x, y] = size(data_inc);    
[x1, ~] = size(data_com);   

a1 = zeros(0);
for i = 1:x
    t = 1;
    for j = 1:y
        if isnan(data_inc(i, j))
            a1(i, t) = j;  
            A2 = a1;   
            t = t + 1;
        end
    end    
end
 
I = zeros(0);
D = zeros(0);
for i = 1:x
    a1 =A2;
    A = data_inc(i,:);      
    B = data_com(1:x1,:);   
    t2 = 1;
    a2 = zeros(0);
    for j = 1:size(a1, 2)     
        if a1(i,j) == 0
            a2(t2,1) = j;
            t2 = t2 + 1;
        end
    end
    a1(:,a2) =[];
    B(:,a1(i,:)) = [];  
    A(:,a1(i,:)) = [];  
    dist=zeros(x1,1);
    for i2=1:x1
        dist(i2,:)=norm(B(i2,:)-A);
    end
    [Y,I2]=sort(dist);   
    I{i, 1} = I2(1:KK,:);     
    D{i, 1} = Y(1:KK,:);    
end
end
