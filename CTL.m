%------------------------------------------------------------------------------------------------------------------------
% Contributed by Zhe Liu
% Ref:
% Credal Transfer Learning with Multi-estimation for Missing Data (IEEE ACCESS 2021).
%------------------------------------------------------------------------------------------------------------------------
clear
clc

%% input dataset (source domain with complete data and target domain with incomplete data)
load ionosphere23.mat

num_K = 5;  % KNNs
phi = 0.1;  % meta-cluster threshold
% 5% one-to-one patterns pairs
n = round(length(D)*0.05);  
D_bridge = D(1:n, :);
T_bridge = T(1:n, :);
D(1:n, :) = [];  % source domain
T(1:n, :) = [];  % target domain

%% multi-estimation
[index, dist] = Find_KNN(T_bridge(:, 1:end-1), T(:, 1:end-1), num_K);  

D_KNN = zeros(0);
for i = 1:size(index, 1)
    DC_data = D_bridge(index{i, 1}, :);  
    D_KNN{i, 1} = DC_data;
end

multi_estimates = [];
for i = 1:size(D_KNN, 1)
    multi_estimates = [multi_estimates; D_KNN{i, 1}];
end

%% classification
test_data = [D_bridge(:, 1:end-1); multi_estimates(:, 1:end-1)];
test_label = [D_bridge(:, end); T(:, end)];

[L, k]=Knnclassify(test_data, D(:, 1:end-1), D(:,end), 5,'euclidean','nearest');

L1 = L(1:size(D_bridge, 1), :);
k1 = k(1:size(D_bridge, 1), :);
L(1:size(D_bridge, 1), :) = [];
k(1:size(D_bridge, 1), :) = [];
LL=[]; kk = [];
for i = 1:num_K:size(L, 1)
    L2 = zeros(0); k2 = zeros(0);
    k2{1, 1} = k(i:i+num_K-1, :);
    L2{1, 1} = L(i:i+num_K-1, :);
    LL = [LL; L2];
    kk = [kk; k2];
end

%% discounted classification results and global fusion
w = zeros(0);
for i = 1:size(dist,1)   
    D3 = [];
    sum1 = 0;
    for j = 1:num_K
        D1 = exp(-dist{i, 1}(j,1));
        sum1 = sum1 + D1;  
        D3 = [D3; D1];
    end
    for j = 1:num_K
        D2 = exp(-dist{i, 1}(j,1))/max(D3);
        w(i,j) = D2;
    end
end

discount = zeros(0);
mass = zeros(0);
for i = 1:size(D, 1)
    in = zeros(0); mass1 = zeros(0); W0 = zeros(0); 
    discount =  LL{i, 1};  
    m_conflict = repmat(w(i, :)', 1, size(discount, 2)) .* discount;  
    mvide = ones(size(m_conflict, 1), 1) - sum(m_conflict, 2);
    mass_conflict = [m_conflict mvide];  
    
    
    for j = 1:max(T(:,end))
        in{j, 1} = find(kk{i, 1}(:, 1) == j);
    end
    t = 1; p = 1; a = zeros(0); b = zeros(0);
    for j = 1:size(in, 1)
        if length(in{j, 1}) > 1
            a(t, 1) = j;
            t = t + 1;
        elseif length(in{j, 1}) == 1
            b(p, 1) = j;
            p = p + 1;
        end
    end
    for j = 1:size(a, 1)
        mass1(j, :) = Average_fusion(mass_conflict(in{a(j, 1), 1}, :));
    end
    if isempty(b)
        max_conflict = mass1;
    else
        max_conflict = [mass1; mass_conflict([in{b, 1}], :)];
    end
    
    if size(max_conflict, 1) == 1
        mass(i, :) = max_conflict;
    else
        max_mass_conflict = zeros(0);
        col = zeros(0);
        for j = 1:size(max_conflict, 1)
            [max_mass_conflict(j, 1), col(j, 1)] = max(max_conflict(j, 1:end-1));
        end
        [max_mass_sort, Q] = sort(max_mass_conflict, 'descend');  
        if abs(max_mass_sort(1) - max_mass_sort(2)) <= phi
            mass_com = max_conflict(1, :);
            for k = 1:size(max_conflict, 1) - 1
                for j = 1:size(max_conflict, 2)
                    if j < size(max_conflict, 2)
                        mass_com(1, j) = mass_com(1, j) * (max_conflict(k+1, j)+max_conflict(k+1, end)) + mass_com(1, end) * max_conflict(k+1, j);
                    else
                        mass_com(1, j) = (max_mass_sort(1) * max_mass_sort(2));
                    end  
                end
            end
            mass(i, :) = mass_com ./ sum(mass_com);
        else
            mass(i, :) = DS_combination(max_conflict); 
%             mass_com = max_conflict(1, :);
%             for k = 1:size(max_conflict, 1) - 1
%                 for j = 1:size(max_conflict, 2)
%                     if j < size(max_conflict, 2)
%                         mass_com(1, j) = mass_com(1, j) * (max_conflict(k+1, j)+max_conflict(k+1, end)) + mass_com(1, end) * max_conflict(k+1, j);
%                     else
%                         mass_com(1, j) = mass_com(1, end) * max_conflict(k + 1, end);
%                     end  
%                 end
%             end
%             mass(i, :) = mass_com ./ sum(mass_com);
        end          
    end
end

%% error rate and imprecision
mass = mass';
[~, I] = max(mass);
I = I';
in0 = find(I == max(T(:, end))+1);
I(in0, :) = [];
data_te1 = T(:, end);
data_te1(in0, :)=[];
new_label = [k1; I];
true_label = [test_label(1:size(D_bridge, 1), :); data_te1];
err = 0;
for i = 1:size(new_label, 1)
    if new_label(i) ~= true_label(i, end)
        err = err + 1;
    end
end

err1 = err*100;
result(1, 1) = err1/length(test_label);  % error rate
result(1, 2) = length(in0)*100/length(test_label);  % imprecision rate










