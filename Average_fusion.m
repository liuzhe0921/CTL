function mass = Average_fusion(mass)

mass_com = mass(1,:);
for i = 1:size(mass, 1) - 1
    for j = 1:size(mass, 2) 
         mass_com(1,j) = mass_com(1, j) * mass(i + 1, j);
    end
end
m_k = sum(mass_com);

mass = mass_com ./ m_k;
