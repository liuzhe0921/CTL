function mass = DS_combination(mass)

mass_com = mass(1,:);
for i = 1:size(mass, 1) - 1
    for j = 1:size(mass, 2) - 1
         mass_com(1,j) = mass_com(1, j) * (mass(i + 1, j) + mass(i + 1, end)) + mass(i + 1, j) * mass(i, end);
    end
    mass_com(1, end) = mass_com(1, end) * mass(i + 1, end);
end
m_k = sum(mass_com);

mass = mass_com ./ m_k;
