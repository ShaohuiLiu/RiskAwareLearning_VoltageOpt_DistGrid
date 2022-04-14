% generate incidence matrix and vector of line admittances
function [Mtilde, y] = nodeEdgeIncidence(Ybus)
    [n,temp] = size(Ybus);
    L = n-1;
    Mtilde = zeros(n,L);
    y = zeros(L,1);
    
    ell = 1;
    for i = 1:n-1
       for j = i+1:n
          if Ybus(i,j) ~= 0
             Mtilde(i,ell) = 1;
             Mtilde(j,ell) = -1;
             y(ell) = -Ybus(i,j);
             ell = ell + 1;
          end
       end
    end
end