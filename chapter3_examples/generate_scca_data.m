function [X,Y] = generate_scca_data

n = 50;
p = 100;
q = 150;

for i = 1:p
    X(:,i) = normrnd(0,1,[n,1]);
end

for j = 1:q
    Y(:,j) = normrnd(0,1,[n,1]);
end




end
