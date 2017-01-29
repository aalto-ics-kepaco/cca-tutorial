function [X,Y] = generate_data_pmd

n = 50;
p = 100;
q = 150;

for i = 1:p
    X(:,i) = normrnd(0,1,[n,1]);
end

y1 = X(:,3) + normrnd(0,0.08,[n 1]);
y2 = X(:,1) + normrnd(0,0.07,[n 1]);
y3 = -X(:,4) + normrnd(0,0.05,[n 1]); 

for j = 1:q-3
    Y2(:,j) = normrnd(0,1,[n,1]);
end
Y = [y1 y2 y3 Y2];

end
