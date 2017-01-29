function [X,Y] = generate_kernel_data2(n)

%n = 150; %hardoon 10000
p = 7;
q = 8;

for i = 1:p
    X(:,i) = normrnd(0,1,[n,1]);
end

y1 = exp(X(:,3)) + normrnd(0,0.4,[n,1]);
y2 = X(:,1).^3 + normrnd(0,0.2,[n,1]);
y3 = -X(:,4) + normrnd(0,0.3,[n,1]);

for j = 1:q-3
    Y2(:,j) = normrnd(0,1,[n,1]);
end
Y = [y1 y2 y3 Y2];

end
