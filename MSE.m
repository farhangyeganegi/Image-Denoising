function y = MSE( A, B );
% Description: This code calculates the normalized mean-square error
% between two images.
% Input
% A: input image one.
% B: input image two.
% A and B should have the same size.

l = size(A);
k = prod(l);
y = 0;
for i=1:l(1)
    for j=1:l(2)
        y = y + ((A(i,j)-B(i,j))^2)/k;
    end
end
end

