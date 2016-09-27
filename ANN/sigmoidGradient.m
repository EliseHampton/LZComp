function g = sigmoidGradient(z)
% SIGMOIDGRADIENT returns the gradient of the sigmoid function
% evaluated at z
% g = SIGMOIDGRADIENT(z) computes the gradient of the sigmoid function
% evaluated at z. 

g = zeros(size(z));

temp1 = sigmoid(z);

g = temp1.*(1-temp1);

end
