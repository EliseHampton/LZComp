function [J] = CostFunctionT(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
% COSTFUNCTIONT Implements the neural network cost function for a
% two hidden layer neural network which performs classification
% [J] = COSTFUNCTONT(nn_params, input_layer_size, hidden_layer_size, num_labels, ...
% X, y, lambda) computes the cost of the neural network. 

% Reshape nn_params back into the parameters Theta1, Theta2, and Theta3
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):(((hidden_layer_size * (input_layer_size + 1)) + (hidden_layer_size*(1+hidden_layer_size))))), ...
                 hidden_layer_size, (hidden_layer_size + 1));

Theta3 = reshape(nn_params((1 + ((hidden_layer_size * (input_layer_size + 1)) + (hidden_layer_size*(1+hidden_layer_size)))):end), ...
                 num_labels, (hidden_layer_size + 1));

% number of examples
m = size(X, 1);
         
% values to be returned
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));
Theta3_grad = zeros(size(Theta3));

%implement the changes in y to be vectors
y_matrix = eye(num_labels)(y,:);

%add X0's
X = [ones(m, 1) X];
% Compute second layer
[z2] = Theta1*X';
[a2] = sigmoid(z2);
%add a0's
a2 = [ones(1, m); a2];
% Compute third layer
[z3] = Theta2*a2;
[a3] = sigmoid(z3);
% Compute output layer
a3 = [ones(1, m); a3];
[z4] = Theta3*a3;
[a4] = sigmoid(z4);

% Calculate cost function
temp1 = -y_matrix.*log(a4)';
temp2 = -(1-y_matrix).*log(1-a4)';
J = sum((1/m)*sum(temp1+temp2));

%sum theta1, theta2, and theta3 and add to J with lambda/2m
sum_theta1 = sum(sum(Theta1(:,2:size(Theta1,2)).^2));
sum_theta2 = sum(sum(Theta2(:,2:size(Theta2,2)).^2));
sum_theta3 = sum(sum(Theta3(:,2:size(Theta3,2)).^2));

J = J + (lambda/(2*m))*(sum_theta1 + sum_theta2 + sum_theta3);



end
