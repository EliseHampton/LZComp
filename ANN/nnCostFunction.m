function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
% NNCOSTFUNCTION Implements the neural network cost function for a
% two hidden layer neural network which performs classification
% [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
% X, y, lambda) computes the cost and gradient of the neural network. The
% parameters for the neural network are "unrolled" into the vector
% nn_params and need to be converted back into the weight matrices. 
%

% Reshape nn_params back into the parameters Theta1, Theta2, and Theta3
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):(((hidden_layer_size * (input_layer_size + 1)) + (hidden_layer_size*(1+hidden_layer_size))))), ...
                 hidden_layer_size, (hidden_layer_size + 1));

Theta3 = reshape(nn_params((1 + ((hidden_layer_size * (input_layer_size + 1)) + (hidden_layer_size*(1+hidden_layer_size)))):end), ...
                 num_labels, (hidden_layer_size + 1));

% number of examples
m = size(X, 1);
         
% values to return 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));
Theta3_grad = zeros(size(Theta3));

y_matrix = eye(num_labels)(y,:);

%add X0's
X = [ones(m, 1) X];

% Compute first layer
[z2] = Theta1*X';
[a2] = sigmoid(z2);
%add a0's
a2 = [ones(1, m); a2];
% Compute second layer
[z3] = Theta2*a2;
[a3] = sigmoid(z3);
a3 = [ones(1, m); a3];
% Compute third layer or output layer
[z4] = Theta3*a3;
[a4] = sigmoid(z4);

%calculate cost function
temp1 = -y_matrix.*log(a4)';
temp2 = -(1-y_matrix).*log(1-a4)';
J = sum((1/m)*sum(temp1+temp2));

%sum theta1 and 2 and add to J
sum_theta1 = sum(sum(Theta1(:,2:size(Theta1,2)).^2));
sum_theta2 = sum(sum(Theta2(:,2:size(Theta2,2)).^2));
sum_theta3 = sum(sum(Theta3(:,2:size(Theta3,2)).^2));

J = J + (lambda/(2*m))*(sum_theta1 + sum_theta2 + sum_theta3);

% Calculate the gradient of the cost function
delta_4 = a4 - y_matrix';
Delta_3 = delta_4*a3';

temp = sigmoidGradient(z3);
delta_3 = Theta3(:,2:size(Theta3,2))'*delta_4.*temp;
Delta_2 = delta_3*a2';

temp = sigmoidGradient(z2);
delta_2 = Theta2(:,2:size(Theta2,2))'*delta_3.*temp;
Delta_1 = delta_2*X;

Theta1_grad = (1/m)*Delta_1;
Theta2_grad = (1/m)*Delta_2;
Theta3_grad = (1/m)*Delta_3;



Theta1_grad(:,2:size(Theta1_grad,2)) = Theta1_grad(:,2:size(Theta1_grad,2)) ...
    + (lambda/m)*Theta1(:,2:size(Theta1_grad,2));
Theta2_grad(:,2:size(Theta2_grad,2)) = Theta2_grad(:,2:size(Theta2_grad,2)) ...
    + (lambda/m)*Theta2(:,2:size(Theta2_grad,2));
Theta3_grad(:,2:size(Theta3_grad,2)) = Theta3_grad(:,2:size(Theta3_grad,2)) ...
    + (lambda/m)*Theta3(:,2:size(Theta3_grad,2));


% Unroll gradients for returning
grad = [Theta1_grad(:) ; Theta2_grad(:) ; Theta3_grad(:)];


end
