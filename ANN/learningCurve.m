function [error_train, error_val] = ...
    learningCurve(X, y, Xval, yval, lambda, init_nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels)
% LEARNINGCURVE Generates the train and cross validation set errors needed 
% to plot a learning curve [error_train, error_val] = ...
% LEARNINGCURVE(X, y, Xval, yval, lambda,init_nn_params, input_layer_size, ...
% hidden_layer_size, num_labels) returns the train and
% cross validation set errors for a learning curve. In particular, 
% it returns two vectors of the same length - error_train and 
% error_val. Then, error_train(i) contains the training error for
% i examples (and similarly for error_val(i)).


% Number of training examples
m = size(X, 1);

% values to return
error_train = zeros(m/100, 1);
error_val   = zeros(m/100, 1);


for i = 100:100:m
    % Compute train/cross validation errors using training examples 
    
    options = optimset('MaxIter', 1000);
    lambda = 3;
    costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, X(1:i,:), y(1:i), lambda);
    [nn_params, cost] = fmincg(costFunction, init_nn_params, ...
                                options);

    % now want to compute the cost on the test and CV set with
    
    tmp1 = CostFunctionT(nn_params, ...
                         input_layer_size, ...
                         hidden_layer_size, ...
                         num_labels, X(1:i,:), y(1:i), 0);
    error_train(i/100) = tmp1;
    tmp2 = CostFunctionT(nn_params, ...
                         input_layer_size, ...
                         hidden_layer_size, ...
                         num_labels, Xval, yval, 0);
    error_val(i/100) = tmp2;

           
end




% -------------------------------------------------------------

% =========================================================================

end
