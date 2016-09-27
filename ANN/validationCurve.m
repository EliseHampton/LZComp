function [lambda_vec, error_train, error_val,error_test] = ...
    validationCurve(X, y, Xval, yval,Xtest,ytest,init_nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels)
% VALIDATIONCURVE Generate the train and validation errors needed to
% plot a validation curve to select best lambda
% [lambda_vec, error_train, error_val] = ...
% VALIDATIONCURVE(X, y, Xval, yval,Xtest,ytest,init_nn_params, ...
% input_layer_size, hidden_layer_size, num_labels ) returns the train
% and validation errors (in error_train, error_val)
% for different values of lambda. 
%

% Selected values of lambda for testing! don't change this
lambda_vec = [0 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10]';

% variables to return
error_train = zeros(length(lambda_vec), 1);
error_val = zeros(length(lambda_vec), 1);


%loop through the different lambdas and calculate the cost function
for i = 1:length(lambda_vec)
    lambda = lambda_vec(i);
    
    options = optimset('MaxIter', 1000);
    
    costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, X, y, lambda);
    [nn_params, cost] = fmincg(costFunction, init_nn_params, options);
    
    tmp1 = CostFunctionT(nn_params, ...
                         input_layer_size, ...
                         hidden_layer_size, ...
                         num_labels, X, y, 0);
    error_train(i) = tmp1;
    tmp2 = CostFunctionT(nn_params, ...
                         input_layer_size, ...
                         hidden_layer_size, ...
                         num_labels, Xval, yval, 0);
    error_val(i) = tmp2;
    
    tmp3 = CostFunctionT(nn_params, ...
                         input_layer_size, ...
                         hidden_layer_size, ...
                         num_labels, Xtest, ytest, lambda);
    error_test(i) = tmp3;
       
end











% =========================================================================

end
