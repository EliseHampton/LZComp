#! /Applications/Octave.app/Contents/Resources/bin/octave -qf
%% ============= Welcome to lzcomp_ann ================= %%
% classifier is the neural net that is used by LZClassify.
% LZClassify is a program to determine the number of components
% that are most probable as output by LZIFU (Ho et al. 2016b)
% Tclassifier uses 2 hidden layers and 15 nodes at each layer
% to calculate the probability that a certain spaxel should 
% be fit with 1, 2, or 3 components.
% classifier is writen in Octave and uses a sigmoid function.
% The minimisation function was originally created by Carl Ramussen
% and adapted by Andrew Ng for the online MOOC in Machine Learning.

% LZComp was writen by Elise Hampton at the Australian National
% University as a part of her PhD (2013-2016) for the SAMI and S7 
% surveys.
% Contact Details: elise.hampton@anu.edu.au

%% ================== Update Log ===================== %%
% Dec 2014 - Testing of ClockWorkv1 completed (83% accuracy)
% Dec 2014 - Rewrite to TheMachine
% Dec 2014 - Test of TheMachine - all good
% Jan 2015 - Comment and add in not-training parameters
% Jan 2015 - train and test variables added; use to tell when want
% to test or train using TheMachine
% 13th Jan 2015 - now outputs predicted values to text files as
% lists
% 19th Jan 2015 - Takes in train and test as command line arguments
% 21st Jul 2016 - Renamed to be better identified when published
% renamed to lzcomp_ann.m

%% Initialization
clear ; close all; clc
arg_list = argv();
train = double(arg_list{4});
testt = double(arg_list{5});

% Read in text files created by createInput.py
X2 = dlmread('input.txt',' ');
m = size(X2, 1);

% if in training mode, also read in cross validation text file and
% test text file
if train == 49
  yay = 'yay'
  y2 = dlmread('output.txt',' ');
  Xval = dlmread('input_cv.txt',' ');
  yval = dlmread('output_cv.txt',' ');
  Xtest = dlmread('input_test.txt',' ');
  ytest = dlmread('output_test.txt',' ');
  ym = size(y2, 1);
  mcv = size(Xval, 1);
  mtest = size(Xtest, 1);
  ymcv = size(yval, 1);
  ymtest = size(ytest, 1);
    
end
 
    
%% Input layers and hidden layers sizes

input_layer_size = size(X2,2);
hidden_layer_size = 15;  % 15 hidden units
num_labels = 3;          % 3 labels, from 1 to 3   
                        

%% initialise the theta values for the 3 hidden layers if in training mode
if train == 49
    
    initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
    initial_Theta2 = randInitializeWeights(hidden_layer_size, hidden_layer_size);
    initial_Theta3 = randInitializeWeights(hidden_layer_size, num_labels);

    % Unroll parameters
    initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:) ; initial_Theta3(:)];

    fprintf('\nTraining Neural Network... \n')

    % the number of iterations can be changed here 50,000 is a decent ball park
    % takes about 20-30 minutes
    options = optimset('MaxIter', 50000);

    %  lambda can be changed here - for SAMI and S7 lambda - 3.0 is best, make sure this is the same value as used when trained so that the results are correct
    lambda = 3.0;

    % for the cost function to be minimized
    costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, X2, y2, lambda);

    % Training happens here
    [nn_params, cost] = fmincg(costFunction, initial_nn_params, options);

    % Obtain Theta1, Theta2, and Theta3 from the training
    Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                     hidden_layer_size, (input_layer_size + 1));

    Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):(((hidden_layer_size * (input_layer_size + 1)) + (hidden_layer_size*(1+hidden_layer_size))))), hidden_layer_size, (hidden_layer_size + 1));

    Theta3 = reshape(nn_params((1 + ((hidden_layer_size * (input_layer_size + 1)) + (hidden_layer_size*(1+hidden_layer_size)))):end), ...
                     num_labels, (hidden_layer_size + 1));

    % Save the Thetas into files that can be accessed later for running
    save theta1.mat Theta1
    save theta2.mat Theta2
    save theta3.mat Theta3
    
end
    
%% If also testing then try out these plots
% It will take a while to do the testing

if testt == 49
    
    % Test how the nueral net goes with different number of
    % training examples
    lambda = 0;
    [error_train, error_val] = ...
        learningCurve(X2, y2, ...
                      Xval, yval, ...
                      lambda,initial_nn_params, ...
                      input_layer_size, ...
                      hidden_layer_size, ...
                      num_labels);
    plot(1:size(error_train,1), error_train, 1:size(error_val,1), error_val);
    title('Learning curve for linear regression')
    legend('Train', 'Cross Validation')
    xlabel('Number of training examples *100')
    ylabel('Error')
    axis([0 (size(error_train,1)+1) 0 12])
    fprintf('Program paused. Press enter to continue.\n');
    %print -deps trainexampls.eps;
    pause; %pause allows the plot to actually show up, if no pause it won't show properly
    
    % Test how the nueral net goes with different regularisation parameters, lambda
    [lambda_vec, error_train, error_val, error_test] = ...
    validationCurve(X2, y2, Xval, yval, Xtest,ytest,initial_nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels);

    close all;
    plot(lambda_vec, error_train, lambda_vec, error_val, lambda_vec, error_test);
    legend('Train', 'Cross Validation','Test');
    xlabel('lambda');
    ylabel('Error');
    fprintf('Program paused. Press enter to continue.\n');
    %print -deps lambda.eps;
    pause; %pause allows the plot to actually show up, if no pause it won't show properly
end

%% If not training then you can get the theta matrices from these
% text files
if train == 48
    
    load theta1.mat
    load theta2.mat
    load theta3.mat
    
end

% calculate the predictions here
[pred1 h3] = predict(Theta1, Theta2, Theta3, X2);

save -ascii pred1.txt pred1
save -ascii prob1.txt h3


if testt == 49
  fprintf('\nTraining Set Accuracy: %f\n', mean(double(pred1 == y2)) ...
					   * 100);
  trues = zeros(5,5);
    for i=1:length(pred1)
   
        trues(pred1(i),y2(i)) = trues(pred1(i),y2(i)) +1;
    
    end
    trues
end

%count out how may are wrong etc in what way in training mode
if train == 49
  fprintf('\nTraining Set Accuracy: %f\n', mean(double(pred1 == y2)) ...
					   * 100);
  trues = zeros(5,5);
    for i=1:length(pred1)
   
        trues(pred1(i),y2(i)) = trues(pred1(i),y2(i)) +1;
    
    end
    trues
  
    [pred h3] = predict(Theta1, Theta2, Theta3, Xtest);
    save -ascii predtest.txt pred
    fprintf('\nTesting Set Accuracy: %f\n', mean(double(pred == ytest)) ...
					    * 100);
    trues = zeros(5,5);
    for i=1:length(pred)
   
        trues(pred(i),ytest(i)) = trues(pred(i),ytest(i)) +1;
    
    end
    trues
    [pred2 h3] = predict(Theta1, Theta2, Theta3, Xval);
    save -ascii predcv.txt pred2
    fprintf('\nCV Set Accuracy: %f\n', mean(double(pred2 == yval)) ...
            * 100);
    trues = zeros(5,5);
    for i=1:length(pred2)
   
        trues(pred2(i),yval(i)) = trues(pred2(i),yval(i)) +1;
    
    end
    trues
end
    

fprintf('\nclassifier done!\n');
