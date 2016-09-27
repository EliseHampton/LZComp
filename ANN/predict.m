function [p h3] = predict(Theta1, Theta2, Theta3, X)
% PREDICT Predict the label of an input given a trained neural network
% p = PREDICT(Theta1, Theta2, Theta3, X) outputs the predicted label of X given the
% trained weights of a neural network (Theta1, Theta2, Theta3)

% number of examples and number of possible layers
m = size(X, 1);
num_labels = size(Theta3, 1);

% return the predictions and the probabilities
p = zeros(size(X, 1), 1);

h1 = sigmoid([ones(m, 1) X] * Theta1');
h2 = sigmoid([ones(m, 1) h1] * Theta2');
h3 = sigmoid([ones(m, 1) h2] * Theta3');

% dummy is actually the probabilities for each possible
% classification for each example which may be useful
[dummy, p] = max(h3, [], 2);

end
