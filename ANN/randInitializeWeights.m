function W = randInitializeWeights(L_in, L_out)
%RANDINITIALIZEWEIGHTS Randomly initialises the weights of a layer
%with L_in (previous layer size)
%incoming connections and L_out (size of this layer) outgoing connections
 
W = zeros(L_out, 1 + L_in);

eps_init = 0.12;

W = (rand(L_out,1+L_in)*2*eps_init)-eps_init;

end
