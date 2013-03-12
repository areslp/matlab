function [machine, err] = train_grbm(X, h, eta, max_iter, weight_cost)
%TRAIN_GRBM Trains a Gaussian Restricted Boltzmann Machine using CD1
%
%   [machine, err] = train_grbm(X, h, eta, max_iter, weight_cost)
%
% Trains a first-order Restricted Boltzmann Machine on dataset X. The RBM
% has h hidden nodes (default = 30). The training is performed by means of
% the contrastive divergence algorithm. The activation functions that
% is applied in the visible and hidden layers are binary stochastic.
% In the training of the RBM, the learning rate is determined by eta 
% (default = 0.1). The maximum number of iterations can be specified 
% through max_iter (default = 10). The variable weight_cost sets the amount 
% of weight decay that is employed (default = 0.0002).
% The trained RBM is returned in the machine struct.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Initialize some variables
    eta = 0.001;
    cd_iter = 1;
    max_iter = 200;
    weight_cost = 0.002;
    initial_momentum = 0.5;
    final_momentum   = 0.9;
    [n, v] = size(X);
    batch_size = 100;
    err = zeros(max_iter, 1);
    machine.type = 'rbm_continuous';
    machine.W    = randn(v, h) * .0001;
    machine.bias_upW   = zeros(1, h);
    machine.bias_downW = zeros(1, v);
    delta_E      = zeros(v, h);
    delta_E_bias = zeros(1, h);
    delta_X_bias = zeros(1, v);
    
    % Main loop
    for iter=1:max_iter
        
        % Set momentum and update batches
        ind = randperm(n);
        if iter <= 5
            momentum = initial_momentum;
        else
            momentum = final_momentum;
        end
        
        % Run for all mini-batches
        for batch=1:batch_size:n          
            if batch + batch_size <= n
            
                % Set values of visible nodes
                vis1 = full(X(ind(batch:min([batch + batch_size - 1 n])),:));
                
                % Compute probabilities for hidden nodes
                hid1 = 1 ./ (1 + exp(-bsxfun(@plus, vis1 * machine.W, machine.bias_upW)));
                
                % Sample states for hidden nodes
                hid_states = hid1 > rand(size(hid1));
                
                % Run contrastive divergence iterations
                for i=1:cd_iter
                    
                    % Compute probabilities for visible nodes
                    vis2 = bsxfun(@plus, hid_states * machine.W', machine.bias_downW);

                    % Compute probabilities for hidden nodes
                    hid2 = 1 ./ (1 + exp(-bsxfun(@plus, vis2 * machine.W, machine.bias_upW)));
                    
                    % Sample states for hidden nodes
                    hid_states = hid2 > rand(size(hid2));                
                end
                
                % Compute the weight updates
                delta_E = momentum * delta_E + (eta ./ batch_size) * ((vis1' * hid1 - vis2' * hid2) - (weight_cost * machine.W));
                delta_E_bias = momentum * delta_E_bias + (eta ./ batch_size) * (sum(hid1, 1) - sum(hid2, 1));
                delta_X_bias = momentum * delta_X_bias + (eta ./ batch_size) * (sum(bsxfun(@minus, vis2, machine.bias_downW), 1) - ...
                                                                                sum(bsxfun(@minus, vis1, machine.bias_downW), 1)); 
                              
                % Update the network weights
                machine.W = machine.W + delta_E;
                machine.bias_upW = machine.bias_upW + delta_E_bias;
                machine.bias_downW = machine.bias_downW + delta_X_bias;
                
                % Sum reconstruction error
                err(iter) = err(iter) + sum((vis1(:) - vis2(:)) .^ 2);
            end
        end 
        
        % Print reconstruction error estimate
        disp(['Iteration ' num2str(iter) ' of ' num2str(max_iter) ' (rec. error = ' num2str(err(iter) ./ n) ')...']);
%         subplot(1, 2, 1), imagesc(reshape(vis1(1,:), [28 28])'); colormap gray, colorbar
%         subplot(1, 2, 2), imagesc(reshape(vis2(1,:), [28 28])'); colormap gray, colorbar
%         drawnow
        
        % Jump to higher CD-n every now and then
        if iter == round(max_iter * .75)
            disp('Switching to CD-3...');
            cd_iter = 3;
            eta = eta / 3;
        end
    end
    