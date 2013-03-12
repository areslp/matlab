% Version 1.000
%
% Code provided by Ruslan Salakhutdinov and Geoff Hinton
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from our
% web page.
% The programs and documents are distributed without any warranty, express or
% implied.  As the programs were written for research purposes only, they have
% not been tested to the degree that would be advisable in any important
% application.  All use of these programs is entirely at the user's own risk.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


function [f, df] = cg_update(VV, Dim, XX)

    l1 = Dim(1);
    l2 = Dim(2);
    l3 = Dim(3);
    l4 = Dim(4);
    l5 = Dim(5);
    l6 = Dim(6);
    l7 = Dim(7);
    l8 = Dim(8);
    l9 = Dim(9);
    N = size(XX, 1);

    % Extract weights back from VV
    w1 = reshape(VV(1:(l1 + 1) * l2), l1 + 1, l2);
    xxx = (l1 + 1) * l2;
    w2 = reshape(VV(xxx + 1:xxx + (l2 + 1) * l3), l2 + 1, l3);
    xxx = xxx + (l2 + 1) * l3;
    w3 = reshape(VV(xxx + 1:xxx + (l3 + 1) * l4), l3 + 1, l4);
    xxx = xxx + (l3 + 1) * l4;
    w4 = reshape(VV(xxx + 1:xxx + (l4 + 1) * l5), l4 + 1, l5);
    xxx = xxx + (l4 + 1) * l5;
    w5 = reshape(VV(xxx + 1:xxx + (l5 + 1) * l6), l5 + 1, l6);
    xxx = xxx + (l5 + 1) * l6;
    w6 = reshape(VV(xxx + 1:xxx + (l6 + 1) * l7), l6 + 1, l7);
    xxx = xxx + (l6 + 1) * l7;
    w7 = reshape(VV(xxx + 1:xxx + (l7 + 1) * l8), l7 + 1, l8);
    xxx = xxx + (l7 + 1) * l8;
    w8 = reshape(VV(xxx + 1:xxx + (l8 + 1) * l9), l8 + 1, l9);

    % Evaluate points 
    XX = [XX ones(N, 1)];
    w1probs = 1 ./ (1 + exp(-XX * w1));         w1probs = [w1probs ones(N,1)];
    w2probs = 1 ./ (1 + exp(-w1probs * w2));    w2probs = [w2probs ones(N,1)];
    w3probs = 1 ./ (1 + exp(-w2probs * w3));    w3probs = [w3probs ones(N,1)];
    w4probs = w3probs * w4;                     w4probs = [w4probs ones(N,1)];
    w5probs = 1 ./ (1 + exp(-w4probs * w5));    w5probs = [w5probs ones(N,1)];
    w6probs = 1 ./ (1 + exp(-w5probs * w6));    w6probs = [w6probs ones(N,1)];
    w7probs = 1 ./ (1 + exp(-w6probs * w7));    w7probs = [w7probs ones(N,1)];
    XXout   = 1 ./ (1 + exp(-w7probs * w8));

    % Compute gradients
    f = (-1 / N) * sum(sum(XX(:,1:end - 1) .* log(XXout) + (1 - XX(:,1:end - 1)) .* log(1 - XXout)));
    IO = (1 / N) * (XXout - XX(:,1:end - 1));
    Ix8 = IO;
    dw8 = w7probs'*Ix8;

    Ix7 = (Ix8 * w8') .* w7probs .* (1 - w7probs); 
    Ix7 = Ix7(:,1:end - 1);
    dw7 = w6probs' * Ix7;

    Ix6 = (Ix7*w7') .* w6probs .* (1 - w6probs); 
    Ix6 = Ix6(:,1:end - 1);
    dw6 = w5probs' * Ix6;

    Ix5 = (Ix6 * w6') .* w5probs .* (1 - w5probs); 
    Ix5 = Ix5(:,1:end - 1);
    dw5 = w4probs' * Ix5;

    Ix4 = (Ix5 * w5');
    Ix4 = Ix4(:,1:end - 1);
    dw4 = w3probs' * Ix4;

    Ix3 = (Ix4 * w4') .* w3probs .* (1 - w3probs); 
    Ix3 = Ix3(:,1:end - 1);
    dw3 = w2probs' * Ix3;

    Ix2 = (Ix3 * w3') .* w2probs .* (1 - w2probs); 
    Ix2 = Ix2(:,1:end - 1);
    dw2 = w1probs' * Ix2;

    Ix1 = (Ix2 * w2') .* w1probs .* (1 - w1probs); 
    Ix1 = Ix1(:,1:end - 1);
    dw1 = XX' * Ix1;

    % Return gradients
    df = [dw1(:)' dw2(:)' dw3(:)' dw4(:)' dw5(:)' dw6(:)' dw7(:)' dw8(:)']';
