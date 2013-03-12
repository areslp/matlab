function [f,df]=hill_obj(x,dims,ii,dd,pars);
% 
% computes the objective function and gradient of the non-convex formulation of MVU. 
%
% copyright by Kilian Q. Weinberger, 2006
%
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    X=reshape(x,dims);
    [df,f]=computegr(X,ii(:,1),ii(:,2),dd);
    df=df(:).*8;

    if(pars.eta>0) % incorporate the trace term
      df=df-pars.eta.*x;
      f=f-pars.eta.*sum(sum(x.^2));
    end;

