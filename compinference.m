function [ theta, pi, mu, alpha, alphan, delta ] = compinference( theta, pi, mu, alpha, phi, alphan, v, scaling)
%compinference This function implements the MCMC inference model for
%zero-tree compressive sensing.  It only runs a single loop, so iterations
%must be run in the code outside the function
%   All these values are ainput and updated in the function:
%   theta - The estimate of the theta vector
%   pi - The estimate of the pi values.  This is a 3x1 vector with
%   [pi_r, pi^0_s, pi^1_s] as the elements
%   mu - The mu value
%   alpha - The alpha value
%   alphan - the noise variance
%   These values are input only for use in the update:
%   phi - The matrix of basis functions
%   v - The values from the compressive sampling of the image
%   scaling is the scaling information for breaking up theta (REMOVE?)


% Initialization


% Loop over theta (s and i) - or do block

% Draw for theta(s,i)

% Update alpha(s,i)

% Update mu(s,i)

% Draw and calculate intermediate pi value = A

% calculate pi = A/(1+A) (verify) to update pi

% draw for alpha(s)

% draw for pi(r)

% draw for pi^0(s)

% draw for pi^1(s)

% draw for alphan



end

