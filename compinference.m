function [ theta, pi, mu, alpha, alphan, delta ] = compinference( theta, pi, mu, alpha, phi, alphan, v, scaling)
%compinference This function implements the MCMC inference model for
%zero-tree compressive sensing.  It only runs a single loop, so iterations
%must be run in the code outside the function
%   All these values are ainput and updated in the function:
%   theta - The estimate of the theta vector - Nx1
%   pi - The estimate of the pi values.  This is a Nx1 vector with
%   the values of pi for each theta value
%   as the vectors for each variable
%   mu - The Nx1 vector of mu values
%   alpha - The Nx1 vector alpha values
%   alphan - the noise variance
%   These values are input only for use in the update:
%   phi - The matrix of basis functions
%   v - The values from the compressive sampling of the image
%   scaling is a vector with the s and i values in the corresponding
%   locations.  It is Nx2 with s values in the first column and i values in
%   the second


% Initialization
N = length(theta);

% Block process for s and i

% Draw for theta(s,i)
pick = rand(N,1); % draws for which distribution to use
check =  pick < pi; % creates a vector of 1s and zeros
% use the 1s to use the non-zero distribution and the zeros to pick the zero distribution
theta = check.*(mu+diag(1./alpha)*randn(N,1));
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

