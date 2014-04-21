function [ theta, pi, pi_s,  mu, alpha, alphas, alphan, delta ] = compinference( theta, pi, pi_s, mu, alpha, alphas, phi, alphan, v, scaling)
%compinference This function implements the MCMC inference model for
%zero-tree compressive sensing.  It only runs a single loop, so iterations
%must be run in the code outside the function
%   All these values are ainput and updated in the function:
%   theta - The estimate of the theta vector - Nx1
%   pi - The estimate of the pi values.  This is a Nx1 vector with
%   the values of pi for each theta value as the vectors for each variable.
%   pi_s - The values of pi for different values of s and different parents.  This is Nx1.
%   mu - The Nx1 vector of mu values
%   alpha - The Nx1 vector alpha values
%   alphas - a vector of the values of alpha for different s-levels
%   alphan - the noise variance
%   These values are input only for use in the update:
%   phi - The matrix of basis functions
%   v - The values from the compressive sampling of the image
%   scaling is a vector with the s and i values in the corresponding
%   locations.  It is Nx2 with s values in the first column and i values in
%   the second


% Initialization
N = length(theta);
smax = max(scaling(:,1));


% Block process for s and i

% Draw for theta(s,i)
pick = rand(N,1); % draws for which distribution to use
check =  pick < pi; % creates a vector of 1s and zeros
% use the 1s to use the non-zero distribution and the zeros to pick the zero distribution
theta = check.*(mu+diag(1./alpha)*randn(N,1));
% Update alpha(s,i)
for iter = 1:N
    alpha(iter) = alpha(iter)+alphan*phi(:,iter).'*phi(:,iter);
end
% Update mu(s,i)
for iter = 1:N
    sum = 0;
    for jiter = 1:N
        if jiter ~=iter
            sum = sum + phi(:,jiter)*theta(jiter);
        end
    end
    vhat = v-sum;
    mu = 1/alpha(iter)*alphan*phi(:,iter).'*vhat;
end
% Draw and calculate intermediate pi value = A
% calculate pi = A/(1+A) (verify) to update pi
for iter = 1:N
    temp1 = pi_s(iter)*1/sqrt(alpha(scaling(iter,1)))*randn(1);
    temp2 = (1-pi_s(iter))*(mu(iter)+1/sqrt(alpha(scaling(iter,1)))*randn(1));
    temp3 = temp1/temp2;
    pi(iter) = temp3/(1+temp3);
end
% draw for alpha(s)

% draw for pi(r)

% draw for pi^0(s)

% draw for pi^1(s)

% assign value locations for the different pi into pi_s

% draw for alphan



end

