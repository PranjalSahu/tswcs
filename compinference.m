function [ theta, pi, pi_s,  mu, alpha, alphas, alphan ] = compinference( theta, pi, pi_s, mu, alpha, alphas, phi, alphan, v, scaling, gamcoeffs, Ms, betacoeffs)
%compinference This function implements the MCMC inference model for
%zero-tree compressive sensing.  It only runs a single loop, so iterations
%must be run in the code outside the function
%   All these values are ainput and updated in the function:
%   theta - The estimate of the theta vector - Nx1
%   pi - The estimate of the pi values.  This is a Nx1 vector with
%   the values of pi for each theta value as the vectors for each variable.
%   (pi_tilda)
%   pi_s - The values of pi for different values of s and different
%   parents.  (Nx1)
%   mu - The Nx1 vector of mu values
%   alpha - The Nx1 vector alpha values (alpha_tilda)
%   alphas - a vector of the values of alpha for different s-levels
%   alphan - the noise variance
%   These values are input only for use in the update:
%   phi - The matrix of basis functions
%   v - The values from the compressive sampling of the image
%   scaling is a vector with the s and i values in the corresponding
%   locations.  It is Nx2 with s values in the first column and the parent
%   locations in the 2nd column
%   gamcoeffs - 4x1 vector with [a0, b0, c0, d0] in it for the gamma
%   distributions
%   Ms is the passed in sizes of the different scales in order.
%   betacoeffs - 8x1 vector with [e0r,f0r, e0s0,f0s0,e0s1,f0s1,e0sc,explicit_coeffs] in it for
%   the beta distributions


% Initialization
N = length(theta);
M = length(v);
smax = max(scaling(:,1));

% Block process for s and i

% Draw for theta(s,i)
pick = rand(N,1); % draws for which distribution to use
check =  pick < pi; % creates a vector of 1s and zeros
% use the 1s to use the non-zero distribution and the zeros to pick the zero distribution
theta = check.*(mu+randn(N,1)./sqrt(alpha));
theta(1:betacoeffs(8)) = v(1:betacoeffs(8));

% Update alpha(s,i)
alpha = alphas(scaling(:,1)+1) + alphan*sum(phi.*phi,1).';

% Update mu(s,i)
phi_times_theta = bsxfun(@times, phi, theta.');
phi_times_theta_j = bsxfun(@minus, sum(phi_times_theta,2), phi_times_theta);
vhat = bsxfun(@minus, v, phi_times_theta_j);
mu = alphan./alpha.*sum(phi.*vhat,1).';

% Draw and calculate intermediate pi value = A
% calculate pi = A/(1+A) (verify) to update pi
temp1 = pi_s./sqrt(alphas(scaling(:,1)+1)).*randn(N,1); % +1 to account for s = 0
temp2 = (1-pi_s).*(mu+1./sqrt(alpha).*randn(N,1));
temp3 = temp1./temp2;
pi = temp3./(1+temp3);
pi(1:4) = 1;

temp3 = 0;
cin = 2+gamcoeffs(3);
for iter = 1:4
    temp3 = temp3+theta(iter)^2;
    din = temp3*.5+gamcoeffs(4);
    alphas(iter) = gamrnd(cin,din);
end

% draw for alpha(s)
count = 4;
for iter = 2:smax+1
    temp1 = 0;
    temp3 = 0;
    for jiter = 1:Ms(iter)
        if theta(count+jiter) == 0
            temp2 = 0;
        else
            temp2 = 1;
        end
        temp1 = temp1+temp2;
        temp3 = temp3+theta(count+jiter)^2;
    end
    count = count + Ms(iter);
    cin = temp1*.5+gamcoeffs(3);
    din = temp3*.5+gamcoeffs(4);
    alphas(iter) = gamrnd(cin,din);
end
% draw for pi(sc)
% pi_s(1) = 1; %betarnd(betacoeffs(7),betacoeffs(8));
% draw for pi(r)
temp1 = 0;
temp2 = 0;
count = Ms(1);

temp1 = nnz(theta(count+(1:Ms(2)+1)));
temp2 = Ms(2) - temp1;

er = temp1 + betacoeffs(1);
fr = temp2 + betacoeffs(2);
pi_s(5) = betarnd(er,fr);
% draw for pi^0(s) and draw for pi^1(s)
count = 16;
temppi = zeros(smax-1,2);
for iter = 1:smax-1
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for jiter = 1:Ms(iter+2)
        if theta(count+jiter) == 0
            if theta(scaling(count+jiter,2)) == 0
                temp2 = temp2 +1;
            else
                temp4 = temp4+1;
            end
        else
            if theta(scaling(count+jiter,2)) == 0
                temp1 = temp1 +1;
            else
                temp3 = temp3+1;
            end
        end
    end
    e0 = betacoeffs(3) + temp1;
    f0 = betacoeffs(4) + temp2;
    e1 = betacoeffs(5) + temp3;
    f1 = betacoeffs(6) + temp4;
    temppi(iter,1) = betarnd(e0*Ms(iter+2),f0*Ms(iter+2));
    temppi(iter,2) = betarnd(e1*Ms(iter+2),f1*Ms(iter+2));
end
% assign value locations for the different pi into pi_s
% pi_s(2:4) = pi_s(1);
pi_s(6:16) = pi_s(5);

stemp = 2;
stemp2 = stemp;
count = Ms(stemp+1);
start = sum(Ms(1:stemp));

for iter = Ms(stemp2+1):N
    if theta(scaling(iter,2)) == 0
        pi_s(iter) = temppi(stemp-1,1);
    else
        pi_s(iter) = temppi(stemp-1,2);
    end
    if iter == start+count
        start = start+count;
        count = count*4;
        stemp = stemp+1;
    end
end

% draw for alphan
v_minus_phi_theta = v - phi*theta;
alphan = gamrnd(gamcoeffs(1)+M/2,gamcoeffs(2)+v_minus_phi_theta.'*v_minus_phi_theta/2);
end

