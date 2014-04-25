function [pi, pi_s, mu, alpha, alphas, phi, alphan] = initialize( theta, v, phi, scaling, gamcoeffs, Ms, betacoeffs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Initialization
N = length(theta);
M = length(v);
smax = max(scaling(:,1));

pi = zeros(N,1);
pi_s = zeros(N,1);
mu = zeros(N,1);
alpha = zeros(N,1);
alphas = zeros(smax+1,1);

% Block process for s and i

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
end



