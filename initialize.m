function [pi, pi_s, mu, alpha, alphas, phi, alphan] = initialize( theta, v, phi, scaling, gamcoeffs, Ms, betacoeffs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Initialization
N = length(theta);
smax = max(scaling(:,1));


% Block process for s and i

% draw for alpha(s)
count = 0;
for iter = 1:smax+1
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
pi_s(1) = 1; %betarnd(betacoeffs(7),betacoeffs(8));
% draw for pi(r)
temp1 = 0;
temp2 = 0;
count = 0;
for jiter = 1:Ms(1)
    if theta(count+jiter) == 0
        temp2 = temp2+1;
    else
        temp1 = temp1+1;
    end
end
count = count + Ms(1);
er = temp1 + betacoeffs(1);
fr = temp2 + betacoeffs(2);
pi_s(5) = betarnd(er,fr);
% draw for pi^0(s) and draw for pi^1(s)
count = 16;
temppi = zeros(smax+1,2);
for iter = 2:smax
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for jiter = 1:Ms(iter)
        if theta(count+jiter) == 0
            if theta(scaling(count+jiter,2)) == 0
                temp3 = temp3 +1;
            else
                temp4 = temp4+1;
            end
        else
            if theta(scaling(count+jiter,2)) == 0
                temp1 = temp1 +1;
            else
                temp2 = temp2+1;
            end
        end
    end
    e0 = betacoeffs(3) + temp1;
    f0 = betacoeffs(4) + temp2;
    e1 = betacoeffs(5) + temp3;
    f1 = betacoeffs(6) + temp4;
    temppi(iter,1) = betarnd(e0*Ms(iter),f0*Ms(iter));
    temppi(iter,2) = betarnd(e1*Ms(iter),f1*Ms(iter));
end
% assign value locations for the different pi into pi_s
for iter = 2:4
    pi_s(iter) = pi_s(1);
end
for iter = 6:16
    pi_s(iter) = pi_s(5);
end
stemp = 2;
count = 48;
start = 16;
for iter = 17:N
    if theta(scaling(iter,2)) == 0
        pi_s(iter) = temppi(stemp,1);
    else
        pi_s(iter) = temppi(stemp,2);
    end
    if iter == start+count
        start = start+count;
        count = count*4;
        stemp = stemp+1;
    end
end
% draw for alphan
alphan = gamrnd(gamcoeffs(1)+N/2,gamcoeffs(2)+(v-phi*theta).'*(v-phi*theta)/2);

% Update alpha(s,i)
for iter = 1:N
    alpha(iter) = alphas(scaling(iter,1)+1)+alphan*phi(:,iter).'*phi(:,iter); % +1 to account for s = 0
end
% Update mu(s,i)
sum = 0;
for jiter = 1:N
    sum = sum + phi(:,jiter)*theta(jiter);
end
for iter = 1:N
    tempsum = sum - phi(:,iter)*theta(iter);
    vhat = v-tempsum;
    mu(iter) = 1/alpha(iter)*alphan*phi(:,iter).'*vhat;
end
% Draw and calculate intermediate pi value = A
% calculate pi = A/(1+A) (verify) to update pi
for iter = 1:N
    temp1 = pi_s(iter)*1/sqrt(alphas(scaling(iter,1)+1))*randn(1);% +1 to account for s = 0
    temp2 = (1-pi_s(iter))*(mu(iter)+1/sqrt(alpha(iter))*randn(1));
    temp3 = temp1/temp2;
    pi(iter) = temp3/(1+temp3);
end

end



