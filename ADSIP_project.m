% ADSIP Project Script
% This is the top level script to run compressive sensing with wavelet
% transforms, zero-tree structure exploitation, and Bayesian estimation
% using Gibbs sampling
clear all;

% Initialize values
smax = 6; % set for scaling calculation below; if not used elsewhere, move there
M0 = 4; % set for scaling calculation below; if not used elsewhere, move there
Msvec = [4 3*4.^(1:smax)];
N = sum(Msvec);
gamcoeffs = 10^-6*[1 1 1 1];
explicit_coeffs = 4;
betacoeffs = [.9*Msvec(2) .1* Msvec(2) 1/N 1-1/N .5 .5 1 explicit_coeffs];
%number of elements
% Read in picture
% Picture Size
% Create Psi Matrix (Wavelet Transform)
% Create Permutation Matrices
% Create Phi Matrix (Sampling Matrix)
phi = sampling_matrix(2000,sum(Msvec),explicit_coeffs);
% pperm = randperm(sum(Msvec),2000);
% Phitemp = zeros(N);
% jiter = 1:N;
% Phitemp(1,:) = 1/sqrt(N)*cos(3.14159*(1-1)*(2*jiter-1)/(2*N)).';
% iter = 2:N;
% Phitemp(2:end,:) = sqrt(2)/sqrt(N)*cos(3.14159*bsxfun(@times,(iter-1).',(2*jiter-1))/(2*N));
% phi = Phitemp(pperm,:);

% Add option here to directly sample "DC" coefficients?
% QUESTION: Do we want to be able to run both simultaneously?
% Combine Wavelet and Sampling matrices
% Create scaling (parent and s-level vectors
% NOTE: MAY NEED TO MOVE SOME OF THE VALUES HERE UP TOP
N = M0*4^smax;
scaling = zeros(N,2);
s = 0;
Mtot = M0;
Ms = M0*3/4;
stepa = 0;
for iter = 1:(N/4) % No children for last level
    scaling(iter,1) = s;
    if s == 0
        scaling(iter,2) = 0;
    elseif s == 1
        scaling(iter,2) = 0;
        num = stepa+2*(row-1)+skip*(col-1);
        scaling(num,2) = iter;
        scaling(num+1,2) = iter;
        scaling(num+skip/2,2) = iter;
        scaling(num+skip/2+1,2) = iter;
        if row ==2^s
            row = 1;
            col = col+1;
        else
            row = row+1;
        end
    else
        num = stepa+2*(row-1)+skip*(col-1);
        scaling(num,2) = iter;
        scaling(num+1,2) = iter;
        scaling(num+skip/2,2) = iter;
        scaling(num+skip/2+1,2) = iter;
        if row ==2^s
            row = 1;
            col = col+1;
        else
            row = row+1;
        end
    end
    if iter == Mtot
        s = s+1;
        Ms = Ms*4;
        stepa = iter+Ms+1;
        Mtot = Mtot*4;
        row = 1;
        col = 1;
        skip = 2^(s+2);
    elseif iter == Mtot-2*Ms/3
        stepa = stepa+Mtot;
        row = 1;
        col = 1;
    elseif iter == Mtot-Ms/3
        stepa = stepa+Mtot;
        row = 1;
        col = 1;
    end
end
scaling(N/4+1:end,1) = 6;
% Other Initialization?

% Main algorithm
% Sample Image
load test_images
[V, Psi, P] = multilevel_haar(test_image{4},1);
theta_true = P*V(:);
figure
plot(theta_true)
v = phi*theta_true;%+5*randn(2000,1);

% Initialize values for bayesian model (theta, pi, alpha, etc, se inputs to compinference
theta = ([200*ones(60,1);20*ones(N-60,1)].*randn(N,1)).*(rand(N,1) < 0.2);
%theta = theta_true;
theta(1:explicit_coeffs) = v(1:explicit_coeffs);
figure
plot(theta)
[theta, pi, pi_s, mu, alpha, alphas, phi, alphan] = initialize(theta, v, phi, scaling, ...
                                                  gamcoeffs, Msvec, ...
                                                  betacoeffs);

% Loop for Bayesian model (use compinference call here)


L = 200;
mse = zeros(L,1);
h_waitbar = waitbar(0,'Bayesian Inference...');

for l = 1:L
    waitbar(l/L,h_waitbar,'Bayesian Inference...');

    old_theta = theta;
    [theta, pi, pi_s,  mu, alpha, alphas, alphan] = compinference(theta, pi, pi_s, ...
                                                      mu, alpha, alphas, phi, alphan, ...
                                                      v, scaling, gamcoeffs, Msvec, ...
                                                      betacoeffs);
    mse(l) = norm(old_theta-theta);
end

close(h_waitbar);
h_waitbar = waitbar(0,'Sampling Theta...');
samplestot = 100;
thetatot = 0;
for l = 1:samplestot
    waitbar(l/samplestot,h_waitbar,'Sampling Theta...');

    % Draw for theta(s,i)
    pick = rand(N,1); % draws for which distribution to use
    check =  pick < pi; % creates a vector of 1s and zeros
    % use the 1s to use the non-zero distribution and the zeros to pick the zero distribution
    theta = check.*(mu+randn(N,1)./sqrt(alpha));
    theta(1:betacoeffs(8)) = v(1:betacoeffs(8));
    
    %[theta, pi, pi_s,  mu, alpha, alphas, alphan] = compinference(theta, pi, pi_s, ...
    %                                                  mu, alpha, alphas, phi, alphan, ...
    %                                                  v, scaling, gamcoeffs, Msvec, ...
    %                                                  betacoeffs);
    thetatot = thetatot+theta;
end
close(h_waitbar);
theta = thetatot/samplestot;
finalerror = norm(theta_true-theta)
figure(1); clf;
subplot(2,2,1);
showme(reshape(P.'*significant(theta,0.05),128,[]));
subplot(2,2,2);
showme(reshape(transform(Psi,P.'*theta),128,[]));
subplot(2,2,3);
showme(reshape(P.'*theta,128,[]));
subplot(2,2,4);
plot(mse);

% Inverse transform converged estimate
% Reorder with permutations (may be combined with above step)
% Generate Picture comparison (original vs our reconstructed)
% Calculate actual metric





