% ADSIP Project Script
% This is the top level script to run compressive sensing with wavelet
% transforms, zero-tree structure exploitation, and Bayesian estimation
% using Gibbs sampling

% Initialize values
gamcoeffs = 10^-6*[1 1 1 1];
betacoeffs = [.9*12 .1*12 1/sum(Msvec) 1-1/sum(Msvec) .5 .5 1 1];
smax = 6; % set for scaling calculation below; if not used elsewhere, move there
M0 = 4; % set for scaling calculation below; if not used elsewhere, move there
Msvec = [4 3*4.^(1:smax)];
%number of elements
% Read in picture
% Picture Size
% Create Psi Matrix (Wavelet Transform)
% Create Permutation Matrices
% Create Phi Matrix (Sampling Matrix)
phi = sampling_matrix(6000,sum(Msvec));
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
% Other Initialization?

% Main algorithm
% Sample Image
v = 2*randn(6000,1);
% Initialize values for bayesian model (theta, pi, alpha, etc, se inputs to compinference
theta = randn(sum(Msvec),1);
[pi, pi_s, mu, alpha, alphas, phi, alphan] = initialize(theta, v, phi, scaling, gamcoeffs, Msvec, betacoeffs);
% Loop for Bayesian model (use compinference call here)
[theta, pi, pi_s,  mu, alpha, alphas, alphan] = compinference(theta, pi, pi_s, mu, alpha, alphas, phi, alphan, v, scaling, gamcoeffs, Msvec, betacoeffs);
% Inverse transform converged estimate
% Reorder with permutations (may be combined with above step)
% Generate Picture comparison (original vs our reconstructed)
% Calculate actual metric





