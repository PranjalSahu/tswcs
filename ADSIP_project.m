% ADSIP Project Script
% This is the top level script to run compressive sensing with wavelet
% transforms, zero-tree structure exploitation, and Bayesian estimation
% using Gibbs sampling

% Initialize values
%gamcoeffs
%betacoefs
smax = 6; % set for scaling calculation below; if not used elsewhere, move there
M0 = 4; % set for scaling calculation below; if not used elsewhere, move there
%number of elements
% Read in picture
% Picture Size
% Create Psi Matrix (Wavelet Transform)
% Create Permutation Matrices
% Create Phi Matrix (Sampling Matrix)
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
for iter = 1:(N/4) % No children for last level
    scaling(iter,1) = s;
    if s == 0
        scaling(iter,2) = 0;
    elseif s == 1
        scaling(iter,2) = 0;
        num = step+2*(row-1)+skip*(col-1);
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
        num = step+2*(row-1)+skip*(col-1);
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
        step = iter+Ms+1;
        Mtot = Mtot*4;
        row = 1;
        col = 1;
        skip = 2^(s+2);
    elseif iter == Mtot-2*Ms/3
        step = step+Mtot;
        row = 1;
        col = 1;
    elseif iter == Mtot-Ms/3
        step = step+Mtot;
        row = 1;
        col = 1;
    end
end
% Other Initialization?

% Main algorithm
% Initialize values for bayesian model (theta, pi, alpha, etc, se inputs to compinference
% Sample Image
% Loop for Bayesian model (use compinference call here)
% Inverse transform converged estimate
% Reorder with permutations (may be combined with above step)
% Generate Picture comparison (original vs our reconstructed)
% Calculate actual metric





