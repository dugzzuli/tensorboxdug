clear all; close all;

addpath Function_SPC
addpath plotting_function

%%% make synthetic data
N = 50;

lin = linspace (0, 2, N);
[x, y, z] = meshgrid (lin, lin, lin);

sig = 0.3;
c1 = exp( -1/sig* ((x-1.3).^2 + (y-.3).^2 + (z-.3).^2) ) ;
c2 = exp( -1/sig* ((x-.3).^2 + (y-1.3).^2 + (z-.3).^2) ) ;
c3 = exp( -1/sig* ((x-.3).^2 + (y-.3).^2 + (z-1.3).^2) ) ;
c4 = exp( -1/sig* ((x-1.4).^2 + (y-1.4).^2 + (z-1.4).^2) ) ;
sig2 = 10.0;
X0 = 2*(c1 + c2 + c3 + c4) + 0.5*exp( - 1/sig2* ((x-.3).^2 + (y-.3).^2 + (z-.3).^2));

%%% make missing entry
II = size(X0);
N  = prod(II);
missing_rate = 0.8;

idd = (randperm(N) > N*0.8);
Q   = reshape(idd,II);
T   = zeros(II);
T(Q)= X0(Q);

%% hyperparameters and run SPC-TV

TVQV    = 'tv';        % 'tv' or 'qv' ;
rho     = [0.01 0.01 0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
K       = 10;          % Number of components which are updated in one iteration.
SNR     = 50;          % error bound
nu      = 0.01;        % threshold for R <-- R + 1.
maxiter = 10000;       % maximum number of iteration
tol     = 1e-7;        % tolerance
out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.

[Xtv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);


%% hyperparameters and run SPC-QV

TVQV    = 'qv';        % 'tv' or 'qv' ;
rho     = [1.0 1.0 1.0]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
K       = 10;          % Number of components which are updated in one iteration.
SNR     = 50;          % error bound
nu      = 0.01;        % threshold for R <-- R + 1.
maxiter = 10000;       % maximum number of iteration
tol     = 1e-7;        % tolerance
out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.

[Xqv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);


