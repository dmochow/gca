function opts = setDefaultOptions
%OPTIONS = SETDEFAULTOPTIONS
% prepare an options structure to feed into runSimulatedVar with default
% fields for all required options

opts.K=3; % number of latent sources
opts.L=3; % max lag in AR model
opts.D=4; % number of sensors
opts.N=5000; % number of time points
opts.sigma_inn=1; % standard deviation of innovation process
opts.P=2; % # of pairs to search for
opts.maxIter=50; % max iterations in grouped coordinate descent

% generate the VAR process
opts.r1=0.9; opts.theta1=40/120*2*pi;
opts.r2=0.7; opts.theta2=10/120*2*pi;
opts.r3=0.8; opts.theta3=50/120*2*pi;

% create the VAR system matrix here
B(:,:,1)=[2*opts.r1*cos(opts.theta1) 0 0; -0.356 2*opts.r2*cos(opts.theta2) 0; 0 -0.3098 2*opts.r3*cos(opts.theta3) ]; % lag 1
B(:,:,2)=[-opts.r1.^2 0 0; 0.7136 -opts.r2.^2 0; 0 0.5 -opts.r3.^2];  % lag 2
B(:,:,3)=[0 0 0; -0.356 0 0; 0 -0.3098 0];  % lag 3

opts.B=B;
end

