function stats = runSimulatedVar(opts)
%STATS=RUNSIMULATEDVAR(OPTS)
% run an analysis of GCA on simulated VAR data with parameters specified by
% options structure opts
%
% updated on 03.06.22 to add call to runGcaTrAlt + measuring G between all
% sensor pairs, order dependent

if nargin<1
    opts=setDefaultOptions;
end
    
K=opts.K; % number of latent sources
L=opts.L; % max lag in AR model
D=opts.D; % number of sensors
N=opts.N; % number of time points
maxIter=opts.maxIter;
sigma_inn=opts.sigma_inn; % standard deviation of innovation process

P=opts.P; % # of pairs to search for

% generate the VAR process
B=opts.B;

S(:,1:L)=randn(K,L);
for n=L+1:N
    for p=1:L
        S(:,n)=B(:,:,p)*S(:,n-p);
    end
    S(:,n)=S(:,n)+sigma_inn*randn(K,1); 
end

% calculate the true G-causality between all pairs of latent sources
lPairs=combnk(1:K,2);
lPairs=cat(1,lPairs,lPairs(:,2:-1:1)); % order-dependent
nLatentPairs=size(lPairs,1);
G_src=nan(nLatentPairs,1);
G_src_2D=nan(K,K);
for p=1:nLatentPairs
    G_src(p,1)=vanillaG(S(lPairs(p,1),:)',S(lPairs(p,2),:)',L);
    G_src_2D(lPairs(p,1),lPairs(p,2))=G_src(p,1);
end


% generate observations
A=rand(D,K);
X=(A*S).';

% mean center
X = bsxfun(@minus, X, mean(X,1));

% run GCA
[What,Vhat,gcs,gcaStats] = runGca(X,L,P,maxIter,inf);

% normalize solutions to unit norm
What= What./ repmat( sqrt(sum(What.^2)) , D , 1);
Vhat= Vhat./ repmat( sqrt(sum(Vhat.^2)) , D , 1);

%% get the the ideal solution
pinvAt=pinv(A');
Wtrue=pinvAt(:,1);
Vtrue=pinvAt(:,2);

%% recover the time sources of driving and driven signals
yhat=X*What;
zhat=X*Vhat;
cpairs=combnk(1:P,2);
cpairs=cat(1,cpairs,cpairs(:,2:-1:1));
nPairs=size(cpairs,1);
G_gca=nan(nPairs,1);
G_gca_2D=nan(P,P);
for p=1:P
    G_gca(p,1)=vanillaG(yhat(:,p),zhat(:,p),L);
    for q=1:P
        G_gca_2D(p,q)=vanillaG(yhat(:,p),zhat(:,q),L);
    end
end

% G-causality between all observed signals
pairs=combnk(1:D,2);
pairs=cat(1,pairs,pairs(:,2:-1:1)); % order-dependent
nPairs=size(pairs,1);
G_sensor=nan(nPairs,1);
G_sensor_2D=nan(D,D);
for p=1:nPairs
    G_sensor(p,1)=vanillaG(X(:,pairs(p,1)),X(:,pairs(p,2)),L);
    G_sensor_2D(pairs(p,1),pairs(p,2))=G_sensor(p,1);
end

% collect outputs
stats.G_sensor=G_sensor;
stats.G_sensor_2D=G_sensor_2D;

stats.G_gca=G_gca;
stats.G_gca_2D=G_gca_2D;

stats.G_src=G_src;
stats.G_src_2D=G_src_2D;

stats.yhat=yhat;
stats.zhat=zhat;

stats.Wtrue=Wtrue;
stats.Vtrue=Vtrue;
stats.What=What;
stats.Vhat=Vhat;

stats.gcs=gcs; 

stats.X=X;
stats.A=A;
stats.S=S;

stats.opts=opts;

stats.gcaStats=gcaStats;

end

