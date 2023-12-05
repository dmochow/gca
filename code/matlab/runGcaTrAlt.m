function [What,Vhat,G_gca,stats] = runGcaTrAlt(X,L,P,maxIters,maxcond)
%[WHAT,VHAT,G_GCA]=RUNGCATRALT(X,L,P,MAXITERS,MAXCOND)
%   fit Granger Components Analysis model on data X
%   
%   X: 2D array of observations, samples x dimension
%   L: scalar maximum lag parameter in GCA model (defaults to 1)
%   P: scalar parameter controlling the number of pairs to recover
%   (defaults to 3)
%   maxIters: scalar maximum number of iterations in coordinate descent,
%   defaults to 10
%   maxcond: scalar controlling maximum condition number of block
%   covariance matrices of the observations, defaults to 1000 (this is a
%   form of regularization)

optTol=1e-6; % internal stopping parameter

if nargin<5 || isempty(maxcond), maxcond=1000; end %  cov matrices with cond #s > maxcond will be regularized 
if nargin<4 || isempty(maxIters), maxIters=10; end
if nargin<3 || isempty(P), P=3; end
if nargin<2 || isempty(L), L=1; end
if nargin<1 || isempty(X), error('runGcaTrAlt needs at least one argument'); end

[N,D,n_trials]=size(X);
Xog=X; % save for G_gca calculation
if ndims(X)==2
    Xr=X(end:-1:1,:);
elseif ndims(X)==3
    Xr=X(end:-1:1,:,:);
else
    error('X must have 2 or 3 dims');
end

vt=0.01*randn(D,1);
wt=0.01*randn(D,1);
options = optimoptions('fmincon','SpecifyObjectiveGradient',false,'CheckGradients',false,...
    'MaxFunctionEvaluations',10000,'MaxIterations',4000);

Vhat=nan(D,P); % weights of the driven signal component
What=nan(D,P); % weights on the driving signal component

for p=1:P
    
    
    if ndims(X)==2
        [blkR1,Rtilde] = gcBlockCovariance(X,L);
        [blkR1r,Rtilder] = gcBlockCovariance(Xr,L);
    elseif ndims(X)==3
       
        blkR1=zeros(D*L,D*L);
        Rtilde=zeros(D*L,D*L);
        blkR1r=zeros(D*L,D*L);
        Rtilder=zeros(D*L,D*L);
        for tr=1:n_trials           
            tX = X(:,:,tr);
            [tblkR1,tRtilde] = gcBlockCovariance(tX,L);
            blkR1 = blkR1 + tblkR1/n_trials;
            Rtilde = Rtilde + tRtilde/n_trials;
            
            tXr = Xr(:,:,tr);
            [tblkR1r,tRtilder] = gcBlockCovariance(tXr,L);
            blkR1r=blkR1r+tblkR1r/n_trials;
            Rtilder=Rtilder+tRtilder/n_trials;
        end
    else 
        error('X can only be 2 or 3 dimensional');
    end
    
        
    % FIXME blkR1 is not positive definite!
    if cond(blkR1)>maxcond
        blkR1=ridgeCov(blkR1,maxcond);
    end
    
    if cond(blkR1r)>maxcond
        blkR1r=ridgeCov(blkR1r,maxcond);
    end
    
    if cond(Rtilde)>maxcond
        Rtilde=ridgeCov(Rtilde,maxcond);
    end
    
    if cond(Rtilder)>maxcond
        Rtilder=ridgeCov(Rtilder,maxcond);
    end
    
    fvals_v=nan(maxIters,1); % initialize for this pair
    fvals_w=nan(maxIters,1); % initialize for this pair
    fvals_t=nan(maxIters,1); % the 2D function G(v,w)
    
    for i=1:maxIters
        
        % first optimize for v
        funv = @(v) gctr_v(v,wt,blkR1,Rtilde,blkR1r,Rtilder,L) ; 
        [xf,fvalv] = fmincon(funv,0.01*randn(D,1),[],[],[],[],[],[],@(x)unitNormCon(x),options);
        vt=xf;
        
        % now optimize for w, and fix v to be the optimal value from above
        funw = @(w) gctr_w(w,vt,blkR1,Rtilde,blkR1r,Rtilder,L); 
        [xf,fvalw] = fmincon(funw,0.01*randn(D,1),[],[],[],[],[],[],@(x)unitNormCon(x),options);
        wt=xf;
        
        if i>1
            delv = abs(fvals_v(i-1)-fvalv);
            delw = abs(fvals_w(i-1)-fvalw); 
        else
            delv=inf;
            delw=inf;
        end
        
        fvals_v(i,1)=fvalv;
        fvals_w(i,1)=fvalw;
        fvals_t(i,1)=gcc(cat(1,vt,wt),blkR1,Rtilde,L); % may be unnecessary
        
        if abs(delv)<optTol && abs(delw)<optTol
            fvals_v(i+1:end,1)=fvalv; % fill in nans with last value
            fvals_w(i+1:end,1)=fvalw;
            fvals_t(i+1:end,1)=fvals_t(i,1); 
            fprintf('Group coordinate descent converged after %d iterations \n',i);
            break
        end
        
        if i==maxIters
             fprintf('Group coordinate descent reached max iterations \n');
             fprintf('Last delv %0.2e and delw %0.ef \n',abs(delv),abs(delw));
        end
        
        
        fprintf('Last value of the objective function = %0.2f \n', fvals_t(i,1) )
    end
    all_fvals_w(:,p)=fvals_w;
    all_fvals_v(:,p)=fvals_v;
    all_fvals_t(:,p)=fvals_t;
    
    Vhat(:,p)=vt;
    What(:,p)=wt;
    
    % regress out driving signal before moving onto next pair
    if p<P
        if ndims(X)==3
            for tr=1:n_trials
                tX = X(:,:,tr);
                tyhat=tX*wt;
                tyhat_p = tplitz2(tyhat,L); % generate convolution matrix
                H = pinv(tyhat_p)*tX; % regress driving signal onto observations
                tX=tX-tyhat_p*H; % remove driving signal from observations
                tXr=tX(end:-1:1,:);
                
                % fill this trial
                X(:,:,tr)=tX;
                Xr(:,:,tr)=tXr;
                %Xr=X(end:-1:1,:); % update Xr
            end
        else
            yhat=X*wt;
            yhat_p = tplitz2(yhat,L); % generate convolution matrix
            H = pinv(yhat_p)*X; % regress driving signal onto observations
            X=X-yhat_p*H; % remove driving signal from observations
            Xr=X(end:-1:1,:); % update Xr
        end
    end
    
end

% tabulate strengths of causality
G_gca=nan(P,1);
yhat=nan(N,P); % driving signal time series
zhat=nan(N,P); % driven signal time series
for p=1:P
    if ndims(X)==3
        yhat=[];
        zhat=[];
        for tr=1:n_trials
            tX = Xog(:,:,tr);
            tyhat=tX*What(:,p);
            tzhat=tX*Vhat(:,p);
            yhat=cat(1,yhat,tyhat);
            zhat=cat(1,zhat,tzhat);
        end
        G_gca(p,1)=vanillaG(yhat,zhat,L);
    else
        yhat(:,p)=Xog*What(:,p);
        zhat(:,p)=Xog*Vhat(:,p);
        G_gca(p,1)=vanillaG(yhat(:,p),zhat(:,p),L);
    end
end

% return some statistics
stats.fvals_v=all_fvals_v;
stats.fvals_w=all_fvals_w; 
stats.fvals_t=all_fvals_t;
stats.options=options;
stats.datestr=datestr(now);
stats.P=P;
stats.L=L;
stats.maxIters=maxIters;
stats.maxcond=maxcond;

end

function fval = gctr_v(v,wo,blkR1,Rtilde,blkR1r,Rtilder,L)
% wo is fixed, v is the independent variable

x=cat(1,v,wo);
fval1 = gcc(x,blkR1,Rtilde,L);

xr=cat(1,wo,v);
fval2 = gcc(xr,blkR1r,Rtilder,L);

fval=fval1+fval2;

fval=real(fval); % ignore small imaginary component

end

function fval = gctr_w(w,vo,blkR1,Rtilde,blkR1r,Rtilder,L)
% vo is fixed, w is the independent variable

x=cat(1,vo,w);
fval1 = gcc(x,blkR1,Rtilde,L);

xr=cat(1,w,vo);
fval2 = gcc(xr,blkR1r,Rtilder,L);

fval=fval1+fval2;

fval=real(fval); % ignore small iamginary component

end

function fval = gcc(x,blkR1,Rtilde,L)

x=x(:); % TODO: check if x is even here 
D=numel(x)/2;
v=x(1:D);
w=x(D+1:2*D);

% compute G causality at (v,w)
tPhi_f=Phi_f(v,w,blkR1,Rtilde,L);
tPhi_r=Phi_r(v,blkR1,Rtilde,L);

%%% original code
fval = tPhi_f/tPhi_r  - 1; % inverted for fmincon

%%% new code
%fval = log(tPhi_f)-log(tPhi_r); 

end

function [res,r,R,R0] = Phi_f(v,w,blkR1,Rtilde,L)
%returns theoretical expression for Phi_f
D=numel(v);

R0=Rtilde(1:D,1:D);

kronILv=kron(eye(L),v);
kronILw=kron(eye(L),w);
kron1Lv=kron(ones(L,1),v);
kron1Lw=kron(ones(L,1),w);

r=kron(eye(2*L),v).' * kron(eye(2),blkR1) * cat(1,kron1Lv,kron1Lw);

R = cat ( 1 , kron(ones(2,1).',kronILv.') ,  kron(ones(2,1).',kronILw.')  ) * ...
    kron(eye(2),Rtilde) * ...
    blkdiag( kronILv , kronILw)  ;

res=v'*R0*v - r.'*(R\r);

end

function [res,q,Q,R0] = Phi_r(v,blkR1,Rtilde,L)
%returns theoretical expression for Phi_r
v=v(:);
D=numel(v);
R0=Rtilde(1:D,1:D);

kronILv=kron(eye(L),v);
kron1Lv=kron(ones(L,1),v);

q = kronILv.'*blkR1*kron1Lv;
Q = kronILv.'*Rtilde*kronILv;

res=v'*R0*v - q.'*(Q\q);

end

function X=tplitz1(x,K)
%X=TPLITZ1(X,K)
% create a convolution matrix from the column vector x, with K time lags
% instead of using circular shift, fill unavailable samples with 0
% do not include intercept term
x=x(:); 
X=toeplitz(x);
uind = logical ( triu ( ones(size(X)) , 1 ) ); 
X(uind)=0;
X=X(:,2:K+1); % only keep up to lag K, and exclude the present sample!
end

function X=tplitz2(x,K)
%X=TPLITZ2(X,K)
% same as tplitz1 but avoids toeplitz() because of memory limitations
% uses for loop instead
x=x(:); 
N=numel(x);
X=nan(N,K);

for n=1:N
    if n==1
        X(n,:)=zeros(1,K);
    elseif n<=K
        nMissing=K-n+1;
        X(n,:)=[x( n-1: -1 : n-(K-nMissing) ); zeros(nMissing,1)].';
    else
        X(n,:)=x(n-1:-1:n-K).';
    end
end

end

function RR= ridgeCov(R,maxcond)
%RR=RIDGECOV(R,MAXCOND) regularize covariance matrix with ridge regression
if nargin<2, maxcond=1000; end

lambdas=eig(R);
lambda1=max(lambdas);
lambdad=min(lambdas);
delta = (lambda1 - lambdad*maxcond)/(maxcond-1);
RR = R + delta*eye(size(R,1)); % this guarantees that cond(RR)=maxcond
end

