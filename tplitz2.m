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

