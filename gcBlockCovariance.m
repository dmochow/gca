function [blkR1,Rtilde] = gcBlockCovariance(X,L)
%[BLKR1,RTILDE]=GCBLOCKCOVARIANCE(X,L)
%   generate the spatiotemporal covariance matrices needed for gca
%   optimization
%   X is a samples-by-sensors observation array
%   L is a scalar denoting the max lag in the VAR model
[N,D]=size(X);

% Rtilde
Rtilde=zeros(L*D,L*D);
for i=1:L
    for j=1:L      
        %Rij = 1/N*(X'*circshift(X,[(j-i),0]));
        Xl = circshift(X,[i,0]);
        Xr = circshift(X,[j,0]);
        Xl(1:i,:)=0;
        Xr(1:j,:)=0; 
        Rij = 1/N*(Xl.'*Xr);
%         if j>i
%             Xp=circshift(X,[(j-i),0]);
%             Xp(1:j-i,:)=0; % remove circular portion
%         elseif j<i
%             Xp=circshift(X,[(j-i),0]);
%             Xp(end-(j-i):end,:)=0; % remove circular portion
%         else %j=i
%             Xp=X;
%         end
%            
%         Rij = 1/N*(X'*Xp);
        rowidx=(i-1)*D+1:i*D;
        colidx=(j-1)*D+1:j*D;
        Rtilde(rowidx,colidx)=Rij;
    end
end

% R_1:L
blkR1=[];
for l=1:L
    %blkR1=blkdiag(blkR1,1/N*(X'*circshift(X,[l,0])));
    Xp=circshift(X,[l,0]);
    Xp(1:l,:)=0; 
    blkR1=blkdiag(blkR1,1/N*(X'*Xp));
end
end

