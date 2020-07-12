function [Qi,Aij,Bik]=getQAB(N,K,X)
%Compute HMM initial values
Qi=[1;zeros(N-1,1)];
Aij=toeplitz([0.5 zeros(1,N-1)],[0.5 0.5 zeros(1,N-2)]);
Aij(N,N)=1;
u=round(size(X,2)*((0:N+1)/(N+1)));
for k=1:N
    X1=X(:,u(k)+1:u(k+2));
    v=kmeans(transpose(X1),K);
    for n=1:K
        X2=X1(:,v==n);
        Bik(k,n).m=mean(X2,2);
        Bik(k,n).C=diag(var(transpose(X2)));
        Bik(k,n).w=size(X2,2)/size(X1,2);
    end
end