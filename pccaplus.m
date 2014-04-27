nC=3; % number of clusters
s=len(P)

% calculate equilibrium distribution
[W E]=eigs(P',1,'lr'); % left eigenvector
pi=W(:,1);
pi=pi/sum(pi); % normalize as partition of unity

[X D]=eigs(P,nC,'lr'); % right eigenvectors (for conformations)

% Alg3.11 PCCA
%  Alg3.9 Infeasible initial guess
%   Alg3.8 Index mapping

for j=1:nC
    rownorm = sqrt(diag(X*X'))
    ind(j) = find(rownorm == max(rownorm))
    
    if j ~= 1
        X = X / rownorm(ind(j))
    end    
    vt(j) = X(ind(j),:)
    if j == 1
        X = X-repmat(vt(1),nC,1)
    else
        for i=1:len
            X(i,:) = X(i,:) - (vt(j)*X(i,:)')*vt(j)
        end
    end     
end
