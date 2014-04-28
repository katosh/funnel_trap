function chi=pccaplus(P)
s=length(P); % size of membership basis
pi = calcpi(P);
[X D]=eigs(P,min(20,s),'lr'); % right eigenvectors (for conformations)
nC = sum(diag(D)>0.85); % number of clusters
X = X(:,1:nC);
X = orthonormalize(X,pi);
g = guess(X);
f = feasible(g);
A = optimize(f);
chi = X*A;

    function X=orthonormalize(X,pi)
        V=X;
        d=diag(pi);
        for k=1:length(X(1,:))
            for i=1:1:k-1
               X(:,k)=X(:,k)-(X(:,i)'*d*V(:,k))*X(:,i);
            end
            X(:,k)=X(:,k)/sqrt(X(:,k)'*d*X(:,k));  %Normalisieren
        end
    end

    function pi=calcpi(P)
        % calculate equilibrium distribution
        [W E]=eigs(P',1,'lr'); % left eigenvector
        pi=W(:,1);
        pi=pi/sum(pi); % normalize as partition of unity
    end
    

    function Ab0 = guess(X)
        %  Alg3.9 Infeasible initial guess
        %   Alg3.8 Index mapping
        Xorg = X;
        ind=[];
        vt=[];
        for j=1:nC
            rownorm = sqrt(diag(X*X'));
            ind(j) = find(rownorm == max(rownorm),1);

            if j ~= 1
                X = X / rownorm(ind(j));
            end

            vt(j,:) = X(ind(j),:); % v-tilda

            if j == 1
                X = X-repmat(vt(1,:),s,1);
            else
                for i=1:s
                    X(i,:) = X(i,:) - (vt(j,:)*X(i,:)')*vt(j,:);
                end
            end
        end
        %  Alg3.9 infeasible initial guess
        Ab0 = Xorg(ind,:)^-1; % A-bar-0
    end

    function A=feasible(Ab)
        % Alg3.10 feasable transformation matrix
        % Step 1
        for j=2:nC
            Ab(j,1) = -sum(Ab(j,2:end));
        end
        % Step 2
        for i=1:nC
            for l=1:nC
                temp(l) = X(l,2:end)*Ab(2:end,i);
            end
            Ab(1,i) = - min(temp);    
        end
        % Step 3
        A = Ab / sum(Ab(1,:));
    end

    function A = unite(kA,A0)
        A = A0;
        A(2:end,2:end) = kA;
    end
    
    function B=I1(A)
        B=sum(max(X*A)');
    end
    
    function I2(A)
        
    end
    
    function A = optimize(A0)
        PROBLEM.objective = @(A) I1(feasible(unite(A,A0)));
        PROBLEM.x0 = A0(2:end,2:end);
        PROBLEM.options = optimset('MaxIter',100);
        PROBLEM.solver = 'fminsearch';
        A = feasible(unite(fminsearch(PROBLEM),A0));
    end
end
