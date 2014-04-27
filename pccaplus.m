function chi=pccaplus(P)
nC=3; % number of clusters
s=length(P);

pi = calcpi(P);
[X D]=eigs(P,nC,'lr'); % right eigenvectors (for conformations)
A = optimize(feasible(guess(X)));
chi = X*A;

    function X=normalize(X)
        X=X % implement normalization
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
        ind=[];
        vt=[];
        for j=1:nC
            rownorm = sqrt(diag(X*X'));
            ind(j) = find(rownorm == max(rownorm));

            if j ~= 1
                X = X / rownorm(ind(j));
            end

            vt(j) = X(ind(j),:); % v-tilda

            if j == 1
                X = X-repmat(vt(1),nC,1);
            else
                for i=1:len
                    X(i,:) = X(i,:) - (vt(j)*X(i,:)')*vt(j);
                end
            end     
        end

        %  Alg3.9 infeasible initial guess
        Ab0 = X(ind,:)^-1; % A-bar-0
    end

    function A=feasible(Ab)
        % Alg3.10 feasable transformation matrix
        % Step 1
        for j=2:nC
            Ab(j,1) = -sum(Ab(j,:)) + Ab(j,1);
        end
        % Step 2
        for i=1:nC
            for l=1:nC
                temp(l) = X(l,:)*Ab(:,i) - Ab(1,i)*X(l,1);
            end
            Ab(1,i) = - min(temp);    
        end
        % Step 3
        A = Ab / sum(Ab(1,:));
    end
    
    function I1(A)
        I1=sum(max((X*A)'));
    end
    
    function I2(A)
        
    end
    
    function optimize(A0)        
        optimize = feasible(fminsearch(@(A) I1(feasible(A)),A0,optimset('MaxIter',100)));
    end
end
