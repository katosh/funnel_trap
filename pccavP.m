n=3;          %3 metastabile Bereiche, da 3 EW in der Naehe von 1
[V D]=eig(P);
[W,E]=eig(P');
p=W(:,1);
p=p/sum(p);

X=V(:,1:n);
d=diag(p);
for k=1:n
    for i=1:1:k-1
       X(:,k)=X(:,k)-(X(:,i)'*d*V(:,k))*X(:,i);
    end
    X(:,k)=X(:,k)/sqrt(X(:,k)'*d*X(:,k));  %Normalisieren
end

le = 77^2;

Y=X;
ind=zeros(1,n);
norm=zeros(le,1);
step1=true;

for k=1:n
    for j=1:le
        norm(j)=Y(j,:)*Y(j,:)';
    end
    ind(k)=find(norm==max(norm),1);

    if step1
        step1=false;
        Y=Y-repmat(Y(ind(k),:),le,1);
    else
        v=Y(ind(k),:)/sqrt(Y(ind(k),:)*Y(ind(k),:)');
        for i=1:le
            Y(i,:)=Y(i,:)-(v*Y(i,:)')*v;
        end
    end
end

for k=1:n
    B(k,1:n)=X(ind(k),:);
end
A=inv(B);
chi=X*A;
