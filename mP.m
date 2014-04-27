function P=mP(t,Q)
P=expm(t*Q);
norm=sum(P);
norm=norm(1);
P=P/norm;
end
