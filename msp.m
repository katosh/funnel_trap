% calculate probability of each metastable state
[W E]=eigs(P',1,'lr'); % left eigenvector
pi=W/sum(W); % normalize as partition of unity
pi'*chi