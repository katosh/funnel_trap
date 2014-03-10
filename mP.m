if ~exist('t')
    t=1;
end
tic
P=expm(t*Q);
toc
