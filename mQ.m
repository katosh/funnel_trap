%%% the no be area %%%
% format: row,column
nobe=[  1,1;
        2,1;
        6,1;
        7,1];

%%% the no go directions %%%
% format: row,column,direction(clockwise)
nogo=[  % dont go to nobe
        1,2,3;
        2,2,3;
        3,1,4;
        5,1,2;
        6,2,3;
        7,2,3;

        % dont cross the line
        1,7,1;1,8,3; % upper part of line
        1,7,2;2,7,4;
        2,6,1;2,7,3;
        2,6,2;3,6,4;
        3,5,1;3,6,3;
        3,5,2;4,5,4;
        5,5,4;4,5,2; % lower part of line
        5,5,1;5,6,3;
        6,6,4;5,6,2;
        6,6,1;6,7,3;
        7,7,4;6,7,2;
        7,7,1;7,8,3];

%%% lures and stinky spotts
% format: column,row,factor for rates
lures=[ 4,1,0.1;
        4,5,10];

%%% generate Q without limits (nobe & nogo) or diagonal %%%
% base rate a
a = 1;
% there are 77 posible states for a particle in the funnel trap
Q = zeros(77);
% filling Q with naive rates (no respect to limits)
for i=1:7
    for j=1:11
        for k=1:4
            [nz,ns] = neig(i,j,k);  % neighbor field
            if 1 <= nz && nz <= 7 && 1 <= ns && ns <= 11
                Qzeile = Qco(i,j);     % for rate from this field
                Qspalte = Qco(nz,ns);  % for rate to this field
                Q(Qzeile, Qspalte) = a;
            end
        end
    end
end

%%% zero rates for nogo directions %%%
for i=1:length(nogo(:,1))
    % from
    Qzeile = Qco(nogo(i,1),nogo(i,2));
    % to
    [nz,ns] = neig(nogo(i,1),nogo(i,2),nogo(i,3));
    Qspalte = Qco(nz,ns);
    Q(Qzeile, Qspalte) = 0;
end

%%% no rates for nobe area %%%
for i=1:length(nobe(:,1))
    Qzeile = Qco(nobe(i,1),nobe(i,2));
    Q(Qzeile,:) = zeros(1,77);
end

%%% fill diagonal %%%
for i=1:77
    %Q(i,i)=-sum(Q(i,:));
    Q(i,i) = a; % difficult explanation
end

%{
%%% lures and stinky spotts %%%
for i=1:length(lures(:,1))
    factorM=ones(77,77);
    co=Qco(lures(i,1),lures(i,2));
    factorM(co,:)=lures(i,3);
    factorM(:,co)=lures(i,3);
    Q = Q .* factorM;
end
%}

%%% preparing the kronecker sum (Q=Q1+Q2) %%%
I=eye(77);
% rates for Particle 1
Q1 = kron(Q,I);
% rates for Particle 2
Q2 = kron(I,Q);

%%% intreducing the rope / eliminate impossible states %%%
rl = 3; % rope length
cutv = zeros(1,77^2); % cutting vector
% looking for possible combination of states
for i=1:77
    % P1 .. Particle 1
    % P2 .. Particle 2
    % for P1 in state i, P2 can be in the stats listed in P
    P=[i];
    for k=1:rl
        nP=[]; % new Positions to be added later
        for p=P
            nP=[nP find(Q(p,:)~=0)];
        end
        P=[P nP];
    end
    P = sort(unique(P)); % not sure if this is necessary
    range = 1:77 + ((i-1)*77); % state i for P1
    % range(P) are all the posible states of the system
    % adding possible states
    cutv(range(P))=1;
end
%%% kronecker sum %%%
Q = Q1 + Q2;
%%% generate cutter matrix
cut = repmat(cutv,77^2,1) .* repmat(cutv',1,77^2);
% cutting out unwanted rates
Q = Q .* cut;
% fixing the diagonal for mass conservation
for i=1:77^2
    Q(i,i) = Q(i,i)-sum(Q(i,:));
end
% done
Q = sparse(Q);
