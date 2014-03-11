%%% the no be area %%%
nobe=[  1,1;
        1,2;
        1,6;
        1,7];

%%% the no go directions %%%
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

%%% generate Q without limits (nobe & nogo) or diagonal %%%
% base rate a
a = 1;
% there are 77 posible states for a particle in the funnel trap
Q = zeros(77);
% filling Q with naive rates (no respect to limits)
for i=1:7
    for j=1:11
        for k=1:4
            [nz,ns] = neig(i,j,k);  % eighbor field
            if 1 <= nz && nz <= 7 && 1 <= ns && ns <= 11
                Qzeile = Qco(i,j);     % for rate from this field
                Qspalte = Qco(nz,ns);  % for rate to this field
                Q(Qzeile, Qspalte) = a;
            end
        end
    end
end

%%% zero rates for nogo directions %%%
for i=1:length(nogo)
    % from
    Qzeile = Qco(nogo(i,1),nogo(i,2));
    % to
    [nz,ns] = neig(nogo(i,1),nogo(i,2),nogo(i,3));
    Qspalte = Qco(nz,ns);

    Q(Qzeile, Qspalte) = 0;
end

%%% no rates for nobe area %%%
for i=1:length(nobe)
    Qzeile = Qco(nobe(i,1),nobe(i,2));
    Q(Qzeile,:) = zeros(1,77);
end

%%% fill diagonal %%%
for i=1:77
    Q(i,i)=-sum(Q(i,:));
end

%%% preparing the kronecker sum (Q=Q1+Q2) %%%
I=eye(77);
% rates for Particle 1
Q1 = kron(Q,I);
% rates for Particle 2
Q2 = kron(I,Q);

%%% intreducing the rope for Particle 2 %%%
rl = 3; % rope length
for i=1:77 % line number ^= start point
    for j=1:77 % collum number ^= the state we ate going to
        
        
        %%% ----------- states not to go to ----------- %%%
        %{
        ------------------------------------------------- 
        We now want to checkt the posible position for P2
        to go to when P1 goes to position j. We do that by
        checking to which neighbors of j the particle could
        go next by looking at the rates in Q wich are
        restricted to posible movments. Then the process is
        repeated on the neighbors for rl-times.
        This indicates the rope length restricts the
        distance between the two particles to rl movement
        steps. We also use negative rates to also count
        backward steps from P2 to P1 as a possible connection.
        -------------------------------------------------
        %}
        % P holds the posible Positions for P2
        P=[j]; % start with going to j
        for k=1:rl
            nP=[]; % new Positions to be added later
            for p=P
                nP=[nP find(Q(p,:)~=0)];
            end
            P=[P nP];
        end
        P = sort(unique(P));
        % generate going cutter vector (no to go to)
        gocutv = ones(1,77);
        gocutv(P) = 0; % using P as index set


        %%% ---------- states not to be at ----------- %%%
        % P holds the posible Positions for P2
        P=[i]; % start with beeing at i
        for k=1:rl
            nP=[]; % new Positions to be added later
            for p=P
                nP=[nP find(Q(p,:)~=0)];
            end
            P=[P nP];
        end
        P = sort(unique(P));
        % generate going cutter vector (no to go to)
        becutv = ones(1,77);
        becutv(P) = 0; % using P as index set


        %%% generate cutter matrix
        cut = repmat(gocutv,77,1) .* repmat(becutv',1,77);
        % cutting out unwanted rates
        zrange = (i-1)*77+1:i*77; % row range
        srange = (j-1)*77+1:j*77; % colum range
        Q2(zrange,srange) = Q2(zrange,srange) .* cut;
        % fixing the diagonal for mass conservation
        for k=0:76
            z = zrange(1) + k;
            s = srange(1) + k;
            Q2(z,s) = Q2(z,s)-sum(Q2(z,srange));
        end
    end
end

%%% kronecker sum %%%
Q = Q1 + Q2;

%%% checking the diagonal %%%
for i=1:77^2
    if Q(i,i) ~= -sum(Q(i,:)) + Q(i,i)
        display(['Error rates at line ',num2str(i)])
    end
end

Q = sparse(Q);

%%%% cut zero lines %%%
%for i=1:77^2
%    if sum(Q(i,:)<=0.1)<=0.1
%        Q(i,:)=[];
%        Q(:,i)=[];
%    end
%end
