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
        5,5,4;4,5,2;
        5,5,1;5,6,3; % lower part of line
        6,6,4;5,6,2;
        6,6,1;6,7,3;
        7,7,4;6,7,2;
        7,7,1;7,8,3];

%%% generate Q without limits (nobe & nogo) or diagonal %%%
a = 1;
Q = zeros(77);
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
    Qzeile = Qco(nogo(i,1),nogo(i,2));
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
