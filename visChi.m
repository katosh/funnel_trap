% This function is to visulize desity functions written in
% the Matrix chi for the Problem of two linked particles
% in a 11x7 state box.

function visChi(chi)

chi=real(chi');

for i=1:length(chi(:,1))
    figure(i)
    P1 = zeros(1,77);
    P2 = zeros(1,77);
    for j=1:77
        reshaped = reshape(chi(i,:),77,77);
        % probabilities
        P1(j) = sum(reshaped(j,:));
        P2(j) = sum(reshaped(:,j));
    end
    % rescaling for better contrast in plot
    mi = min(min(P1),min(P2));
    ma = max(max(P1),max(P2));
    P1 = (P1 - mi);
    P2 = (P2 - mi);
    if mi ~= ma
        P1 = P1/(ma-mi);
        P2 = P2/(ma-mi);
    end
    % ordering the values into the 11*7 box
    S1 = reshape(P1,11,7)';
    S2 = reshape(P2,11,7)';
    subplot(2,1,1)
    image(S1*255);
    subplot(2,1,2);
    image(S2*255);
    %colorbar(mi:ma);
end
end
