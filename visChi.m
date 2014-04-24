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
        reshaped = reshape(chi(i,:),77,77)
        % probabilities
        P1(j) = sum(reshaped(j,:))
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
    S1 = zeros(7,11);
    S2 = zeros(7,11);
    for z=1:7
        for s=1:11
            S1(z,s) = P1(Qco(z,s));
            S2(z,s) = P2(Qco(z,s));
        end
    end
    subplot(2,1,1)
    image(S1*255);
    subplot(2,1,2);
    image(S2*255);
    %colorbar(mi:ma);
end
end
