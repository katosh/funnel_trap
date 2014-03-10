% This function is to visulize desity functions written in
% the Matrix chi for the Problem of two linked particles
% in a 11x7 state box.

function visChi(chi)

for i=1:length(chi)
    figure(i)
    P1 = zeros(7,11);
    P2 = zeros(7,11);
    for j=1:77
        % index sets
        I1 = (1:77) + ((j-1)*77);
        I2 = (0:76) * 77 + 1;
        % probabilities
        P1(j) = sum(chi(i,I1));
        P2(j) = sum(chi(i,I2));
    end
    % rescaling for better contrast in plot
    mi = min(min(P1),min(P2));
    ma = max(max(P1),max(P2));
    P1 = (P1 - mi)/(ma - mi);
    P2 = (P2 - mi)/(ma - mi);
    subplot(2,1,1)
    image(P1*255);
    subplot(2,1,2);
    image(P2*255);
end
end
