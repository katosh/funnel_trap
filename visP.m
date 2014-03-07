mi = min(min(P))
ma = max(max(P))
VP=(P-mi)/(ma-mi);
me = median(median(P))
image(255*VP);
