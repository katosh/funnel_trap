P=exp(t*Q);
sum(sum(P>0))
mi = min(min(P))
ma = max(max(P))
image(255*P/mi*0.5);
