function testCutv(cutv)

% this is to visulize parts of the cutting vector cutv
square = reshape(cutv,77,77)';

% loop over sime positions of P1 to limit the number of plots
for i=37:41
    pos2 = reshape(square(i,:),11,7)';
    pos2 = pos2 * 255;
    figure(i)
    image(pos2)
    p1row = ceil(i/11);
    p1col = mod(i,11);
    title(['Position for P2 when P1 is in row ',num2str(p1row),' and column ',num2str(p1col)])
end
