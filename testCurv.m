square = reshape(cutv,77,77)';
for i=1:20:77
    pos2 = reshape(square(i,:),11,7)';
    pos2 = pos2 * 255;
    figure(i)
    image(pos2)
end
