function [z,s]=neig(zeile,spalte,richtung)

switch richtung
case 1
    s=spalte+1;
    z=zeile;
case 2
    s=spalte;
    z=zeile+1;
case 3
    s=spalte-1;
    z=zeile;
case 4
    s=spalte;
    z=zeile-1;
otherwise
    err = MException('ResultChk:OutOfRange', ...
            'this is not a direction');
    throw(err)
end
end
