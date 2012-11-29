function [ hash ] = meshHash( points,normals )
%[ hash ] = meshHash( points,normals )
%   returns integer hash code to identify mesh

ph = round(sum(points)*10);
if exist('normals','var')
    nh = round(sum(normals)*10);
    ph = ph+nh;
end
ph = abs(ph);
hstring = [num2str(ph(1)), num2str(ph(2)), num2str(ph(3))];
hash = str2num(hstring);

end
