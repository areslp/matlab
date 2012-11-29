function [ k ] = tripick2( face, i, j )
% [k] = tripick2(face,i,j)  
%   pick [k] != i,j from face
fi = find(face~=i & face~=j);
k = face(fi);
end
