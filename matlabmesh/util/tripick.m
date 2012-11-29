function [ j, k ] = tripick( face, i )
% [j,k] = tripick(face,i)  
%   pick [j,k] != i from face, preserving cw/ccw ordering
   
fi = find(face==i);
face = circshift(face,[0,1-fi]);
j = face(2);
k = face(3);

end
