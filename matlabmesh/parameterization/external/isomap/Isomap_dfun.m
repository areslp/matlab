function out = Isomap_dfun(in)

global Isomap_X

out = Isomap_L2_distance(Isomap_X,Isomap_X(:,in)); 
out = out'; 
