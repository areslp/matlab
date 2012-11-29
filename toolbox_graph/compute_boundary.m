function bv=compute_boundary(face, options)
% modified, since matlab7.1 same name of function and var is not allowed, so change the var boundary to bv(boundary vertex)

% compute_boundary - compute the vertices on the boundary of a 3D mesh
%
%   boundary=compute_boundary(face);
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
verb = getoptions(options, 'verb', 1);

if size(face,1)<size(face,2)
    face=face';
end

nvert=max(max(face));
nface=size(face,1);

A=sparse(nvert,nvert);
for i=1:nface
    if verb
        progressbar(i,nface);
    end
    f=face(i,:);
    A(f(1),f(2))=A(f(1),f(2))+1;
    A(f(1),f(3))=A(f(1),f(3))+1;
    A(f(3),f(2))=A(f(3),f(2))+1;
end
A=A+A';

for i=1:nvert
    u=find(A(i,:)==1);
    if ~isempty(u)
		fprintf(1,'the first boundary point found!\n');
        bv=[i u(1)]; % edge [i,u(1)] is contained only in one face
        break;
    end
end

if ~exist('bv','var')
	bv=[];
	return;
end

s=bv(2);
i=2;
while(i<=nvert)
    u=find(A(s,:)==1);
    if length(u)~=2
        warning('problem in boundary');
    end
    if u(1)==bv(i-1)
        s=u(2);
    else
        s=u(1);
    end
    if s~=bv(1)
        bv=[bv s];
    else
        break;
    end
    i=i+1;
end
       
if i>nvert
    warning('problem in boundary');
end


%%% OLD %%%
function v = compute_boundary_old(faces)

nvert = max(face(:));
ring = compute_vertex_ring( face );

% compute boundary
v = -1;
for i=1:nvert   % first find a starting vertex
    f = ring{i};
    if f(end)<0
        v = i;
        break;
    end
end
if v<0
    error('No boundary found.');
end
bv = [v];
prev = -1;
while true
    f = ring{v};
    if f(end)>=0
        error('Problem in boundary');
    end
    if f(1)~=prev
        prev = v;
        v = f(1);
    else
        prev = v;
        v = f(end-1);
    end
    if ~isempty( find(bv==v) )
        % we have reach the begining of the boundary
        if v~=bv(1)
            warning('Begining and end of boundary doesn''t match.');
        else
            break;
        end
    end
    bv = [bv,v];
end