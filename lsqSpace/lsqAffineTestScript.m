%//////////////////////////////////////////////////////////////////////////////////////////////////////
% To demonstrate usage of some lsq** packaged functions, we draw n random
% 3D points, and plot :
% (1) their least-square-optimal line (in blue),
% (2) their lsq-optimal plane (in grey) -which obviously includes the optimal line,
% (3) their lsq-optimal plane after discarding the 2 worst samples (those with 
% largest deviation). This is the cyan patch.  The samples that survived the 
% discarding have their number explicitly stated in black.
%///////////////////////////////////////////////////////////////////////////////////////////////////////////


%%
n=6;
%%
V = (rand(n,3)-.5) * 20; % generate random 3D point set
clf
%%

[Pl , Dl] = lsqLine(V);

[Pp,Dp] = lsqAffineSpace(V,2) ; 
[Pp2,Dp2, usedsamples] = lsqAffineSpace(V,2,'discardthresh',1);

% visualize points
hold on;
line(V(:,1),V(:,2),V(:,3),'linestyle','none','marker','o');

% visualize line
lDat =[ Pl - 10*Dl; Pl+10*Dl];
line(lDat(:,1),lDat(:,2),lDat(:,3));

% visualize plane
pDat.Vertices = [	Pp + 6*Dp(1,:) - 6*Dp(2,:) ; Pp + 6*Dp(1,:) + 6*Dp(2,:) ; ...
					Pp - 6*Dp(1,:) + 6*Dp(2,:); Pp - 6*Dp(1,:) - 6*Dp(2,:) ];
pDat.Faces = [ 1 2 3 4];
patch(pDat,'FaceAlpha',0.4);

% visualize plane2
pDat2.Vertices = [	Pp2 + 6*Dp2(1,:) - 6*Dp2(2,:) ; Pp2 + 6*Dp2(1,:) + 6*Dp2(2,:) ; ...
						Pp2 - 6*Dp2(1,:) + 6*Dp2(2,:); Pp2 - 6*Dp2(1,:) - 6*Dp2(2,:) ];
pDat2.Faces = [ 1 2 3 4];
patch(pDat2,'FaceAlpha',0.4,'FaceColor','c');

ns = numel(usedsamples);
stringcells = mat2cell(usedsamples,1,ones(ns,1));
ht = text(V(usedsamples,1)+1,V(usedsamples,2),V(usedsamples,3), stringcells);

% full point strings
% stringcells = mat2cell(1:n,1,ones(n,1));
% ht = text(V(:,1)+1,V(:,2),V(:,3), stringcells);

axis vis3d auto equal
view(3)