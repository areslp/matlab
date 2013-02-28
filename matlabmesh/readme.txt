MATLAB mesh toolkit
Ryan Schmidt  rms@dgp.toronto.edu or rms@unknownroad.com
http://www.dgp.toronto.edu/~rms/matlabmesh

Contents:
 - What is this?
 - Installation / Setup
 - Data Structures
 - 2D/3D Geometry
 - Mesh and Point Set Support
 - Mesh Fairing
 - Planar Parameterization
 - Mesh Deformation
 - Miscellany
 - Sample Meshes and Point Sets
 - Mex
 - Credits, Copyright Notices, etc
 - Citations

========================
  WHAT IS THIS??
========================
 
This toolkit is a sandbox for mesh and point set processing in Matlab. 

Matlab is a mixed bag for geometry processing. The built-in libraries provide lots of
great support if you want to do things like solve linear systems and minimize energies, 
and since Matlab is so ubiquitous in machine learning it is easy to get code
for all sorts of complicated things. The interactive environment is also great
for rapid prototyping and experimentation.

Unfortunately, matlab is not particularly suitable for working with large 
arbitrary graphs (ie meshes). Explicit loops and branching are ridiculously slow,
so efficient code often takes the form of cryptic, "write-only" statements that
only matlab gurus will be able to understand.

Even text parsing is abysmally slow. Perhaps the most useful feature of this
package is that the OBJ loader automatically creates binary caches. The first time
you load your OBJ will be slow, after that it will be instantaneous. The package
also supports useful basic mesh operations like:
	- extracting a submesh copy with a map back to the base mesh
	- plotting meshes in various ways
	- picking faces and vertices
	- interrogating meshes (one rings, areas, edges, boundary loops, etc)

I use this package for prototyping and experimentation, so my goal has been
clarity over efficiency. To that end, I have written many convenience functions 
which, in my opinion, make the code vastly more readable. For example,
you can compute the magnitude of vector X using:
     vmag(X)
rather than the more cryptic (and longer):
     sqrt(sum(X.^2))
And vmag() automatically computes per-row magnitudes if X is a list.

Similarly, adding a vector v to one or more points X is
	 vadd(X,v)
instead of
	 X + repmat(v,size(X,1),1)
 
Where possible, the utility functions will automatically handle lists
of elements, arbitrary dimension, and sometimes even ordering (for
exmaple the matrix-vector multiply mulv produces the same output vector 
regardless of whether v is a row or column vector)
 
In terms of higher-order geometry processing, there are implementations of
various algorithms presented in the literature, such as:

	- vertex weights (cotangent, authalic, mean-value, voronoi area-weighting,
					  and my own unpublished techniques based on discrete exponential maps)
	- laplace-beltrami operators and analysis (based on above)
	- planar parameterization (MIPS, standard fixed-boundary linear techniques,
			natural conformal, spectral conformal, variants of dimensionality
			reduction techniques from machine learning, discrete exponential map,
			and my own combinations of these techniques)
	- mesh deformation (Laplacian, Poisson, Rotation-Invariant Coordinates)
	- mesh fairing (implicit laplacian, second-order mean curvature)

See the documentation below for specifics and references.	

Various sample meshes and point sets are included in the \data
subdirectory, and lots of sample code in the various test files
and the misc\ subdirectory (which contains exploratory code that
was the basis for some unpublished work).

There is mex used here and there, and some external dependencies used 
in some of the sample/testing code. You may have to jump through some
hoops to get those to work on your computer. Good luck!
 
 
========================
  INSTALLATION / SETUP
========================

This toolkit contains several subdirectories, and will not work
properly unless they are all in your paths. The best way to deal
with this is to have matlab automatically add the necessary paths
on startup. 

On windows, you can automatically add the paths by inserting
the following commands in the startup file stored at
<Your Documents>\MATLAB\startup.m (create if it does not exist):

   addpath(genpath('c:\path\to\matlab'));

(where 'c:\path\to\matlab' is the path to the package root). 

If you do not have it already, I strongly recommend downloading
Tom Minka's excellent lightspeed package:

http://research.microsoft.com/~minka/software/lightspeed/


===================
  DATA STRUCTURES
===================

Matlab doesn't really do explicit data structures. I use a
struct for meshes (see below) but vectors are just arrays.
The code in this package assumes that:

	- Points and Vectors are (1xd) matrices
	- lists of points/vectors are (Mxd) matrices  (M rows x d columns)

So, stick with that.

==============================
  2D/3D GEOMETRY & UTILITIES
==============================

util\ 

	%% basic geometric computations

	[angleDist.m]
		- angular distance between two 2D points
	[argN.m]
		- return Nth returned argument of f(varargin)
	[cart2sphM.m]
		- cartesian/spherical-coords conversion consistent with Mathematica
	[cart2sphY.m]
	[sph2cartY.m]
		- cartesian/spherical-coords conversion consistent with mathworld, graphics papers, etc
	[clamp.m]
		- clamp value to range [min,max]
	[cross2.m]
		- 2D cross product x1*y2 - x2*y1
	[errnorm.m]
		- compute scaled error using p-norm
	[findKNN.m]
		- find k-nearest nbrs for each input point
	[genpoints.m]
		- generate random points in 2D using varions strategies:
			%       'grid_unit_square' - regular-grid in [0,1]^(dim==2)
			%       'stratified_unit_square' - jittered 'grid_unit_square'
			%       'uniform_unit_square' - uniform random samples in [0,1]^dim		
	[mvmul.m]
		- multiply row-vectors v by matrix M and return row-vector
	[ncross.m]
		- compute normalized cross-product
	[ngrad.m]
		- compute numerical gradient approximation
	[normalize.m]
		- normalize rows of M. If magnitude is 0, skip vector
	[perpdot.m]
		- compute 2D perpendicular vector (-y,x)
	[tangentFrame.m]
		- create orthogonal tangent vectors to normal
	[vadd.m]
		- adds v to each row of M
	[vangle.m]
		- compute angle between vectors A and B
	[vcot.m]
		- compute cotan of angle between n-D vectors A and B
	[vdot.m]
		- returns vector of dot products of rows of M and M2. If M2 is a 
		  single row, expands to # of rows in M
	[vmag.m]
		- compute magnitude of rows in M
	[vmag2.m]
		- compute sqiared magnitude of rows in M
	[vmul.m]
		- matrix/vector multiplication (produces same result whether vector
		  is row or column)

	%% triangle geometry
		  
	[triarea.m]
		- area of triangle defined by 3 points
	[triarea2s.m]
		- signed area of triangle defined by 3 2D points
	[trigrad.m]
		- compute barycentric basis function gradients for a triangle
	[trinormal.m]
		- face normal of triangle
	[triobtuse.m]
		- return 1 if any triangle angles are obtuse
	[tripick.m]
		- pick [j,k] != i from face, preserving cw/ccw ordering
	[tripick2.m]
		- pick [k] != i,j from face
		  
		  
	%% transformations
		  
  	[axisrot.m]
		- construct 3x3 rotation matrix around arbitrary axis
	[valign.m]
		- generate matrix that rotates vector vfrom into vector vto

		
	%% matrix constraints (probably should go elsewhere)
   
	[issymmetric.m]
		- test if matrix is symmetric with epsilon tolerance
	[hardConstrain.m]
		- add hard constraints to linear system via rewriting
	[hardConstrainSolve.m]
		- solve linear system M.X-B with hard constraints
	[softConstrain.m]
		- add soft constraints to linear system
   
   
external\
   dijkstra\ 
     [*] matlab and mex implementations of Dijkstra's algorithm

drawing\           
	[draw*.m] convenience 2D drawing routines for circle, line, point, polyline
	[TESTS.m] testing code

geometry2\         
	[circle2_2pr.m] find circle center from 2 points and radius
	[isect_line2_circle2.m] intersect 2D line with circle 
	[turningNumber.m] compute integer turning number of closed polygon
	
geometry3\
	[sphereNormCoords.m] compute analytic (u,v) normal coordinates on sphere

implicit3\
	[implicitPlane.m] signed distance to plane
	[implicitSphere.m] signed distance to sphere
	[thinPlateSpline.m] solve and evaluate 2D sqr(n)log(n) thin-plate spline interpolant
					    for sparse constraints (current code evaluates over grid)
	[test\*.m] testing code for above

===============================
  MESH AND POINT SET SUPPORT
===============================

mesh\
	
	%% input/output/datastructure creation

	[readMesh.m] 
		- Read OBJ meshes with optional binary caching (much faster on second read)
	[readPoints.m] 
		- Read OBJ-style pointset with optional binary caching
	[writeMesh.m]
		- write mesh datastructure in OBJ format
 
	[makeMesh.m] 
		- create mesh from vertex and face lists, optionally with normals and uv's. Output:
			% mesh.v = Nx3 list of vertices
			% mesh.n = Nx3 list of vertex normals
			% mesh.u = Nx2 list of vertex uv's
			% mesh.f = Mx3 list of triangles
			% mesh.fn = Mx3 list of triangle face normals
			% mesh.vidx = Nx1 list of vertex indices (useful for various things)
			% mesh.fidx = Nx1 list of face indices (useful for various things)
			% mesh.bounds = [min;max;center] bounding box
			% mesh.e = symmetric sparse array of edges
			%             mesh.e(i,j) = 1 if boundary edge, 2 if interior edge
			% mesh.te = sparse array of edge triangles
			%             mesh.e(i,j) = tri1, mesh.e(j,i) = tri2
			%             for boundary tri, one of these will be 0...
			% mesh.loops = boundary loops of mesh
			% mesh.valence = valence of each vertex of mesh
			% mesh.isboundaryv = 0/1 flags for boundary vertices
			% mesh.isboundaryt = 0/1 flags for boundary triangles
 
	[makePointSet.m] 
		- create point set from list of points, with optional normals, uvs, colors. Output:
			% pointset.v = Nx3 list of vertices
			% pointset.n = Nx3 list of vertex normals
			% pointset.u = Nx2 list of vertex uv's
			% pointset.vidx = Nx1 list of vertex indices (useful for various things)
			% pointset.bounds = [min;max;center] bounding box
			% pointset.e = symmetric sparse array of edges
			%             pointset.e(i,j) = |v_j-v_i| if connected, 0 otherwise
			% pointset.valence = valence of each vertex of pointset
  
	[subMesh.m] 
		- create submesh from face roi (also returns vertex and face maps)
	[clipEars.m]
		- remove 'ear' triangles from mesh (ie tris with only 1 nbr)
	[orientMesh.m]
		- generate consistent orientation for 2D mesh triangles	
		
	[meshHash.m]
		- compute integer hash code to identify mesh
				(this is a hack and probably not very good hash...)
		
	%% display
		
	[plotMesh.m] 
		- display mesh data in many ways based on flag set:
			%   'v' = draw vertex points
			%   'e' = draw triangle edges
			%   'b' = draw boundary edges
			%   'n' = draw normals
			%   'f' = draw filled faces
			%   'l' = apply lighting to faces
			%   'u' = use uv-values as vertices
			%   'U' = use uv-values as vertices, and try to consistently-orient/scale plots
			%   'i' = plot vertex indices (numbers)  ** WARNING: very slow 
			%   'O' = pretty-plot for output (thicker lines, etc)
			
	[plotMeshGL.m] 
		- display mesh using external viewer program view3D.exe (in external)
	
	[plotPoints.m]
		- display pointset in various ways based on flag set:
			%   'v' = draw vertex points
			%   'n' = draw normals
			%   'u' = use uv-values as vertices
			%   'U' = use uv-values as vertices, and try to consistently-orient/scale plots
			%   'i' = plot vertex indices (numbers)  ** WARNING: very slow 
			%   'O' = pretty-plot for output (thicker lines, etc)	
	
	[plotPolyline.m]
		- display 2D/3D polyline
	
		   
	%% selection/picking	   
	
	[selectMesh.m] 
		- select vertices or faces of mesh using implicit function
	[pickVertex.m] 
		- interactive selection of vertex of mesh
	[nearestVertex.m] 
		- finds nearest vertex to input point


	%% interrogation

	[estimateNormal.m]
		- estimate normals using given strategy ( face average, area-weighted face avg )
	[faceNormal.m]
		- compute face normals
	[faceArea.m]
		- compute area of mesh faces
	[meshArea.m]
		- compute face area of ROI of mesh
	[faceColors.m]
		- compute mapped color values for mesh faces
	[edgeStats.m]
		- compute min/max/avg edge lengths for given connectivity graph
	[onering.m] 
		- find vertex one-ring vertices and triangles (supports ccw sorting)
	[oneringf.m]
		- return indices of one-ring faces at nVertex
	[oneringv.m]
		- return indices of one-ring vertices at nVertex
	[findBoundaries.m]
		- extract set of boundary loops from estimated set of 
		  boundary vertices, using either 3D or UV coordinates (right?)


	%% approximate/differential properties

	[vertexArea.m]
		- compute area for a vertex on mesh using various methods:
			%     'uniform' - area = 1   
			%     'onering' - sum of areas of one-ring triangles
			%     'voronoi' - voronoi areas from [Meyer02]  (invalid for obtuse triangles)
			%     'mixed'   - mixed voronoi/geometric areas from [Meyer02]		
	[pointArea.m]
		- estimate area for a point of a pointset using various methods:
			%     'uvVoronoi' - voronoi cell area in UV-space
			%     'uvDelArea' - full delaunay area in UV-space
			%     'uvDelAreaR' - delaunay area within radius
			%     'uvDelRing' - delaunay one-ring area in UV-space
		  
	[cotanWeights.m] 
		- one-ring cotangent weights, optinally compute authalic and/or
			triangle area-weighted [Desbrun02, Meyer02, Mullen08]
	[meanvalueWeights.m]
		- one-ring mean value weights [Floater04]
	[laplacian.m] 
		- compute vertex laplacian vectors (using uniform or cotangent weights)
		
	[meanCurv.m] 
		- Mean Curvature estimation techniques, including [Moreton92] Mean Curvature Estimation 
		  w/ [Schneider01] improvement for low-valence vertices, and [Meyer02] discrete cotangent

	   
=====================
  MESH FAIRING
=====================
	   
fairing\
	[fair_LBH.m]  
		- implmentation of [Schneider01] intrinsic mesh fairing
	[fair_Laplacian.m]
		- implicit fairing by solving for zero-length Laplacian vectors
		    (ie minimal mean-curvature surface)
	[fair_LFlow.m]
		- explicit Laplacian fairing (currently using uniform weights)
    [TESTS.m]
	    - tests for above functions
 
=====================
  MESH DEFORMATION
=====================

deformation\
	[deformLaplacian.m] 
		 - solve for laplacian deformation problem [Lipman04], given dense 
		   per-vertex rotations and sparse soft vertex position constraints
	[deformPoisson.m]
		 - solve Poisson deformation problem, given dense per-vertex 
		   rotations and sparse *hard* vertex constraints
	[deformRogInvCoords.m] 
		 - solve Rotation-Invariant Coordinates deformation problem [Lipman05], 
		   given sparse soft vertex position and orientation constraints
	[meshDeform_Lipman04.m]  
		- [Lipman04] Laplacian mesh deformation (iterative rotation estimation)
	[meshDeform_Laplacian.m]
		 - simple Laplacian deformation with soft vertex position constraints,
		   with optional postprocess to align to target normal field
	[meshDeform_Yu04.m]
		 - (partial, possibly incorrect) implementation of [Yu04] Poisson mesh deformation.
			input is mesh, boundary constraint points, and transformed versions
			of input triangles

	TESTS\*.m
		- testing/sample code for above deformations
	
		
=============================
  PLANAR PARAMETERIZATION
=============================		
		
parameterization\

	%% embeddings

	[embedMIPS.m]
		- parameterize mesh using MIPS iterative nonlinear technique [Hormann00]
	[embedDNCP.m]
		- compute natural conformal parameterization (DNCP) from [Desbrun02]
	[embedSCP.m]
		- compute spectral parameterization using various techniques
		  (basically [Mullen02] with some variants...)
	[embedPSCP.m]
		- variant of [Mullen02] that can be applied to point sets
	[embedLLE.m]
		- spectral parameterization based on Locally-Linear Embedding (LLE) [Roweis02]
	[embedLEM.m]
		- spectral parameterization based on Laplacian Eigenmaps (LEM) [Belkin03]
	[embedBoundary.m]
		- embed boundary of mesh in 2D polygon (currently only circle via arc-length)
	[embedInterior.m]
		- solve for parameterization given boundary constraints and weight matrix,
		  with optional additional position constraints
		  
	  
	 %% vertex weights
	  
	[makeEpsBallWeights.m]
		- compute per-vertex weights based on euclidean ball nbrhoods
				(uniform, invdist, invdist2, gaussian heat-kernel)
	[makeExpMapWeights.m]
		- compute per-vertex weights based on geodesic/DEM nbrhoods
			(uniform, uniform on connected tangent-space delaunay,
			 invdist, [Roweis02] optimal, optimal in UV,
			 heat-kernel gaussian, heat-kernel gaussian in UV)
	[makeKNbrWeights.m]
		- compute per-vertex weights based on K-nbrhoods
			(uniform, invdist, invdist2, heat-kernel gaussian)
	[makeOneRingWeights.m]
		- compute per-vertex weights based on mesh one-ring nbrhoods
			(uniform, invdist, invdist^2, [Desbrun02] discrete conformal,
			 [Desbrun02] discrete authalic, [Floater04] mean value,
			 [Roweis02] optimal, optimal in UV) with optional area-weighting
			 (uniform, one-ring area, [Meyer02] mixed voronoi)
			
	%% exponential-map tangent space coordinates
		  
	[expmap.m]
		- compute Discrete Exponential Map (DEM) [Schmidt06] with upwind-average
		  extension [Schmidt10]
	[fastExpMaps.m]
		- compute DEM for each mesh vertex using external binary expmapCL.exe
	[oneringExpMap.m]
		- compute tangent space coords for one-ring of mesh vertex based
		  on angle mapping
		
	%% distortion metrics
		
	[eConformal.m] 
		- compute conformal, dirichlet energy of mesh+vtxweights
	[faceDistortion.m]
		- compute per-face dirichlet or quasi-conformal distortion
	[triDistortion.m] 
		- compute uv angle/area distortion for triangles
	[triStretch.m] 
		- compute [Sander02] per-triangle uv stretch metrics L2,Linf
	[meshDistortion.m]
		- compute area, angle, L2, Linf metrics for mesh
	[FnOptLaplacianScale.m]
		- compute conformal energy based on weightmatrix and paramtype params
	
	
	%% sample/testing code
	
	[TESTS\*.m] various testing code
	

	%% external algorithms (mainly machine-learing dimensionality reduction codes)
	
	external\
		[csdp.exe] constrained semidefinite problem solver
		[lle.m] code from [Roweis02] for LLE
		[ltsa.m] implementation of LTSA (Zhenyue Zhang & Hongyuan Zha, 2004)
					(http://epubs.siam.org/sam-bin/dbq/article/41915)
		[mani.m] manifold learning demonstration GUI
				 by Todd Wittman, Department of Mathematics, University of Minnesota
				 http://www.math.umn.edu/~wittman/mani/index.html
		[LEigenmaps.m] function for Laplacian eigenmap
					   Written by Belkin & Niyogi, 2002
		[HessianLLE.m] implementation of HLLE
					   Written by David Donoho & Carrie Grimes, 2003
		isomap\
			[Isomap.m]   implementation of Isomap algorithm [Tenenbaum00]
				- uses Floyd's algorithm to estimate geodesics
			[IsomapII.m] implementation of sparse landmark-Isomap algorithm  
				- uses Dijkstra's algorithm if available (in util\external)
			[*] support files
		mvu\
			[mvu.m]  implementation of Maximum-Variance Unfolding (MVU)
			[*] support files
		expmapCL\
			- binary command-line DEM computation. supplied pre-compiled for win32,
			  if you need to rebuild or port, download source from:
			  http://www.dgp.toronto.edu/~rms/software/expmapdemo.html
			

  

===================
  MISCELLANY
===================

config\
  [globalConfig.m] 
     global settings for versions/paths/etc
	     'meshVersion': version for mesh cache format
		 'pointsetVersion': same, for pointset cache format
		 'cachePath': place cache files can be stored (current config\cache)
		 'enableMeshCaching': default 1
		 'enableWeightMatrixCaching': default 0

misc\
  [EG09_WeightTests.m] analysis of vertex weight distributions
  [EG10_LNLE.m] various techniques for spectral parameterization of meshes & point clouds
        ( most of the code used in [Schmidt09] )
  [SGP09.m] various embedding strategies using DEM-weights (precursor to files above)
  
  [SGP10.m] analysis of eigenspectrum of laplace-beltrami operators,
            comparison of DEM vs analytic normal coords on hemisphere
  [SGP10_expmap.m] more comparison of DEM vs analytic normal coords on hemisphere
  
  [SGP09_LBTests.m] analysis of various mesh/pointset laplace-beltrami operators on plane
  [SGP10_lb.m] analysis of various mesh/pointset laplace-beltrami operators on sphere
    **NOTE** above files contain analytic l-b functions for various functions on plane/sphere,
	         which I got via mathematica. Code for generating these solutions is in
			 _background\analytic_lb.nb

  [figures\*.m] Code to generate the figures in the tech report [Schmidt09]
       
       

	   
=============================
SAMPLE MESHES AND POINT SETS
=============================

data\
    [*.obj] misc OBJ files used in testing code, etc
	[armadilloman_tiny.obj] 75-vertex version of armadillo
	[bunny*.obj] various versions of bunny mesh
	[hemisphere_graphite_*.obj] various near-regular hemispheres (created in Graphite)
	[sphere_graphite_*.obj] various near-regular spheres (created in Graphite)
	[sphere_semireg_*.obj] various semi-regular spheres (created by smoothing marching cubes)
  	pointsets\
	    [*.pobj]
            various decimated/edited versions of surfel sample files from
            http://graphics.ethz.ch/pointshop3d/download.html (in OBJ format)		
		[*.pl]
		    various perl scripts for converting from .surfel format, decimating point set, etc
			 
			    
===========
  MEX
===========

You may need to re-mex some things to get this stuff
to run on your computer/version-of-matlab:

[util\external\dijkstra]
  mex dijkstra_maxdist.cpp -largeArrayDims

			 


======================================
  CREDITS / COPYRIGHT NOTICES / ETC
======================================

Isomap code    [parameterization\external\isomap\*]
   - Author: Josh Tenenbaum (jbt@psych.stanford.edu)
   - (c) 1998-2000 Josh Tenenbaum
   - Website for updates: http://isomap.stanford.edu\

Fibonacci heap   [util\external\dijkstra\(fibheap.h,dijkstra.cpp)]
   - Copyright (c) 1996 by John Boyer
   - updated by Ryan Schmidt, Feb 2009, to compile under Visual Studio 2005

Dijkstra's algorithm  [util\external\dijkstra\dijkstra.*]
   - Copyright (c) Mark Steyvers, Stanford University, 12/19/00

Dijkstra's algorithm  [util\external\dijkstra\dijk.m]
   - Copyright (c) 1998-2000 by Michael G. Kay
   - Modified by JBT, Dec 2000, to delete paths

MVU Code    [parameterization\external\mvu\*]
   - mvu.m copyright (c) 2004 by Kilian Q. Weinberger   ( http://www.weinbergerweb.net/Downloads/MVU.html )
   - mvu_writespda.m   copyright (c) 2006 Brian Borchers, borchers@nmt.edu
   - csdp.exe copyright (c) Brian Borchers ( http://infohost.nmt.edu/~borchers/csdp.html )

   


======================================
  CITATIONS
======================================

[Belkin02]
	Laplacian Eigenmaps for Dimensionality Reduction and Data Representation
	M. Belkin, P. Niyogi
	Neural Computation, June 2003; 15 (6):1373-1396.
	http://www.cse.ohio-state.edu/~mbelkin/papers/LEM_NC_03.pdf

[Desbrun02]
	Mathieu Desbrun, Mark Meyer, Pierre Alliez
	Intrinsic Parameterizations of Surface Meshes. 
	Comput. Graph. Forum 21(3): (2002)
	http://www.geometry.caltech.edu/pubs/DMA02.pdf

[Floater04]
	M. S. Floater
	Mean value coordinates
	Comp. Aided Geom. Design 20 (2003), 19-27
	http://heim.ifi.uio.no/~michaelf/papers/mean_value.pdf
	
[Hormann00]
	K. Hormann and G. Greiner
	mips: an efficient global parametrization method
	http://www.inf.usi.ch/hormann/papers/Hormann.2000.MAE.pdf

[Lipman04]
	Yaron Lipman, Olga Sorkine, Daniel Cohen-Or, David Levin, Christian Roessl, Hans-Peter Seidel
	Differential Coordinates for Interactive Mesh Editing, SMI 2004
	http://igl.ethz.ch/projects/Laplacian-mesh-processing/Laplacian-mesh-editing/diffcoords-editing.pdf
   
[Lipman05]
	Yaron Lipman, Olga Sorkine, Daniel Cohen-Or, David Levin
	Linear Rotation-Invariant Coordinates for Meshes, ACM SIGGRAPH 2005
	http://igl.ethz.ch/projects/Laplacian-mesh-processing/linear-rotinv-coordinates/index.php
	
[Meyer02]
	MEYER, M., DESBRUN, M., SCHRÖDER, P., AND BARR, A. H. (2002). 
	Discrete Differential-Geometry Operators for Triangulated 2-Manifolds. 
	http://multires.caltech.edu/pubs/diffGeoOps.pdf
	
[Moreton92]
	Henry Packard Moreton and Carlo H. Sequin
	Functional Optimization for Fair Surface Design
	SIGGRAPH 1992
	http://www.cs.berkeley.edu/~sequin/PAPERS/SIGGR92_MVC_MVS.pdf

[Mullen08]
	Patrick Mullen, Yiying Tong, Pierre Alliez, Mathieu Desbrun
	Spectral Conformal Parameterization 
	Symposium of Geometry Processing, 2008
	http://www.geometry.caltech.edu/pubs/MTAD08.pdf

[Roweis02]
	Nonlinear dimensionality reduction by locally linear embedding.
	Sam Roweis & Lawrence Saul.
	Science, v.290 no.5500 , Dec.22, 2000. pp.2323--2326.
	http://www.cs.nyu.edu/~roweis/lle/
   
[Sander02]
	Sander, Pedro V. and Snyder, John and Gortler, Steven J. and Hoppe, Hugues
	Texture mapping progressive meshes, SIGGRAPH 2001
	http://dl.acm.org/citation.cfm?id=383307
	http://research.microsoft.com/en-us/um/people/hoppe/tmpm.pdf   

[Schmidt06]
	Ryan Schmidt, Cindy Grimm, Brian Wyvill. 
	Interactive Decal Compositing with Discrete Exponential Maps (2006). 
	ACM Transactions on Graphics (SIGGRAPH 2006), 25(3), July 2006, pp. 605-613.
	http://www.unknownroad.com/publications/ExpMapSIGGRAPH06.pdf
	http://www.dgp.toronto.edu/~rms/software/expmapdemo.html
	
[Schmidt09]
	Ryan Schmidt, Karan Singh 
	Approximate Conformal Parameterization of Point-Sampled Surfaces (2009)
	Technical Report CSRG-605, Department of Computer Science, University of Toronto 
	http://papervideos.s3.amazonaws.com/PointSetParam09.pdf
	
[Schmidt10]
	Ryan Schmidt.
	Part-Based Representation and Editing of 3D Surface Models (2010). 
	PhD Thesis, University of Toronto
	http://www.dgp.toronto.edu/~rms/pubs/PhDThesis10.html
	
[Schneider01]
	Robert Schneider, Leif Kobbelt
	Geometric Fairing of Irregular Meshes for Free-Form Surface Design
	Computer Aided Geometric Design 18 (4): 359-379, 2001	
	
[Tenenbaum00]
	J. B. Tenenbaum, V. de Silva, J. C. Langford (2000).  A global
	geometric framework for nonlinear dimensionality reduction.  
	Science 290 (5500): 2319-2323, 22 December 2000.  	
	
[Yu04]
	Yu, Yizhou and Zhou, Kun and Xu, Dong and Shi, Xiaohan and Bao, Hujun and Guo, Baining and Shum, Heung-Yeung
	Mesh editing with poisson-based gradient field manipulation, SIGGRAPH 2004
	http://dl.acm.org/citation.cfm?id=1015774
	http://www.cs.uiuc.edu/~yyz/publication/poisson_mesh-sig04.pdf
	