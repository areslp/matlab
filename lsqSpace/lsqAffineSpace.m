function [Point, DirVecs,varargout] = lsqAffineSpace(V, SpaceDim,varargin)
    % [Point, SpanVecs] = lsqAffineSpace(V, SpaceDim)
    %
    % Given an n x m matrix V, whose rows are a set of sample m-vectors,
    % the routine calculates an affine subspace of R^m of dimension SpaceDim, 
    % which optimally represents the samples in a least-square sense. The affine
    % subspace is returned as a Point, and DirVecs, a SpaceDim x m matrix whose
    % rows are the spanning (or orthogonal - see below) directions.
    %
    % [...] = lsqAffineSpace(...,'orthogonal')
    % Fills the returned DirVecs with a basis to the *orthogonal* subspace to the aforementioned
    % optimal one. (i.e., the 'worst' directions with respect to the samples).
    %
    % [...] = lsqAffineSpace(...,'discardworst',k)
    % Iteratively discards the k samples farthest from the optimal affine subspace found.
    %
    % [...] = lsqAffineSpace(...,'discardthresh',th)
    % Iteratively discards samples whose distance from the found affine subspace exceeds
    % th .
    %
    % [... , samplesused] = lsqAffineSpace(...)
    % Returns a vector of indices of the samples that weren't discarded in the process.
    %
    
    
    discardbycount = false;
    discardnum = 0;
    discardbyval = false;
    discardthresh = 0;
    isorth = false;
    
    parseargs(varargin);
        
    samplesused = 1:size(V,1) ;
    
    varargout={};
    
    while true
	  n = numel(samplesused);
	  Y = V(samplesused,:);

	  Point = sum(Y) ./ n; % current center of mass
	  Y = Y - Point(ones(n,1),:);
	  
	  [eigVecs, eigVals] = eig(Y.' * Y);
	  eigVals=diag(eigVals);

	  % for a symmetric matrix, the lapack routines eig uses return values in
	  % ascending eigVals order, but we don't count on that here.
	  [dmp, prm] = sort(eigVals,'descend');
	  eigVecs = eigVecs(:,prm);

	  SpanVecs = eigVecs(:,1:SpaceDim) .' ;
	  OrthVecs = eigVecs(:,SpaceDim+1:end).';
	  
	  if discardbycount || discardbyval 
		% discard worst sample. Repeating calculation in this manner is NOT
		% equivalent to discarding in advance, say, the discardnum farthest vectors from 
		% the first computed affine space.

		OrthProjs = diag(Y * (OrthVecs.' * OrthVecs) * Y.' );  
		% may be suboptimal, especially for large datasets
		[maxval, maxidx] = max(OrthProjs);
		
		if discardbyval % discarding by distance thresh
		    
		    if maxval > discardthresh*discardthresh
			  samplesused(maxidx) = [];
		    else
			  break; % done discarding
		    end % if maxval > discardthresh*discardthresh
		    
		else % ignoring a specified numer of worst samples
		    discardnum = discardnum - 1;
		    if discardnum==-1
			  break;% done discarding
		    end
		    samplesused(maxidx) = [];
	    
		end % discardbyval
		
	  else % if discardbycount || discardbyval
		break;
	  end
	  
    end % while true
    
    if isorth
	  DirVecs = OrthVecs ;
    else
	  DirVecs = SpanVecs;
    end
    
    if nargout==3
	  varargout{1} = samplesused;
    end
    
    %******************************************************************    
    function parseargs(argcells)
	  j=1;
	  while j<=numel(argcells)
		switch lower(argcells{j})
		    case 'orthogonal'
			  isorth = true;
		    case 'discardworst'
			  discardbycount = true;
			  discardnum = argcells{j+1};
			  j=j+1;
		    case 'discardthresh'
			  discardbyval = true;
			  discardthresh = argcells{j+1};
			  if discardthresh<eps(max(V(:)))
				error('The specified discard thresh is too close to zero.');
			  end
			  j=j+1; 
		    otherwise
			  error('Unrecognized option: %s',num2str(argcells{j}) ); % legal also when argcells{j} is a string
		end
		j=j+1;
	  end %while
    end %parseargs

end % lsqAffineSpace