function [ loops ] = findBoundaries( pointset, iboundary, mode, Vu, Vv, parameter  )
%FINDBOUNDARIES Summary of this function goes here
%   Detailed explanation goes here

if strcmp(mode,'uv')
    found_loops = findBoundaries_uv(pointset, iboundary, Vu, Vv);
elseif strcmp(mode,'3D')
    if ~ exist('parameter','var')
        parameter = inf;
    end
    found_loops = findBoundaries_3D(pointset, iboundary, parameter );
else
    error('[findBoundaries] unknown mode\n');
end

% sort loops by size
for k = 1:numel(found_loops)
    loopsize(k) = numel(found_loops{k});
end
[sorted,idx] = sort(loopsize,'descend');

loops = [];
for k = 1:numel(idx)
    if sorted(k) > 10
        loops{k} = found_loops{idx(k)};
    else
        break;
    end
end

end













function [loops] = findBoundaries_3D(pointset, iboundary, threshfactor  )

dist_thresh = mean(nonzeros(pointset.e));

loops = [];
loopk = 1;

curset = iboundary;    
    
while ~isempty(curset) 
    
    vcur = curset(1);
    sorted = [vcur];
    curset = curset(curset~=vcur);

    done = 0;
    while ~isempty(curset) & ~done
        
        
        nbrs = curset;
        nbrdists2 = vmag2(vadd(pointset.v(nbrs,:),-pointset.v(vcur,:)));
        [dsort,idx] = sort(nbrdists2);
        
        nbridx = 1;
        vnbr = nbrs(idx(nbridx));

%         if numel(sorted) > 2
%             next_dist = vmag2(pointset.v(vcur,:) - pointset.v(vnbr,:));
%             vals = find(dsort < next_dist*5);
%             checknbrs = nbrs(idx(vals));
%             if numel(checknbrs) > 1
%                 vprev = normalize(pointset.v(vcur,:) - pointset.v(sorted(end-1),:));
%                 vnbrdot = vdot( normalize(vadd(pointset.v(checknbrs,:),-pointset.v(vcur,:))), vprev );
%                 [asort,aidx] = sort(vnbrdot,'descend');
%                 vnbr = checknbrs(aidx(1));
%             end
%         end
        

        next_dist = vmag(pointset.v(vcur,:) - pointset.v(vnbr,:));
        cur_loop_close_dist = vmag(pointset.v(vcur,:) - pointset.v(sorted(1),:));
        if ( numel(sorted) > 5 & cur_loop_close_dist < next_dist )
            fprintf('loop size: %d\n', numel(sorted));
            loops{loopk} = sorted';
            loopk = loopk+1;
            done = 1;
            continue;
        end
        
        if (next_dist > dist_thresh*threshfactor)
            fprintf('loop size: %d\n', numel(sorted));
            loops{loopk} = sorted';
            loopk = loopk+1;
            done = 1;
            continue;            
        end
        
        sorted = [sorted,vnbr];
        vcur = vnbr;
        curset = curset(curset~=vcur);       

        if ( numel(sorted) > 5 & vmag(pointset.v(vcur,:) - pointset.v(sorted(1),:)) < dist_thresh ) 
            fprintf('loop size: %d\n', numel(sorted));
            loops{loopk} = sorted';
            loopk = loopk+1;
            done = 1;
        end
        
    end
end

if ~isempty(sorted)
    fprintf('loop size: %d\n', numel(sorted));
    loops{loopk} = sorted';
    loopk = loopk+1;
end

end





function [loops] = findBoundaries_uv(pointset, iboundary, Vu, Vv  )

loops = [];
loopk = 1;

curset = iboundary;    
    
while ~isempty(curset) 
    
    vcur = curset(1);
    sorted = [vcur];
    curset = curset(curset~=vcur);

    done = 0;
    while ~isempty(curset) & ~done
        nbrdists2 = Vu(vcur,:).^2 + Vv(vcur,:).^2;
        nbrs = find(nbrdists2 > 0)';
        uv = [full(Vu(vcur,nbrs))', full(Vv(vcur,nbrs))'];
        nbrdists2 = vmag2(uv);
        [dsort,idx] = sort(nbrdists2);

        found = 0;
        for k = 1:numel(idx)
            vnbr = nbrs(idx(k));
            isb = find(curset==vnbr);
            if ~isempty(isb)
                found = vnbr;
                sorted = [sorted, vnbr];
                break;
            end
        end
        if ~found
            fprintf('loop size: %d\n', numel(sorted));
            loops{loopk} = sorted';
            loopk = loopk+1;
            done = 1;
        end

        vcur = vnbr;
        curset = curset(curset~=vcur);
    end

    if ~isempty(sorted)
        fprintf('loop size: %d\n', numel(sorted));
        loops{loopk} = sorted';
        loopk = loopk+1;
    end

end

% sort by loop size/length ?

end
