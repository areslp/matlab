function n = fitNormal(data, show_graph)
%FITNORMAL - Fit a plane to the set of coordinates
%
%For a passed list of points in (x,y,z) cartesian coordinates,
%find the plane that best fits the data, the unit vector
%normal to that plane with an initial point at the average
%of the x, y, and z values.
%
% Input:
% @data
%   value - Matrix composed of of N sets of (x,y,z) coordinates
%           with dimensions Nx3
%   type  - Nx3 matrix
% @show_graph
%   value   - Option to display plot the result
%   default - false
%   type    - logical
%
% Return:
% @n
%   value - Unit vector that is normal to the fit
%           plane and initial point the average of the
%           @data.y, @data.y, @data.z
%   type  - 3x1 vector
	
	if nargin == 1
		show_graph = false;
	end
	
	for i = 1:3
		X = data;
		X(:,i) = 1;
		
		X_m = X' * X;
		if det(X_m) == 0
			can_solve(i) = 0;
			continue
		end
		can_solve(i) = 1;
		
		% Construct and normalize the normal vector
		coeff = (X_m)^-1 * X' * data(:,i);
		c_neg = -coeff;
		c_neg(i) = 1;
		coeff(i) = 1;
		n(:,i) = c_neg / norm(coeff);
		
	end
	
	if sum(can_solve) == 0
		error('Planar fit to the data caused a singular matrix.')
		return
	end
	
	% Calculating residuals for each fit
	center = mean(data);
	off_center = [data(:,1)-center(1) data(:,2)-center(2) data(:,3)-center(3)];
	for i = 1:3
		if can_solve(i) == 0
			residual_sum(i) = NaN;
			continue
		end
		
		residuals = off_center * n(:,i);
		residual_sum(i) = sum(residuals .* residuals);
		
	end
	
	% Find the lowest residual index
	best_fit = find(residual_sum == min(residual_sum));
	
	% Possible that equal mins so just use the first index found
	n = n(:,best_fit(1));
	
	if ~show_graph
		return
	end
	
	L=plot3(data(:,1),data(:,2),data(:,3),'ro','Markerfacecolor','r'); % Plot the original data points
	hold on;
	set(get(L, 'Parent'),'DataAspectRatio',[1 1 1],'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);
	
	norm_data = [mean(data); mean(data) + n'];
	
	% Plot the original data points
	L=plot3(norm_data(:,1),norm_data(:,2),norm_data(:,3),'b-','LineWidth',3);
	set(get(get(L,'parent'),'XLabel'),'String','x','FontSize',14,'FontWeight','bold')
	set(get(get(L,'parent'),'YLabel'),'String','y','FontSize',14,'FontWeight','bold')
	set(get(get(L,'parent'),'ZLabel'),'String','z','FontSize',14,'FontWeight','bold')
	title(sprintf('Normal Vector: <%0.3f, %0.3f, %0.3f>',n),'FontWeight','bold','FontSize',14)
	grid on;
	axis square;
	hold off;
end