%TESTING SCRIPT FOR FITNORMAL
%
%This script tests that the normal vector for the 
%best fit plane can be properly calculated by
%the fitNormal function.
%
%A test normal vector is chosen first and then points
%created that would like on its plane.  The points are
%passed to fitNormal and the resultant unit vector
%should match up with the initial known vector.  Special
%cases such as when the normal vector is collinear with
%the global coordinate axes are specifically tested.
%
%Author: Daniel Robert Couture
%e-mail address: daniel.robert.couture@gmail.com
%Release: 2
%Release date: Aug 11, 2012
function t_fitNormal()
	
	disp('Testing fitNormal')
	
	for i = 1:10
		pause(1)
		% Create a normal vector
		if i == 1
			disp('Testing yz planar points')
			test_normal = [1 0 0]';
		elseif i == 2
			disp('Testing xz planar points')
			test_normal = [0 1 0]';
		elseif i == 3
			disp('Testing xy planar points')
			test_normal = [0 0 1]';
		elseif i == 4
			disp('Testing yz planar points')
			test_normal = [-1 0 0]';
		elseif i == 5
			disp('Testing xz planar points')
			test_normal = [0 -1 0]';
		elseif i == 6
			disp('Testing xy planar points')
			test_normal = [0 0 -1]';
		else
			disp('Testing random skew planar points')
			test_normal = rand(3,1);
		end
		
		test_normal = test_normal / norm(test_normal);
		[x,y] = generate_unit_vectors(test_normal);
		
		% Generate random points in the plane
		data = [];
		for i = 1:50
			pt = x * (rand() * 2 - 1);
			pt = pt + y * (rand() * 2 - 1);
			data(i,:) = pt';
		end
		
		% Let's see what we get
		test_result = fitNormal(data, true);
		
		checkResults(test_normal,test_result,1e-10)
	end
	
	disp(sprintf('\nCheck for fit with noise'))
	for i = 1:6
		pause(1)
		% Create a normal vector
		if i == 1
			disp('Testing yz noisy planar points')
			test_normal = [1 0 0]';
		elseif i == 2
			disp('Testing xz noisy planar points')
			test_normal = [0 1 0]';
		elseif i == 3
			disp('Testing xy noisy planar points')
			test_normal = [0 0 1]';
		else
			disp('Testing random skew noisy planar points')
			test_normal = rand(3,1);
		end
		
		test_normal = test_normal / norm(test_normal);
		[x,y] = generate_unit_vectors(test_normal);
		
		% Generate random points in the plane
		data = [];
		for i = 1:50
			pt = x * (rand() * 2 - 1);
			pt = pt + y * (rand() * 2 - 1);
			pt = pt + test_normal * (rand() * 0.2 - 0.1);
			data(i,:) = pt';
		end
		
		% Let's see what we get
		test_result = fitNormal(data, true);
		
		checkResults(test_normal,test_result,0.001)
	end
	
end

function [x,y] = generate_unit_vectors(test_normal)
		% We need to create a local x and y vectors in the plane
		% Start with a random vector
		x = zeros(3,1);
		% While loop to keep going in case randomly generated
		% vectors are collinear to the test normal vector
		while sum(x) == 0
			x = cross(test_normal,rand(3,1));
		end
		% Normalize the local x vector
		x = x / norm(x);
		% Cross the test normal with the local x to generate a local y
		% and have 3 orthogonal vectors for testing.
		y = cross(test_normal,x);
		% Normalize the local y vector
		y = y / norm(y);
		
		% Verify that the three local vectors for testing are orthogonal.
		if abs(dot(test_normal,x)) > 1e-8
			error('Verify that the test normal is orthogonal to the plane local x axis');
		end
		if abs(dot(test_normal,y)) > 1e-8
			error('Verify that the test normal is orthogonal to the plane local y axis');
		end
		if abs(dot(x,y)) > 1e-8
			error('Verify that the plane local x and y axes are orthogonal');
		end
end

function checkResults(test_normal,test_result,allowed_error)
	% Verify that we get back a unit vector
	if norm(test_result) < 1 - 1e-8 && norm(test_result) > 1 + 1e-8
		error('fitNormal returned a unit vector')
	end
	
	% Verify that the result is nearly identical to the original test normal
	% vector +- some rounding error.
	if abs(dot(test_normal, test_result)) > (1 + allowed_error) || abs(dot(test_normal, test_result)) < (1 - allowed_error)
		error('fitNormal returned an appropriate vector')
	end
	
	disp(' *Passed*')
end
