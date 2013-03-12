% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

function plot12n(name,data)
    ld=length(data(1,:));
    if ld==2
        hf=figure;
        set(hf,'name',name,'NumberTitle','off');
%         if handles.islb
%             scatter(data(:,1),data(:,2),5,handles.labels);
%         else
            plot(data(:,1),data(:,2),'.r');
%         end
        title(name);
    else
        if ld==1
            hf=figure;
            set(hf,'name',name,'NumberTitle','off');
            %if handles.islb
                %scatter(data(:,1),zeros(length(data(:,1)),1),5,handles.labels);
            %else
                plot(data(:,1),zeros(length(data(:,1)),1),'.k');
            %end
            title(name);
        else
            %if handles.islb
               % scattern('Result of dimensionality reduction',data,handles.labels);
            %else
                plotn(name,data);
            %end
        end
        
    end