% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

function scatter12n(name,data,labels)

ld=length(data(1,:));
    if ld==2
        hf=figure;
        set(hf,'name',name,'NumberTitle','off');
        %if handles.islb
            scatter(data(:,1),data(:,2),5,labels);
%         else
%             scatter(data(:,1),data(:,2),5);
%         end
        title(name);
    else
        if ld==1
            hf=figure;
            set(hf,'name',name,'NumberTitle','off');
            %if handles.islb
                scatter(data(:,1),zeros(length(data(:,1)),1),5,labels);
%             else
%                 scatter(data(:,1),zeros(length(data(:,1)),1),5);
%             end
            title(name);
        else
            %if handles.islb
                scattern(name,data,labels);
%             else
%                 scattern('Result of dimensionality reduction',data);
%             end
        end
        
    end