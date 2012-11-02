function [subf] = get_subfeature(f,type)
% type
% 1 Curvature
% 2 AGD
% 3 SDF
% 4 SC
% 5 SpinImage

% features=features(395:395+18-1,:); % AGD
% features=features(395,:); % AGD
% features=features(413:413+72-1,:); % SDF
% features=features(413,:); % SDF
% features=features(1:64-1,:); % Curvature
% features=features(65:65+60-1,:); % PCA
% features=features(7,:); % 平均曲率
% features=features(8,:); % scale1高斯曲率
% features=features(24,:); % scale2 高斯曲率
% features=features(40,:); % scale3 高斯曲率
% features=features(56,:); % scale4 高斯曲率
% features=features(85:85+270-1,:); % SC   %%%TODO::::85是哪里来的。。。 应该是125
% features=features(85:85+30-1,:); % scale1 SC
% features=features(509:509+100-1,:); % SpinImage

switch type
case 1
    subf=f(8,:); % scale1高斯曲率
case 2
    subf=f(395,:);
case 3
    subf=f(413,:);
case 4
    subf=f(125:125+270-1,:);
case 5
    subf=f(509:509+100-1,:);
otherwise
    subf=0;
end
end
