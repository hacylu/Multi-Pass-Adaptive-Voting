%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for showing the centriod of nuclei in the mask on the image
% LshowCrossfromBWonIM(mask,IM,figNo)
% Input:
%   -im     given image
%   -mask   mask indicates the ROI
%
% (c) Edited by Cheng Lu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  3rd, Aug, 2011
% If you have any problem feel free to contact me.
% Please address questions or comments to: hacylu@yahoo.com

% Terms of use: You are free to copy,
% distribute, display, and use this work, under the following
% conditions. (1) You must give the original authors credit. (2) You may
% not use or redistribute this work for commercial purposes. (3) You may
% not alter, transform, or build upon this work. (4) For any reuse or
% distribution, you must make clear to others the license terms of this
% work. (5) Any of these conditions can be waived if you get permission
% from the authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LshowCrossfromBWonIM(mask,IM,figNo,str_title)
% subplot('position',[0 0 1 1]); % Use the whole window
if nargin<3
    show(IM); hold on;
    str_title='';
else if nargin<4
        %     scrsz = get(0,'ScreenSize');
        figure(figNo);
        %     set(gcf, 'Position',[ 1 1 scrsz(3)/2 scrsz(4)]);
        %     figure(figNo);
        %     imshow(IM,[],'InitialMagnification',100);
        imshow(IM);
        hold on;
        str_title='';
    else
        figure(figNo);
        imshow(IM);
        hold on;
    end
end
%% get contour of the mask
ss=regionprops(mask,'Centroid');
for i=1:length(ss)
    curB=ss(i).Centroid;
    plot(curB(1),curB(2),'+b','LineWidth',3);
end
axis('image'); axis('off');
title(str_title);
hold off;
end