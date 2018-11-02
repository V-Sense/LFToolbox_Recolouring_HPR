function [axis_limits] = DisplayPoints3D(Model, Scene, sampling, axis_limits)
%%=====================================================================
%% $RCSfile: DisplayPoints3D.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 21:34:29 +0000 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================

%set(gca,'FontSize',16,'FontName','Times','FontWeight','bold');

screensize = get(0,'ScreenSize');
sz = [576, 1024];
figure(1)%[ ceil((screensize(3)-sz(2))/2), ceil((screensize(4)-sz(1))/2), sz(2), sz(1)]);
subplot(1,2,1);%1,[0.01  0.4850 0.3200 .47]); 

   
scatter3(Model(:,1),Model(:,2),Model(:,3),5, lab2rgb(Model), 'filled');%'r+', 'MarkerSize', 8, 'LineWidth',1.5);
axis equal;
view(88,11);
xlim([0, 255]);
ylim([-100, 100]);   
zlim([-100, 100]);

subplot(1,2,2);%1,[0.3400  0.4850 0.3200 .47]); 
scatter3(Scene(:,1),Scene(:,2),Scene(:,3),5, lab2rgb(Scene), 'filled');
%plot3(Scene(:,1),Scene(:,2),Scene(:,3),'bo', 'MarkerSize', 8, 'LineWidth',1.5);
axis equal;
view(88,11);
xlim([0, 255]);
ylim([-100, 100]);   
zlim([-100, 100]);

% if (nargin<3)
% %    axis_limits = determine_border(Model, Scene);
%     sampling = 0;
% end
% 
% m = size(Model,1);
% if (sampling>0)
%     for i=1:sampling:m
%         text(Model(i,1), Model(i,2), Model(i,3), [' \leftarrow',sprintf('%d',i)]);
%     end
% end
% 
% m = size(Scene,1);
% if (sampling>0)
%     for i=1:sampling:m
%         text(Scene(i,1), Scene(i,2), Scene(i,3), [' \leftarrow',sprintf('%d',i)]);
%     end
% end
% 
% if (nargin<4)
%     axis_limits = determine_border(Model, Scene);
% end 

%pbaspect([1,1,1]);