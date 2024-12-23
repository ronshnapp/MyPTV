clear all;
close all;
clc;

figure
target_file_all = load(['target_file']);
scatter3(target_file_all(:,1),target_file_all(:,2),target_file_all(:,3),'o','filled','MarkerEdgeColor','k','MarkerFaceColor','b'); hold on
daspect([1 1 1])

% floors = unique(target_file_all(:,1));
% for i = 1:length(floors)
%     N(i) = length(find(target_file_all(:,1)==floors(i)));
% end