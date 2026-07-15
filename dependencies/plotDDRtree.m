function plotDDRtree(projection, PTree, edges, label, label_level, FaceColor)
%% ========================================================================
% Plot a principal tree
%
%--------------------------------------------------------------------------
% Input
%   projection  : projection of data points in the sapce of reduced dimension
%   PTree       : data points of the principal tree
%   edges       : edges of the principal tree
%   label       : annotations of data points
%   label_level : annotation levels
%   FaceColor   : facecolor of data points
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================
if nargin<6
    % FaceColor=cbrewer('qual', 'Set1',9);
    FaceColor=defaultColor(length(label_level));
end
%% initializations
Y = label;
U = unique(Y(~isnan(Y)));
n_class = length(U);

%plot smaples
figure,
hold on

for i=1:n_class
    % scatter3(projection(1,Y==U(i)), projection(2,Y==U(i)), projection(3,Y==U(i)), 70, FaceColor(U(i),:), 'filled', 'MarkerFaceAlpha', 0.85,'MarkerEdgeColor','k');
    % scatter3(projection(1,Y==U(i)), projection(2,Y==U(i)), projection(3,Y==U(i)), 70, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.85,'MarkerEdgeColor','k');
    % scatter3(projection(1,Y==U(i)), projection(2,Y==U(i)), projection(3,Y==U(i)), 70, FaceColor(U(i),:), 'filled', 'MarkerFaceAlpha', 0.85);
    plot3(projection(1,Y==U(i)), projection(2,Y==U(i)), projection(3,Y==U(i)),...
        'o','MarkerFaceColor',FaceColor(U(i),:),'MarkerSize',6,'MarkerEdgeColor','k');
end
TreePoints = PTree(1:3,:);

% plot principal tree
[m,n] = size(edges);
for i=1:m
    for j=1:n
        if edges(i,j)~=0
            plot3([TreePoints(1,i), TreePoints(1,j)], [TreePoints(2,i), TreePoints(2,j)],...
                [TreePoints(3,i), TreePoints(3,j)],'-k','linewidth',4);
        end
    end
end
xlabel('DDR1')
ylabel('DDR2')
zlabel('DDR3')
if length(label_level)<20
    legend(label_level);
end
grid;
set(gca,'FontSize',14);
view(-122,14)
ax = gca; % Get current axes
ax.LineWidth = 1.5; % Set the axis line width to 2
grid on;
end