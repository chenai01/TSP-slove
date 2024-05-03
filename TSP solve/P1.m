
% 示例距离矩阵
% distances = [
%     0, 10, 15, 20;
%     10, 0, 35, 25;
%     15, 35, 0, 30;
%     20, 25, 30, 0
% ];
% distances= [
%     0   5   8   4   7   6   3   2   9   1;
%     5   0   6   2   8   9   4   3   7   5;
%     8   6   0   3   9   5   7   6   4   8;
%     4   2   3   0   6   4   5   2   8   3;
%     7   8   9   6   0   3   6   5   2   7;
%     6   9   5   4   3   0   8   7   6   4;
%     3   4   7   5   6   8   0   2   3   5;
%     2   3   6   2   5   7   2   0   4   6;
%     9   7   4   8   2   6   3   4   0   9;
%     1   5   8   3   7   4   5   6   9   0
%     ];
% distances=[ 0,	2,  5,  7;
%             2,	0,	8,	3;
%             5,	8,	0,	1;
%             7,	3,	1,	0];

% distances= [
%     0   3   7   2   5   6   1   4;
%     3   0   4   5   6   7   2   3;
%     7   4   0   6   2   3   5   8;
%     2   5   6   0   4   3   2   5;
%     5   6   2   4   0   1   3   6;
%     6   7   3   3   1   0   4   7;
%     1   2   5   2   3   4   0   3;
%     4   3   8   5   6   7   3   0
%     ];
% distances=[0,10,inf,inf,6,8;
%          10,0,13,inf,inf,12;
%          inf,13,0,23,inf,20;
%          inf,inf,23,0,7,6;
%          6,inf,inf,7,0,9;
%          8,12,20,6,9,0];
for i= 1:13
    distances=floor(10*rand(i,i));

    tic
    [shortest_path, shortest_distance] = tsp_bruteforce(distances);
    shortest_path=[shortest_path,1];
    % fprintf('最短路径: %s\n', num2str(shortest_path));
    % fprintf('最短距离: %d\n', shortest_distance);
    elapsed_time = toc;
    disp([num2str(i)   'Elapsed time: ' num2str(elapsed_time) ' seconds']);
end
function [shortest_path, shortest_distance] = tsp_bruteforce(distances)
n = size(distances, 1);
cities = 2:n;
shortest_path = [];
shortest_distance = Inf;

% 生成所有可能的路径
all_permutations = perms(cities);
temp=ones(size(all_permutations,1),1);
all_permutations=[temp,all_permutations];
% 初始化保存每种情况路径和距离的矩阵
all_paths = zeros(size(all_permutations));
all_distances = zeros(size(all_permutations, 1), 1);

% 计算每条路径的总长度
for i = 1:size(all_permutations, 1)
    total_distance = 0;
    path = all_permutations(i,:);
    for j = 1:(n - 1)
        total_distance = total_distance + distances(path(j), path(j + 1));
    end
    % 回到起点城市
    total_distance = total_distance + distances(path(n), path(1));

    % 更新最短路径和距离
    all_paths(i, :) = path;
    all_distances(i) = total_distance;
    if total_distance < shortest_distance
        shortest_distance = total_distance;
        shortest_path = path;
    end
end

    % % 显示每种情况的路径和距离
    % fprintf('所有可能的路径和距离：\n');
    % for i = 1:size(all_paths, 1)
    %     fprintf('路径：%s，距离：%d\n', num2str(all_paths(i, :)), all_distances(i));
    % end

end
