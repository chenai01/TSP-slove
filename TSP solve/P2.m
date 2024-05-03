
% 示例用法
% distances = [0, 10, 15, 20;
%              10, 0, 35, 25;
%              15, 35, 0, 30;
%              20, 25, 30, 0];
% 
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

for i= 22:22
    distances=floor(10*rand(i,i));

    tic
    [tour, minCost,dp] = tsp_dp(distances);
    fprintf('Optimal tour: %s\n', num2str(tour));
    fprintf('Minimum cost: %d\n', minCost);
    elapsed_time = toc;
    disp([num2str(i)   'Elapsed time: ' num2str(elapsed_time) ' seconds']);
end
function [optTour, minCost,dp] = tsp_dp(distances)
    % Number of cities
    n = size(distances, 1);
    
    % Set up dynamic programming table
    dp = inf(n, 2^n);
    dp(1, 1) = 0;
    
    % Initialize parent array to store optimal path
    parent = zeros(n, 2^n);
    
    % Fill up dynamic programming table
    for set = 1:(2^n)-1
        for current = 1:n
            if bitget(set, current) == 0
                continue;
            end
            
            prev_set = set - 2^(current-1);
            if prev_set == 0
                prev_set = 1;
            end
            
            for prev = 1:n
                if bitget(prev_set, prev) == 0
                    continue;
                end
                
                if dp(current, set) > dp(prev, prev_set) + distances(prev, current)
                    dp(current, set) = dp(prev, prev_set) + distances(prev, current);
                    parent(current, set) = prev;
                end
            end
        end
    end
    
    % Find the optimal tour
    minCost = inf;
    for i = 2:n
        if minCost > dp(i, 2^n-1) + distances(i, 1)
            minCost = dp(i, 2^n-1) + distances(i, 1);
            last_city = i;
        end
    end
    
    % Reconstruct the optimal tour
    optTour = zeros(n, 1);
    set = 2^n-1;
    optTour(n) = 1;
    for i = n-1:-1:1
        optTour(i) = parent(last_city, set);
        set = set - 2^(last_city-1);
        last_city = optTour(i);
    end

end
