function MatingPool = MatingSelection(Obj,RefPoint,Range)
% The mating selection of AR-MOEA

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group
% You are free to use the PlatEMO for research purposes. All publications
% which use this platform or any code in the platform should acknowledge
% the use of "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and
% Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, IEEE Computational Intelligence Magazine, 2017, in press".
%--------------------------------------------------------------------------

    %% Calculate the degree of violation of each solution

    %% Calculate the fitness of each feasible solution based on IGD-NS
        % Calculate the distance between each solution and point
        N = size(Obj,1);
        Distance    = CalDistance(Obj-repmat(Range(1,:),N,1),RefPoint);%N*NR
        Convergence = min(Distance,[],2);%每个solution与参考点的最小夹角
        [dis,rank]  = sort(Distance,1);%按列排列, 对参考点而言
        % Calculate the fitness of noncontributing solutions
        Noncontributing = true(1,N);
        Noncontributing(rank(1,:)) = false;
        METRIC   = sum(dis(1,:)) + sum(Convergence(Noncontributing));
        %所有参考点的最小夹角之和+所有非贡献解的最小夹角之和
        fitness  = inf(1,N);
        fitness(Noncontributing) = METRIC - Convergence(Noncontributing);
        % Calculate the fitness of contributing solutions
        for p = find(~Noncontributing)
            temp = rank(1,:) == p;
            noncontributing = false(1,N);
            noncontributing(rank(2,temp)) = true;%除去该贡献点后，新的贡献点index
            noncontributing = noncontributing & Noncontributing;%新的贡献点是否是原来的非贡献点
            fitness(p) = METRIC - sum(dis(1,temp)) + sum(dis(2,temp)) - sum(Convergence(noncontributing));
            %除去该贡献点后，METRIC调整
        end


    %% Combine the fitness of feasible solutions with the fitness of infeasible solutions
    Fitness = fitness;
    
    %% Binary tournament selection
    CV=zeros(N,1);
    MatingPool = TournamentSelection(2,ceil(N/2)*2,CV,-Fitness);
end