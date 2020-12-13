function Distance = CalDistance(PopObj,RefPoint)
% Calculate the distance between each solution to each adjusted reference
% point

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group
% You are free to use the PlatEMO for research purposes. All publications
% which use this platform or any code in the platform should acknowledge
% the use of "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and
% Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, IEEE Computational Intelligence Magazine, 2017, in press".
%--------------------------------------------------------------------------

    %Algorithm4
    N  = size(PopObj,1);
    NR = size(RefPoint,1);

    %% Adjust the location of each reference point
    %if unique(PopObj)==1
        %PopObj=PopObj+(2*rand(size(PopObj,1),size(PopObj,2))-1)*1e-06;
    %end
    index=find(sum(PopObj.^2,2)==0);
    if ~isempty(index)
        PopObj(index,:)=PopObj(index,:)+1e-06*rand(size(PopObj(index,:)));
    end
    Cosine = 1 - pdist2(PopObj,RefPoint,'cosine');%种群中个体与参考点的余弦值 N*NR
    NormR  = sqrt(sum(RefPoint.^2,2)); %NR*1
    NormP  = sqrt(sum(PopObj.^2,2)); %N*1
    d1     = repmat(NormP,1,NR).*Cosine; %N*NR
    d2     = repmat(NormP,1,NR).*sqrt(1-Cosine.^2); %N*NR
    [~,nearest] = min(d2,[],1); %1*NR
    RefPoint    = RefPoint.*repmat(d1(N.*(0:NR-1)+nearest)'./NormR,1,size(RefPoint,2));
    
    %% Calculate the distance between each solution to each point
    Distance = pdist2(PopObj,RefPoint);
end