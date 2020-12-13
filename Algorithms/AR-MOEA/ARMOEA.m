function ARMOEA(Global)
% <algorithm> <A-G>
% An Indicator Based Multi-Objective Evolutionary Algorithm with Reference
% Point Adaptation for Better Versatility

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group
% You are free to use the PlatEMO for research purposes. All publications
% which use this platform or any code in the platform should acknowledge
% the use of "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and
% Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, IEEE Computational Intelligence Magazine, 2017, in press".
%--------------------------------------------------------------------------
    %% Parameter setting
    W          = UniformPoint(Global.N,Global.M);
    NI    = 11*Global.D-1;
    Data=Global.data;
    P    = INDIVIDUAL(repmat(Global.upper-Global.lower,NI,1).*Data+repmat(Global.lower,NI,1));
    tr_x=P.decs;tr_y=P.objs;
    [tr_xx,ps]=mapminmax(tr_x');tr_xx=tr_xx';
    [tr_yy,qs]=mapminmax(tr_y');tr_yy=tr_yy';
    Params.ps=ps;Params.qs=qs;
    timestore=[];Ke=3;RatioOld=[];delta=0.05;wmax=20;
    
    while Global.NotTermination(P, timestore)
        %% Model
        t1=clock;
        if isempty(timestore)
            [net, Params]=trainmodel(tr_xx, tr_yy, Params);
        else
            net=updatemodel(tr_xx, tr_yy, Params, net);
        end
        timestore=[timestore;etime(clock, t1)];
           
        %% Generate the sampling points and random population
        popsize=50;
        PopDec=repmat(Global.upper-Global.lower,popsize,1).*rand(popsize,Global.D)+repmat(Global.lower,popsize,1);
        [PopObj, PopMSE]=Estimate(PopDec, net, Params, Global.M);
        [Archive,RefPoint,Range, Ratio] = UpdateRefPoint(PopObj,W,[]);
        if isempty(RatioOld)
            RatioOld=Ratio;
        end
        
        w=1;
        %% Start the interations
        while w < wmax
            MatingPool = MatingSelection(PopObj,RefPoint,Range);%PopulationµÄindex
            OffspringDec  = Global.VariationDec(PopDec(MatingPool,:));
            [OffspringObj, OffspringMSE] = Estimate(OffspringDec, net, Params, Global.M);
            [Archive,RefPoint,Range, Ratio] = UpdateRefPoint([Archive;OffspringObj],W,Range);
            MediatePopDec=[PopDec;OffspringDec];
            MediatePopObj=[PopObj;OffspringObj];
            MediatePopMSE=[PopMSE;OffspringMSE];
            [Index,Range]       = EnvironmentalSelection(MediatePopObj,RefPoint,Range,popsize);
            PopDec=MediatePopDec(Index,:);
            PopObj=MediatePopObj(Index,:);
            PopMSE=MediatePopMSE(Index,:);
            w=w+1;
        end
        flag=RatioOld-Ratio< delta;
        PopNew=IndividualSelect(PopDec, PopObj, PopMSE, Ke, flag);
        RatioOld=Ratio;
        New       = INDIVIDUAL(PopNew);n=length(New);
        P        = [P,New];
        [tr_x, tr_y]=SelectTrainData(P, NI, n);
        [tr_xx,ps]=mapminmax(tr_x');tr_xx=tr_xx';
        [tr_yy,qs]=mapminmax(tr_y');tr_yy=tr_yy';
        Params.ps=ps;Params.qs=qs;     
    end