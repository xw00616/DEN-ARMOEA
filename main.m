
% The source code of AR-MOEA

% This code is extracted from PlatEMO, which can be downloaded at
% http://bimk.ahu.edu.cn/index.php?s=/Index/Software/index.html

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group
% You are free to use the PlatEMO for research purposes. All publications
% which use this platform or any code in the platform should acknowledge
% the use of "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and
% Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, IEEE Computational Intelligence Magazine, 2017, in press".
%--------------------------------------------------------------------------
    clear;clc;
    cd(fileparts(mfilename('fullpath')));%从当前文件的完整路径中提取出路径和文件名，然后切换当前工作目录到matlab
    addpath(genpath(cd));%将该路径下所有文件夹都添加到搜索路径中
    
    %Problems={'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7','ZDT'};
        Problems={'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9',...
                       'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7'};
    itr=20;D=20;M=3;FE=11*D-1+120;
    load ini20.mat;
    
    for Prob = 3:4%length(Problems)
        clear r;
        Problem=Problems{Prob};
        k = find(~isstrprop(Problem,'digit'),1,'last');
        a=[Problem(1) Problem(k+1) 'DARMOEA.mat'];
        for i=1:itr
            disp(sprintf('%u / %u loop begin', i, itr));
            data=ini(i).chrom;
            Global = GLOBAL('-algorithm',@ARMOEA,'-problem',str2func(Problem),'-D',D,'-M',M,...
                                          '-evaluation',FE,'-data',data);
            Global.Start();
            Population=Global.result{end,2};
            tr_x=Population.decs;tr_y=Population.objs;
            r(i).ch=[tr_x,tr_y];
            r(i).traintime=Global.traintime;
            r(i).runtime=Global.runtime;
            save(a,'r');
        end
    end
