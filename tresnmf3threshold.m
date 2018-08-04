function [config, store, obs] = tresnmf3threshold(config, setting, data)                   
% tresnmf3threshold THRESHOLD step of the expLanes experiment trafficEstimationNMF         
%    [config, store, obs] = tresnmf3threshold(config, setting, data)                       
%      - config : expLanes configuration state                                             
%      - setting   : set of factors to be evaluated                                        
%      - data   : processing data stored during the previous step                          
%      -- store  : processing data to be saved for the other steps                         
%      -- obs    : observations to be saved for analysis                                   
                                                                                           
% Copyright: <userName>                                                                    
% Date: 04-Aug-2018                                                                        
                                                                                           
% Set behavior for debug mode                                                              
if nargin==0, trafficEstimationNMF('do', 3, 'mask', {2 1 1 1 1, ...
        2 7 3 1:2 1, 1 1 0 1 2, 1 1 1 1 1, 1:2 1 1 1 1}); return; else store=[]; obs=[]; end
                                                                                           
estimator =  data.estimator;
levels = cell(1,length(estimator));

if strcmp(setting.type,'nmf') && strcmp(setting.nmfType,'threshold')
    distanceMethod = setting.distanceMethod;
    displayDistance = setting.displayDistance;
    methodThreshold = setting.methodThreshold;
    
    valueThresholdCheck = 1;
    switch  methodThreshold
        case 'firm'
            threshold(1) = setting.thresholdFirmHigh;
            threshold(2) = setting.thresholdFirmLow;
            if threshold(1) < threshold(2)
                valueThresholdCheck = 0;
            end
        otherwise
            threshold = setting.threshold;
    end
    
    [soundMix] = cutSpectrogramEXP(estimator{1}.W0,setting);
    W0 = soundMix.W(1:soundMix.ind,:);
    [F,K] = size(W0);
    
    if strcmp(setting.displayDistance,'sigmoid')
        lambda = 2;
    end
    
    if valueThresholdCheck == 1
        for ii = 1:length(estimator)
            W = estimator{ii}.W(1:F,:);
            H = estimator{ii}.H;
            
            switch distanceMethod
                case 'cosine'
                    dist = sum(W.*W0)./(sqrt(sum(W.^2,1)).*sqrt(sum(W0.^2,1)));
                    dist(isnan(dist)) = 0;
                    [~, order] = sort(dist,'descend');
                    
                case 'none'
                    order = 1:K;
            end
            
            Wn = W(:,order);
            Hn = H(order,:);
            
            switch distanceMethod
                case 'none'
                    Wtraffic = Wn;
                    Htraffic = Hn;
                    
                otherwise
                    dist = dist(order);
                    
                    if strcmp(displayDistance,'sigmoid')
                        dist = 1./(1+exp(-lambda*dist));
                    end
                    vec = zeros(1,K);
                    
                    switch methodThreshold
                        case 'hard'
                            vec(dist>threshold) = 1;
                            vec(dist<=threshold) = 0;
                        case 'firm'
                            vec(dist>threshold(1)) = 1;
                            vec(dist<=threshold(1) & dist>threshold(2)) = 2;
                            vec(dist<=threshold(2)) = 3;
                            
                            mid = dist(vec==2);
                            if ~isempty(vec(vec==1))
                                vec(vec==1) = 1;
                            end
                            if ~isempty(mid)
                                vec(vec==2) = (mid-threshold(2))/(mid(1)-threshold(2));
                            end
                            if ~isempty(vec(vec==3))
                                vec(vec==3) = 0;
                            end
                    end
                    Wtraffic = Wn.*repmat(vec,F,1);
                    Htraffic = Hn;
            end
            
            [LpTrafficEstimate,LeqTrafficEstimate] = estimationLpEXP(Wtraffic*Htraffic,setting);
            
            levels{ii}.LpTrafficEstimate = LpTrafficEstimate;
            levels{ii}.LeqTrafficEstimate = LeqTrafficEstimate;
            levels{ii}.LpTraffic = estimator{ii}.LpTraffic;
            levels{ii}.LeqTraffic = estimator{ii}.LeqTraffic;
            levels{ii}.LpGlobal = estimator{ii}.LpGlobal;
            levels{ii}.LeqGlobal = estimator{ii}.LeqGlobal;
            levels{ii}.cost = estimator{ii}.cost;
        end
    else
        for ii = 1:length(estimator)
            levels{ii}.LpTrafficEstimate{1} = nan;
            levels{ii}.LeqTrafficEstimate(1) = nan;
            levels{ii}.LpTraffic = estimator{ii}.LpTraffic;
            levels{ii}.LeqTraffic = estimator{ii}.LeqTraffic;
            levels{ii}.LpGlobal = estimator{ii}.LpGlobal;
            levels{ii}.LeqGlobal = estimator{ii}.LeqGlobal;
            levels{ii}.cost = estimator{ii}.cost;
        end
        
    end
    store.levels = levels;
else
    for ii = 1:length(estimator)
        levels{ii}.LeqTrafficEstimate = estimator{ii}.LeqTrafficEstimate;
        levels{ii}.LpTrafficEstimate = estimator{ii}.LpTrafficEstimate;
        levels{ii}.LpTraffic = estimator{ii}.LpTraffic;
        levels{ii}.LeqTraffic = estimator{ii}.LeqTraffic;
        levels{ii}.LpGlobal = estimator{ii}.LpGlobal;
        levels{ii}.LeqGlobal = estimator{ii}.LeqGlobal;
        levels{ii}.cost = estimator{ii}.cost;
    end
    store.levels = levels;
end