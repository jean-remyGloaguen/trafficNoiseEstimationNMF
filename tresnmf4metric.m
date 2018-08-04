function [config, store, obs] = tresnmf4metric(config, setting, data)                      
% tresnmf4metric METRIC step of the expLanes experiment trafficEstimationNMF               
%    [config, store, obs] = tresnmf4metric(config, setting, data)                          
%      - config : expLanes configuration state                                             
%      - setting   : set of factors to be evaluated                                        
%      - data   : processing data stored during the previous step                          
%      -- store  : processing data to be saved for the other steps                         
%      -- obs    : observations to be saved for analysis                                   
                                                                                           
% Copyright: <userName>                                                                    
% Date: 04-Aug-2018                                                                        
                                                                                           
% Set behavior for debug mode                                                              
if nargin==0, trafficEstimationNMF('do', 4, 'mask', {2 1 1 1 1, ...
        2 7 3 1:2 1, 1 1 0 1 2, 1 1 1 1 1, 1:2 1 1 1 1}); return; else store=[]; obs=[]; end
                                                                                           
levels = data.levels;
pref = setting.p0;

LeqGlobal = zeros(length(levels),1);
LeqTrafficEstimate = LeqGlobal;
LeqTrafficEstimatedB = LeqGlobal;
LeqTrafficExact = LeqGlobal;
LeqTrafficExactdB = LeqGlobal;
LpTrafficExactTot = [];
LpTrafficEstTot = [];

%% EXTRACTION LEVELS
for ii = 1:length(levels)
    LpTrafficExact = levels{ii}.LpTraffic;
    LpTrafficExact(LpTrafficExact==0) = pref;
    
    LpTrafficEstimate = levels{ii}.LpTrafficEstimate{1};
    LpTrafficEstimate(LpTrafficEstimate==0) = pref;
    
    LpGlobal = levels{ii}.LpGlobal;
    LpGlobal(LpGlobal==0) = pref;
       
    LeqTrafficEstimate(ii) = levels{ii}.LeqTrafficEstimate(1);
    LeqTrafficExact(ii) = levels{ii}.LeqTraffic;
    LeqGlobal(ii) = levels{ii}.LeqGlobal;
    
    if size(LpGlobal,2)>size(LpTrafficEstimate,2)
        LpGlobal = LpGlobal(:,1:size(LpTrafficEstimate,2));
    elseif size(LpGlobal,2)<size(LpTrafficEstimate,2)
        LpTrafficEstimate = LpTrafficEstimate(:,1:size(LpGlobal,2));
    end
    
    
    %% sound level in dB
    LpTrafficExactdB = 20*log10(LpTrafficExact/pref);
    LpTrafficEstimatedB = 20*log10(LpTrafficEstimate/pref);
    LeqTrafficEstimatedB(ii) = 20*log10(LeqTrafficEstimate(ii)/pref);
    LeqTrafficExactdB(ii) = 20*log10(LeqTrafficExact(ii)/pref);
    LeqGlobaldB(ii) = 20*log10(LeqGlobal(ii)/pref);
    
    LpTrafficExactTot = [LpTrafficExactTot LpTrafficExactdB];
    LpTrafficEstTot = [LpTrafficEstTot LpTrafficEstimatedB];
   
    obs.cost(ii) = levels{ii}.cost(end);

end
mae = sum(abs(LeqTrafficEstimatedB-LeqTrafficExactdB))/length(LeqTrafficEstimatedB);
for ii = 1:length(LpTrafficEstTot)/60
    LeqTrafficExactMin(ii) = 10*log10(1/60*sum(10.^(LpTrafficExactTot((ii-1)*60+1:ii*60)./10))); 
    LeqTrafficEstMin(ii) = 10*log10(1/60*sum(10.^(LpTrafficEstTot((ii-1)*60+1:ii*60)./10)));
end
mae60 = sum(abs(LeqTrafficEstMin-LeqTrafficExactMin))/length(LeqTrafficEstMin);

ec = std(abs(LeqTrafficEstimatedB-LeqTrafficExactdB));

obs.LeqGlobal = 20*log10(LeqGlobal./pref);
obs.LeqTrafficExa = LeqTrafficExactdB;
obs.LeqTrafficEst = LeqTrafficEstimatedB;
obs.mae = mae;
obs.mae60 = mae60;
obs.ec = ec;