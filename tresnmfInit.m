function [config, store] = tresnmfInit(config)                              
% tresnmfInit INITIALIZATION of the expLanes experiment trafficEstimationNMF
%    [config, store] = tresnmfInit(config)                                  
%      - config : expLanes configuration state                              
%      -- store  : processing data to be saved for the other steps          
                                                                            
% Copyright: <userName>                                                     
% Date: 04-Aug-2018                                                         
                                                                            
if nargin==0, trafficEstimationNMF(); return; else store=[];  end           
