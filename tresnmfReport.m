function config = tresnmfReport(config)                                  
% tresnmfReport REPORTING of the expLanes experiment trafficEstimationNMF
%    config = tresnmfInitReport(config)                                  
%       config : expLanes configuration state                            
                                                                         
% Copyright: <userName>                                                  
% Date: 04-Aug-2018                                                      
                                                                         
if nargin==0, trafficEstimationNMF('report', 'rhv'); return; end         
                                                                         
config = expExpose(config, 't');                                         
