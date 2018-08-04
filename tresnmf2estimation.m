function [config, store, obs] = tresnmf2estimation(config, setting, data)                  
% tresnmf2estimation ESTIMATION step of the expLanes experiment trafficEstimationNMF       
%    [config, store, obs] = tresnmf2estimation(config, setting, data)                      
%      - config : expLanes configuration state                                             
%      - setting   : set of factors to be evaluated                                        
%      - data   : processing data stored during the previous step                          
%      -- store  : processing data to be saved for the other steps                         
%      -- obs    : observations to be saved for analysis                                   
                                                                                           
% Copyright: <userName>                                                                    
% Date: 04-Aug-2018                                                                        
                                                                                           
% Set behavior for debug mode                                                              
if nargin==0, trafficEstimationNMF('do', 2, 'mask', {2 1 1 1 1, ...
        2 7 3 1:2 1, 1 1 0 1 2, 1 1 1 1}); return; else store=[]; obs=[]; end
                                                                                           
dataset = setting.dataset;
type = setting.type;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DICTIONARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(type,'nmf')
    dictionary.W = data.dictionary.W;
    dictionary.frequency = data.dictionary.frequency;
    dictionary.indTraffic = data.dictionary.indTraffic;
else
    dictionary = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BASELINE GLOBALE ERROR + FILTRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch dataset
    case 'ambience'
        creationSceneDir = strcat(config.inputPath, dataset, filesep, setting.aType, filesep);
        
    case 'SOUR'
        creationSceneDir = strcat(config.inputPath, dataset, filesep, setting.sType, filesep);

end

files = dir(strcat(creationSceneDir,filesep,'*.wav'));
globalName = cell(1,length(files)/2);
ind = 1;

for ii = 1:length(files)
    if ~isempty(strfind(files(ii).name,'traffic'))
        globalName{ind} = files(ii).name(1:end-12);
        ind = ind+1;
    end
end

numberScene = length(globalName)-floor((setting.sceneSelect/100*length(globalName)));
switch dataset
    case 'ambience' 
        TIR = setting.TIR;
        
    case 'SOUR'
        LeqRef = 2e-5*10.^([69 70 73 76]./20);
        amb = {'park','quietStreet','noisyStreet','veryNoisyStreet'};
        ind = contains(amb,setting.sType)==1;
end

estimator = cell(numberScene,1);

if strcmp(setting.nmfType,'threshold')
    estimator{1}.W0 = dictionary.W;
end

if isempty(config.sequentialData)
    sequentialData = cell(numberScene,1);
else
    sequentialData = config.sequentialData;
end


for ii = 1:numberScene
    
    fileTraffic = audioread(strcat(creationSceneDir,globalName{ii},'_traffic.wav'));
    fileRest = audioread(strcat(creationSceneDir,globalName{ii},'_interfering.wav'));
    
    if strcmp(dataset,'ambience')        
        %% modif SNR
        A = rms(fileTraffic);
        B = rms(fileRest);
        if B ~= 0
            SNR_temp = 20*log10(A/B);
            facteur = 10.^((TIR-SNR_temp)/20);
            fileTraffic = facteur.*(fileTraffic);  % fichier perturbateur modifi?*
        end 
    end

    fileTraffic(fileTraffic == 0) = eps;
    fileTot = fileRest+fileTraffic;
    
    if strcmp(dataset,'SOUR')
        LeqGlobal = rms(fileTot);
        weightLevel = LeqRef(ind)/LeqGlobal;
        fileTot = weightLevel*fileTot;
        fileTraffic = weightLevel*fileTraffic;
    end
     
    
    [Vtraffic] = audio2SpectrogramEXP(fileTraffic',setting);
    [Lp,Leq] = estimationLpEXP(Vtraffic,setting);
    estimator{ii}.LpTraffic =  Lp{1};
    estimator{ii}.LeqTraffic = Leq(1);
    
    [V,Vlinear] = audio2SpectrogramEXP(fileTot',setting);
    [Lp,Leq] = estimationLpEXP(Vlinear,setting);
    estimator{ii}.LpGlobal = Lp{1};
    estimator{ii}.LeqGlobal = Leq(1);
    
    switch type
        case 'filter'
            [LeqFiltre,LpFiltre] = lowPassFilter(Vlinear,setting);
            estimator{ii}.LeqTrafficEstimate = LeqFiltre;
            estimator{ii}.LpTrafficEstimate = LpFiltre;
            sequentialData{ii} = 0;
            estimator{ii}.cost = 0;
            
        case 'nmf'
            [NMF,sequentialData{ii}] =...
                NMFestimationEXP(V,dictionary,setting,sequentialData{ii});
            
            estimator{ii}.H = NMF.H;
            estimator{ii}.cost = NMF.cost(end);
            
            if strcmp(setting.nmfType,'threshold')
                estimator{ii}.W = NMF.W;
            else
                estimator{ii}.LeqTrafficEstimate = NMF.LeqTrafficEstimate;
                estimator{ii}.LpTrafficEstimate = NMF.LpTrafficEstimate;
            end
    end 
end

config.sequentialData = sequentialData;
store.estimator = estimator;
