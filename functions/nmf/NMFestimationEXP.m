function [NMF,sequentialData] = NMFestimationEXP(V,dictionary,setting,sequentialData)

[soundMix] = cutSpectrogramEXP(dictionary.W,setting);
W = soundMix.W(1:soundMix.ind,:);
V = V{1}(1:soundMix.ind,:);
binTemp = size(V,2);

switch setting.nmfType
    case 'supervised'
        
        if isempty(sequentialData)
            iteration = setting.iteration;
            rng(soundMix.seed)
            H = rand(size(W,2),binTemp);
            
        else
            % continuing step of the sequential run
            iteration = setting.iteration-sequentialData.numberIteration;
            H = sequentialData.H;
            
        end

        NMF = algo_nmfSupervisedEXP(H,W,V,iteration,setting);        
        NMF.cost = betadivEXP(V,NMF.Vap,setting.beta,NMF.H,setting.sparsity,setting.smoothness);
        
        [NMF] = separationTrafficEXP(NMF,soundMix,dictionary,setting);
        

    case 'semi-supervised'
        
        Wrand = setting.SS_sizeWrand;
        if isempty(sequentialData)
            iteration = setting.iteration;
            rng(soundMix.seed); Y =  rand(size(W,1),Wrand);
            rng(soundMix.seed); H = rand(size(W,2),binTemp);
            rng(soundMix.seed); Z = rand(Wrand,binTemp);
            
        else
            % continuing step of the sequential run
            iteration = setting.iteration-sequentialData.numberIteration;
            H = sequentialData.H;
            Z = sequentialData.Z;
            Y = sequentialData.Y;
            
        end
        NMF = algo_nmfSemiSupervisedEXP(H,W,Y,Z,V,iteration,setting);
        NMF.cost = betadivEXP(V,NMF.Vtot,setting.beta,NMF.H,setting.sparsity,setting.smoothness);
        
        sequentialData.Y = NMF.Y;
        sequentialData.Z = NMF.Z;
        
        [NMF] = separationTrafficEXP(NMF,soundMix,dictionary,setting);

    case 'threshold'
        W0 = W(1:soundMix.ind,:);      % initial dictionary
        
        if isempty(sequentialData)
            iteration = setting.iteration;
            W = W0;                             % first iteration W = W0
            rng(soundMix.seed)
            H = rand(size(W,2),binTemp);    % H initiate randomly
        else
            % continuing step of the sequential run
            iteration = setting.iteration-sequentialData.numberIteration;
            H = sequentialData.H;
            W = sequentialData.W;
        end
        %
        NMF = algo_nmfUnsupervisedEXP(H,W,V,iteration,setting);
        NMF.cost = betadivEXP(V,NMF.Vap,setting.beta,NMF.H,setting.sparsity,setting.smoothness);
        sequentialData.W = NMF.W;
end

sequentialData.numberIteration = setting.iteration;
sequentialData.H = NMF.H;