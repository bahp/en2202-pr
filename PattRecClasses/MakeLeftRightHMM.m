%hmm=MakeLeftRightHMM(nStates,pD,obsData,lData);%Create, initialize and train a Hidden Markov Model (a HMM object)%in a simple standardized way,%to conform with a given set of training data sequences.%%The HMM will have a first-order left-right structure,%allowing transitions only to the nearest following state.%%Input:%nStates=   desired number of HMM states%pD=        a single object of some probability-distribution class%obsData=   matrix with concatenated finite-duration training sequences.%           Observed vector samples are stored column-wise.%lData=     vector with lengths of training sub-sequences.%lData(r)=  length of r:th training sequence.%           sum(lData) == size(obsData,2)%%Result:%hmm=       the resulting Hidden Markov Model object%%Arne Leijon 2004-11-23 tested%            2011-08-04, minor fix to use HMM/initfunction hmm=MakeLeftRightHMM(nStates,pD,obsData,lData)if nStates<=0	error('Number of states must be >0');end;if nargin <4%just one single sequence    lData=size(obsData,2);end;%Left-right MarkovChain sub-object with finite duration:D=mean(lData);%average total sequence lengthD=D./nStates;%average state durationmc=initLeftRight(MarkovChain,nStates,D);hmm=init(HMM(mc,pD),obsData,lData);%hmm=InitLeftRightHMM(nStates,pD,obsData,lData);%OLDhmm=train(hmm,obsData,lData,5,log(1.01));%standard training