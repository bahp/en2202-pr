function S=rand(mc,T)
%S=rand(mc,T) returns a random state sequence from given MarkovChain object.
%
%Input:
%mc=    a single MarkovChain object
%T= scalar defining maximum length of desired state sequence.
%   An infinite-duration MarkovChain always generates sequence of length=T
%   A finite-duration MarkovChain may return shorter a sequence,
%   if the END state was reached before T samples.
%
%Result:
%S= integer row vector with random state sequence,
%   NOT INCLUDING the END state,
%   even if encountered within T samples
%If mc has INFINITE duration,
%   length(S) == T
%If mc has FINITE duration,
%   length(S) <= T
%
%---------------------------------------------
%Code Authors:
%---------------------------------------------

S=zeros(1,T);%space for resulting row vector
nS=mc.nStates;

d = DiscreteD(mc.InitialProb);
S(1) = d.rand(1);

if ~mc.finiteDuration()
  for i = 2:T
    d = DiscreteD(mc.TransitionProb(S(i-1), :));
    S(i) = d.rand(1);
  end
else
  i = 2;
  stop = false;
  while i <= T && ~stop
    d = DiscreteD(mc.TransitionProb(S(i-1), :));
    r = d.rand(1);
    if r == nS+1
      stop = true;
      S = S(1:i-1);
    else
      S(i) = r;
      i = i+1;
    end
  end
end


