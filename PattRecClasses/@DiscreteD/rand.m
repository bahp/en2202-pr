function R=rand(pD,nData)
%R=rand(pD,nData) returns random scalars drawn from given Discrete Distribution.
%
%Input:
%pD=    DiscreteD object
%nData= scalar defining number of wanted random data elements
%
%Result:
%R= row vector with integer random data drawn from the DiscreteD object pD
%   (size(R)= [1, nData]
%
%------------------------------------------------------
%Code Authors: Bernard Hernandes and Fernando Iglesias
%------------------------------------------------------

if numel(pD)>1
    error('Method works only for a single DiscreteD object');
end;

if nData < 30 % simple heuristic to choose vectorized/non-vectorized
  % Non-vectorized implementation
  cs = cumsum(pD.ProbMass);
  s = rand([1 nData]);
  R = zeros(1, nData); % memory pre-allocation
  for i = 1:nData
   R(i) = find(cs > s(i), 1);
  end
else
  % Vectorized implementation
  cs = cumsum(pD.ProbMass);
  s = rand([1 nData]);
  CS = repmat(cs,1,nData);
  S = repmat(s,length(cs),1);
  R = CS > S;
  [~, R] = max(R);
end