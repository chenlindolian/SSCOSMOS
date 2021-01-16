function weight = AverageEchoWeight(GREMag,TE)
    
  Magtemp = zeros(size(GREMag));
  weight = zeros(size(GREMag));
  
  for nEchos = 1:length(TE)
    Magtemp(:,:,:,nEchos) = TE(nEchos).*GREMag(:,:,:,nEchos);
  end
  Magtotal = sum(Magtemp,4);
  e = 0.00001;
  
  for nEchos = 1:length(TE)
    weight(:,:,:,nEchos) = Magtemp(:,:,:,nEchos)./ (Magtotal+e);  
  end
end