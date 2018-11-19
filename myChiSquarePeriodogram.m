function [per,slevel] = myChiSquarePeriodogram(ts)

  samp=1.0; % sampling time interval of 1 h
  MaxP=42;  % maximum period to estimate
  MinP=14;  % minimum period to estimate
  duration=size(ts,1);

  percent=zeros(MaxP,1,'double');
  cycle=zeros(MaxP,1,'double');

  Qp=zeros(MaxP,1,'double');
  for Period=MinP:MaxP,
    K=floor(duration/Period);
    N=K*Period;
    M=mean(ts(1:N));
    var(ts(1:N));

    Mh=zeros(Period,1,'double');
    for i=1:Period,
      for j=0:K-1,
        Mh(i)=Mh(i)+ts(i+j*Period);
      end
      Mh(i)=Mh(i)/K;
    end

    for i=1:Period,
      Qp(Period)=Qp(Period)+(Mh(i)-M)^2;
    end
    clear Mh
    Qp(Period)=K*Qp(Period)/var(ts(1:N)); %Chi-square statistics
    [tail,percentile] = myChiSquareDistribution(Period-1,Qp(Period));
    percent(Period)=percentile;
  end

  minQp=zeros(MaxP,1,'double');
  for Period=MinP:MaxP,
    if Period==MinP, 
      minQp(Period)=0; 
    else
      minQp(Period)=minQp(Period-1);
    end
    [tail,percentile] = myChiSquareDistribution(Period-1,minQp(Period));
    while percentile <= 0.99,
      minQp(Period)=minQp(Period)+0.01;
      [tail,percentile] = myChiSquareDistribution(Period-1,minQp(Period));
    end 
  end

  maxpower=0;
  for k = MinP:MaxP,
    if maxpower<(Qp(k)/minQp(k)) & 0<minQp(k),
      maxpower=Qp(k)/minQp(k);
      per=samp*k;
      [tail,percentile] = myChiSquareDistribution(k-1,Qp(k));
      slevel=1.0-percentile;
    end
  end

  %figure 
  %for k = 1:MaxP,
  %  cycle(k) = samp*k;
  %end
  %plot(cycle,Qp,'b',cycle,minQp,'r')
  %xlabel('Period [days]')
  %ylabel('Chi-Square Periodogram')

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [tail,percentile] = myChiSquareDistribution(v,w)

  %Asks the user for the relevant input
  %v = input('Degrees of freedom: ');
  %w = input('Chi-square: '); 

  % Calculates the denominator product
  r = 1;
  for i = v : -2 : 2
    r = r*i;
  end
  % Calculates the numerator product
  k = w^(floor((v+1)/2))*exp(-w/2)/r; 

  % If degrees of freedom are odd, then uses the pi factor
  if floor(v/2) == v/2
    j = 1;
  else
    j = sqrt(2/(w*pi));
  end
  l = 1;
  m = 1;
  v = v + 2;
  m = m*w/v; 

  % Summation factor
  while m >= 1e-7,
    l = l + m;
    v = v + 2;
    m = m*w/v;
  end 

  tail=1-j*k*l;
  percentile=j*k*l;

end
