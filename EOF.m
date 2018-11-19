%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical orthogonal function analysis of SCN slice data
% 2018/11/19 by Isao Tokuda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

transient=40;  % transient time to be removed from analysis
skip=4;        % number of points to skip in the analysis

%%%%%%%%%%%%% Reading the movie data  %%%%%%%%%%
% Sample can be downloaded from 
% https://www.ebi.ac.uk/biostudies/files/S-BSST212/20130208_PER2LUC_A_SNR10.avi
% (sampling time interval of 1 h) 

fileinfo = VideoReader('20130208_PER2LUC_A_SNR10.avi', 'tag', 'myreader1');
vidFrames = read(fileinfo);
                       
for k = 1 : fileinfo.NumberOfFrames
  mov(k).cdata = vidFrames(:,:,:,k);
  mov(k).colormap = [];
end

DataLength=fileinfo.NumberOfFrames-transient;

cell=0;
for i=1:skip:fileinfo.Height-1,
  for j=1:skip:fileinfo.Width-1,
    for n=1:DataLength,
      z(n,1)=cast(mov(transient+n).cdata(i,j,1),'double');
    end
    %%%%%%%%%%%%% Pick up pixels having large amplitudes  
    if std(z)>30,
      cell=cell+1;
      for t=1:1:DataLength,
        ts(cell,t)=cast(mov(transient+t).cdata(i,j,1),'double');
      end
      ylocation(cell)=i;
      xlocation(cell)=j;
    end
  end
end


%%%%%%%%% Period analysis by Chi-square periodogram %%%%%%%%% 

for i=1:1:cell,
  [cycle,slevel] = myChiSquarePeriodogram(transpose(ts(i,:)));
  if slevel < 0.1 
    period(i)=cycle;
  else
    period(i)=0;
  end
end


%%%%%%%%% Computation of acrophase by Cosinor fitting %%%%%%%%% 

per=24.0;
z=zeros(DataLength,1);
for t=1:1:DataLength,
  for i=1:1:cell,
    z(t,1)=z(t,1)+ts(i,t);
  end
end
[phi,R] = Cosinor4(z,per);
meanphase = per*phi/(2*pi); 

for i=1:1:cell,
  [phi,R] = Cosinor4(ts(i,:),per);
  if R>0.8,
    phase(i) = per*phi/(2*pi) - meanphase;
    if phase(i) > 0.5*per,
      phase(i) = phase(i)-per;
    end
    if phase(i) < -0.5*per,
      phase(i) = phase(i)+per;
    end
  else
    phase(i) = -per;
  end
end


%%%%%%%%% Detrend and normalize data %%%%%%%%%%%
ts = transpose(ts);

samp=1.0; %[hour]
Fs = 1.0/samp;
t = 1:size(ts,1);
t = t/Fs;
t = transpose(t);

for i=1:size(ts,2)
  [fdata fcor ftype(i,:)] = myTrendPol2([t ts(:,i)],2);
  ts(:,i) = fdata(:,2);
end  

for i=1:size(ts,2),
  ts(:,i) = ts(:,i) - mean(ts(:,i));
  ts(:,i) = ts(:,i)/std(ts(:,i));
end

%%%%%%%%% Principle component analysis %%%%%%%%%%%
X=ts;
V = cov(X);
[COEFF, LATENT, EXPLAINED] = pcacov(V);


%%%%%%%%% Plot eigenvalues %%%%%%%%%%%
figure
plot(EXPLAINED(1:10),'bo-','MarkerSize',3)
xlabel('Components')
ylabel('Eigenvalues [%]')


%%%%%%%%% Plot eigenmodes %%%%%%%%%%%
meanX=mean(X,1);
for i=1:size(X,1),
  X(i,:)=X(i,:)-meanX;
end

Y=X*COEFF;
for i=1:size(Y,2),
  Y(:,i) = Y(:,i) - mean(Y(:,i));
  Y(:,i) = Y(:,i)/std(Y(:,i));
end
figure
plot(t,Y(:,1),'r-',t,Y(:,2),'g-')
xlabel('Time [h]')
ylabel('Eigenmode')

%%%%%%%%% Drawing slice locations of cells for 1st and 2nd modes %%%%%%%%%%%
portion=0.4;
for i=1:1:fileinfo.Height,
  for j=1:1:fileinfo.Width,
  C(i,j,1)=cast(0,'uint8');
  C(i,j,2)=cast(0,'uint8');
  C(i,j,3)=cast(0,'uint8');
  end
end

[Y,I]=sort(abs(COEFF(:,1)));
for i=1:1:portion*0.01*EXPLAINED(1)*cell,
  n=I(cell+1-i);
  for j=-1:1:1,
    for k=-1:1:1,
      C(ylocation(n)+j,xlocation(n)+k,1)=255;
    end
  end
end

[Y,I]=sort(abs(COEFF(:,2)));
for i=1:1:portion*0.01*EXPLAINED(2)*cell,
  n=I(cell+1-i);
  for j=-1:1:1,
    for k=-1:1:1,
      C(ylocation(n)+j,xlocation(n)+k,1)=0;
      C(ylocation(n)+j,xlocation(n)+k,2)=255;
    end
  end
end

figure
image(C)
xlabel('X-axis')
ylabel('Y-axis')


%%%%%%%%% Histogram of cellular periods for 1st and 2nd modes %%%%%%%%%%%
[Y,I]=sort(abs(COEFF(:,1)));
per1=[];
for i=1:1:portion*0.01*EXPLAINED(1)*cell,
  per1=[per1 period(I(cell+1-i))];
end
[Y,I]=sort(abs(COEFF(:,2)));
per2=[];
for i=1:1:portion*0.01*EXPLAINED(2)*cell,
  per2=[per2 period(I(cell+1-i))];
end
[Y,I]=sort(abs(COEFF(:,3)));

figure
x = 15:1:35;
elements1 = hist(per1,x);
elements1 = 100*hist(per1,x)/(size(ts,2));
elements2 = hist(per2,x);
elements2 = 100*hist(per2,x)/(size(ts,2));
y=[elements1; elements2]; 
b = bar(x,transpose(y),'stacked','r','g');
b(1).FaceColor = 'r'; 
b(2).FaceColor = 'g'; 
xlim([20 30])
xlabel('Period [h]')
ylabel('Histogram [%]')

%%%%%%%%% Histogram of cellular acrophases for 1st and 2nd modes %%%%%%%%%%%
[Y,I]=sort(abs(COEFF(:,1)));
phase1=[];
for i=1:1:portion*0.01*EXPLAINED(1)*cell,
  phase1=[phase1 phase(I(cell+1-i))];
end
[Y,I]=sort(abs(COEFF(:,2)));
phase2=[];
for i=1:1:portion*0.01*EXPLAINED(2)*cell,
  n=I(cell+1-i);
  phase2=[phase2 phase(I(cell+1-i))];
end

figure
x = -12:1:12;
elements1 = hist(phase1,x);
elements1 = 100*hist(phase1,x)/(size(ts,2));
elements2 = hist(phase2,x);
elements2 = 100*hist(phase2,x)/(size(ts,2));
y=[elements1; elements2]; 
b = bar(x,transpose(y),'stacked','r','g');
b(1).FaceColor = 'r';
b(2).FaceColor = 'g'; 
xlim([-6 6])
xlabel('Phase [h]')
ylabel('Histogram [%]')

%%%%%%%%% Time traces for cells showing 1st mode %%%%%%%%%%%
figure
[Y,I]=sort(abs(COEFF(:,1)));
for i=1:1:portion*(EXPLAINED(1)/100)*cell,
  n=I(cell+1-i);
  plot(t,ts(:,n),'r-');
  hold on;
end
xlabel('Time [h]')
ylabel('Bioluminescence signals (1st mode)')

%%%%%%%%% Time traces for cells showing 2nd mode %%%%%%%%%%%
figure
[Y,I]=sort(abs(COEFF(:,2)));
for i=1:1:portion*(EXPLAINED(2)/100)*cell,
  n=I(cell+1-i);
  plot(t,ts(:,n),'g-');
  hold on;
end
xlabel('Time [h]')
ylabel('Bioluminescence signals (2nd mode)')

