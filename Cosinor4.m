function [max_phi,max_R] = Cosinor4(activity,Tcycle)

  vt=zeros(3,1,'double');
  for l=1:length(activity),
    vt(1) = vt(1)+activity(l)*sin(2.0*pi*l/Tcycle);
    vt(2) = vt(2)+activity(l)*cos(2.0*pi*l/Tcycle);
    vt(3) = vt(3)+activity(l);
  end
  mx=zeros(3,3,'double');
  for l=1:length(activity),
    mx(1,1) = mx(1,1)+sin(2.0*pi*l/Tcycle)*sin(2.0*pi*l/Tcycle);
    mx(1,2) = mx(1,2)+sin(2.0*pi*l/Tcycle)*cos(2.0*pi*l/Tcycle);
    mx(1,3) = mx(1,3)+sin(2.0*pi*l/Tcycle);
    mx(2,2) = mx(2,2)+cos(2.0*pi*l/Tcycle)*cos(2.0*pi*l/Tcycle);
    mx(2,3) = mx(2,3)+cos(2.0*pi*l/Tcycle);
    mx(3,3) = mx(3,3)+1.0;
  end
  mx(2,1)=mx(1,2); mx(3,1)=mx(1,3); mx(3,2)=mx(2,3);
  w=inv(mx)*vt;
  Amp=sqrt(w(1)*w(1)+w(2)*w(2));
  if 0<=w(1),
    phi=asin(w(2)/Amp);
  else
    phi=pi-asin(w(2)/Amp);
  end

  SS_error=0.0; SS_regress=0.0; SS_T=0.0;
  for l=1:length(activity),
    y(l)=w(1)*sin(2.0*pi*l/Tcycle)+w(2)*cos(2.0*pi*l/Tcycle)+w(3);
    y(l)=Amp*sin(phi+2.0*pi*l/Tcycle)+w(3);
    SS_error=SS_error+(activity(l)-y(l))^2;
    SS_regress=SS_regress+(y(l)-mean(activity))^2;
    SS_T=SS_T+(activity(l)-mean(activity))^2;
  end
  F=(SS_regress/2)/(SS_error/(length(activity)-2-1));

  max_F=F; 
  max_cycle=Tcycle; 
  max_phi=phi; 
  max_R=1.0-(SS_error/SS_T);


