function varargout = nufft(varargin)

%Author: J. Selva
%Date: April 2011.
%
%nufft computes the FFT at arbitrary frequencies using barycentric interpolation.
%
%For the interpolation method, see
%
%J. Selva, 'Design of barycentric interpolator for uniform and nonuniform sampling
%  grids', IEEE Trans. on Signal Processing, vol. 58, n. 3, pp. 1618-1627, March 2010.
%
%Usage:
%
% -First call nufft with the vector of samples as argument,
%
%       st = nufft(a); %a is the vector of samples.
%
%    The output is an struct. The field st.aF is the DFT of vector a, sampled with spacing
%    1/st.K. 
%
%    The DFT is defined using
%
%        A_k = \sum_{n=m_1}^{m_1+M-1} a_n e^{-j 2 \pi n f}
%
%    where m_1=-ceil(M/2) and M the number of samples. 
%
% -Then call nufft with sintax 
%
%        out = nufft(st,f,nd,method)
%
%   where
%
%     f: Vector of frequencies where the DFT is evaluated. Its elements must follow
%        abs(f(k))<=0.5
%
%     nd: Derivative order. nufft computes derivatives up to order nd. At the output 
%         out(p,q) is the (q-1)-order derivative of the DFT at frequency f(p). The first
%         column of out contains the zero-order derivatives, i.e, the values of the DFT at
%         frequencies in vector f. 
%     method: If omitted, it is 'baryVec'. Four methods have been implemented:
%
%         +'baryLoop': The output is interpolated using the barycentric method and a loop
%           implementation.
%
%         +'baryVec': The output is interpolated using the barycentric method and a
%          vectorized implementation.
%
%         +'directLoop': The output is computed using the DFT sum directly and a loop
%           implementation.
%  
%         +'directVec': The output is computed using the DFT sum directly and a vectorized
%           implementation.

 
if ~isstruct(varargin{1})
  st = [];
  st.a = varargin{1};
  st.a = st.a(:);
  
  st.M = length(st.a);
  st.m1 = -ceil(st.M/2);
  st.K =  2^ceil(log2(st.M)+1);
  st.aF = fft([st.a(-st.m1+1:st.M); zeros(st.K-st.M,1); st.a(1:-st.m1)]);
  
  errb = 10.^-13; %This is the interpolation accuracy. You can set any larger value, and ...
      %this would reduce st.P. The number of interpolation samples is 2*st.P+1.
  
  st.T = 1/st.K;
  st.B = -2*st.m1;
  st.P = ceil(acsch(errb)/(pi*(1-st.B*st.T)));
  st.vt = MidElementToEnd((-st.P:st.P)*st.T);
  
  st.alpha =  MidElementToEnd(baryWeights(st.T,st.B,st.P));
  
  varargout = {st};
  
  return
end

[st,f] = deal(varargin{1:2});
  
nd = 0;
   
if nargin > 2
  nd = varargin{3};
end
  
method = 'baryVec';
  
if nargin > 3
  method = varargin{4};
end
  
Nf = length(f);
out = zeros(Nf,nd+1);

switch method
    
  case 'baryLoop' %Interpolated method using loops

    for k = 1:length(f)

      x = f(k);
      
      n = floor(x/st.T+0.5);
      u = x - n * st.T;

      da = MidElementToEnd(st.aF(1+mod(n-st.P:n+st.P,st.K)).');

      out(k,:) = DerBarycentricInterp3(st.alpha,da,st.vt,u,nd);

    end

  case 'baryVec' %Vectorized interpolated method
    
    f = f(:);
    Nf = length(f);
    
    n = floor(f/st.T+0.5);
    u = f - n * st.T;
    
    pr = [-st.P:-1 , 1:st.P , 0];
    
    ms = st.aF(1+mod(n(:,ones(1,2*st.P+1)) + pr(ones(Nf,1),:),st.K));
    
    if length(f) == 1
      ms = ms.';
    end
    
    out = DerBarycentricInterp3Vec(st.alpha,ms,st.vt,u,nd);
    
  case 'directLoop' % Direct method using loops
      
    for k = 1:length(f)
      out(k,1) = 0;
	
      for r = st.m1:st.m1+st.M-1
	out(k,1) = out(k,1)+exp(-1i*2*pi*r*f(k))*st.a(r-st.m1+1);
      end
	
      for kd = 1:nd
	out(k,kd+1) = 0;
	
	for r = st.m1:st.m1+st.M-1
	  out(k,kd+1) = out(k,kd+1) + ...
	      ((-1i*2*pi*r).^kd .* exp(-1i*2*pi*r*f(k)))*st.a(r-st.m1+1);
	end
      
      end
    end  
    
  case 'directVec' %Vectorized direct method
	
    for k = 1:length(f)
      out(k,1) = exp(-1i*2*pi*(st.m1:st.m1+st.M-1)*f(k))*st.a;
      
      for kd = 1:nd
	out(k,kd+1) = ...
	    ((-1i*2*pi*(st.m1:st.m1+st.M-1)).^kd .* ...
	    exp(-1i*2*pi*(st.m1:st.m1+st.M-1)*f(k)))*st.a;
      end
      
    end
      
end
 
varargout = {out};

function v = MidElementToEnd(v)

ind = ceil(length(v)/2);
v = [v(1:ind-1),v(ind+1:end),v(ind)];

function v = APPulse(t,B,TSL)

v = real(sinc(B*sqrt(t.^2-TSL^2)))/real(sinc(1i*pi*B*TSL));

function w = baryWeights(T,B,P)

vt = (-P:P)*T;
g = APPulse(vt,1/T-B,T*(P+1));
gam = gamma(vt/T+P+1) .* gamma(-vt/T+P+1) .* g;

N = length(vt);
LD = ones(1,N);

for k = 1:N
  LD(k) = prod(vt(k)-vt(1:k-1))* prod(vt(k)-vt(k+1:N));
end

w = gam ./ LD;
w = w / max(abs(w));

function out = DerBarycentricInterp3(alpha,s,t,tau,nd)

vD = 0;
Nt = length(t);
LF = [1,zeros(1,Nt-1)];
out = zeros(1,nd+1);

for k = 1:Nt-1
  vD = vD*(tau-t(k))+alpha(k)*LF(k);
  LF(k+1) = LF(k)*(tau-t(k));
end

vD = vD * (tau-t(Nt)) + alpha(Nt) * LF(Nt);

z = s;

for kd = 0:nd

  z1 = z-z(end); cN = 0;

  for k = 1:Nt-1
    cN = cN * (tau-t(k))+z1(k)*alpha(k)*LF(k);
  end
  cN = cN/vD;
  
  ztau = z(end)+(tau-t(end))*cN;
  out(kd+1) = ztau;
  
  if kd < nd
    z = [ (z(1:end-1)-ztau)./(t(1:end-1)-tau) , cN ] * (kd+1);
  end
  
end

function out = DerBarycentricInterp3Vec(alpha,zL,t,tauL,nd)

NtauL = size(tauL,1);
LBtau = 200; %Block size for vectorized processing. Change this to adjust the memory 
             %usage. 

NBlock = 1+floor(NtauL/LBtau);

Nt = length(t);

out = zeros(NtauL,nd+1);

for r = 0:NBlock-1
  ind1 = 1+r*LBtau;
  ind2 = min([(r+1)*LBtau,NtauL]);
  
  Ntau = ind2-ind1+1;
  z = zL(ind1:ind2,:);
  tau = tauL(ind1:ind2);
  
  vD = zeros(Ntau,1);
  
  LF = [1,zeros(1,Nt-1)];
  LF = LF(ones(Ntau,1),:);

  for k = 1:Nt-1
    vD = vD .* (tau-t(k))+alpha(k)*LF(:,k);
    LF(:,k+1) = LF(:,k) .* (tau-t(k));
  end

  vD = vD .* (tau-t(Nt)) + alpha(Nt) * LF(:,Nt);

  for kd = 0:nd

    pr = z(:,end); z1 = z-pr(:,ones(1,Nt)); cN = zeros(Ntau,1);

    for k = 1:Nt-1
      cN = cN .* (tau-t(k)) + z1(:,k) .* alpha(k) .* LF(:,k);
    end
    cN = cN ./ vD;
  
    ztau = z(:,end)+(tau-t(end)) .* cN;
    out(ind1:ind2,kd+1) = ztau;
  
    if kd < nd
      pr1 = ztau;
      pr1 = z(:,1:end-1) - pr1(:,ones(1,Nt-1));
      
      pr2 = t(1:end-1);
      pr2 = pr2(ones(Ntau,1),:)-tau(:,ones(1,Nt-1));
      
      z = [pr1 ./ pr2, cN] * (kd+1);
      
    end
  end
  
end