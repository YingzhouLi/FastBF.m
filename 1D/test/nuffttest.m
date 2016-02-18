%Author: J. Selva. 
%Date: April 2011.

%Define a random vector with 1024 elements. 

M = 1024;
a = rand(M,1)+1i*rand(M,1);

%Call the nufft function. The output is a struct with some variable for later use.

st = nufft(a)

%The field st.aF is the fft of size st.K.

%Define a vector of random frequencies.

f = rand([2*M,1])-0.5;

%Compute the DFT at frequencies in f and the derivatives of first and second order using
%loops and direct evaluation. This takes a while
%

tic; outDL = nufft(st,f,2,'directLoop'); tDL = toc;

%The function returns the DFT values in the first column, the first derivatives in the
%second column, and so on.

outDL(1:10,:)

%Do the same computation using interpolation. This is fast.
%
tic; outIL = nufft(st,f,2); tIL = toc;

%Compare the execution times
%

[tDL,tIL]

%The executions are faster if the code is vectorized. This is an example
%
tic; outDV = nufft(st,f,2,'directVec'); tDV = toc;
tic; outIV = nufft(st,f,2,'baryVec'); tIV = toc;

[tDV, tIV]

%For larger data sizes the difference is more significant. To  check this, let us repeat
%the same experiment with M=1024*5.
%

M = 1024*5;
a = rand(M,1)+1i*rand(M,1);
st = nufft(a);

f = rand([2*M,1])-0.5;

tic; outDV = nufft(st,f,2,'directVec'); tDV = toc; %This line takes a while.
tic; outIV = nufft(st,f,2); tIV = toc;

%These are the timings:

[tDV,tIV]

%Finally, let  us check the accuracy of the interpolated values

max(abs(outDV-outIV))./max(abs(outDV))

