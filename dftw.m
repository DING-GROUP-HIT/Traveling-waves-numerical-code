function [A,B,D,x]=dftw(N)
%x [0 2pi]
dxx=2*pi/(N);x=dxx*((1:N)');
kk=[-N/2+1:N/2]';
C=diag(1i*kk);
A=exp(1i*x*kk');
B=exp(-1i*kk*x')/N;
D=A*(C*B);
%D=real(D);
