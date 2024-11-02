clear all
global O;
global I;
global S;
global Re;
global M;
global D1;
global D3;
global N;
global km;
global ds;
global b;
global theta;
N=1000;
[A,B1,D,x]=dftw(N);
h0=1+0.1*sin(x);
D1=real(D);
D2=real(D^2);
D3=real(D^3);
Re=3.5;
S=3000;
b=0.0;
theta=pi/180*45;
Mc=(1+3*b)*cot(theta)/3/(1+7*b);
M=0.1;
kc=sqrt((3*M*(1+7*b)+3*Re*(2/5+12/5*b+5*b^2+3*b^3))/(1+3*b)-cot(theta))/sqrt(S/3);
km=kc*(1-0.02);
c0=3*(1+2*b);
q0=1+3*b-c0;
I=eye(N);
O=ones(N);
%% psedo-archlength continuation method
[h,c,q]=travelBE(h0,c0,q0,km,N,S,Re,M,b,theta);
h0=h;
c0=c;
q0=q;
hh0=max(h0)-min(h0);
hmax=max(h0);
k0=km;
plot(x/km,h0);

km=km-0.0001;
[h,c,q]=travelBE(h0,c0,q0,km,N,S,Re,M,b,theta);
h1=h;
c1=c;
q1=q;
k1=km;
hh1=max(h1)-min(h1);
hmax1=max(h1);
% plot(x/km,h1,'r');


ds=0.02;
dh=h1-h0;
dc=c1-c0;
dq=q1-q0;
dk=k1-k0;
Nor=norm([dh' dc dq dk]);
%%information for continuation
dh=dh/Nor;
dc=dc/Nor;
dq=dq/Nor;
dk=dk/Nor;

h1=h0+dh*ds;
c1=c0+dc*ds;
q1=q0+dq*ds;
k1=k0+dk*ds;


change=1;
jj=0;
while k0>0.03
    jj=jj+1;
    while change>1e-8
        dh1=k1*D1*h1;
        dh3=k1^3*D3*h1;
        Q1=(h1.^2).*(h1+3*b)+(h1.^2).*(h1+3*b).*(S/3*dh3-cot(theta)*dh1);
        Q2=3*M*(h1+7*b).*(h1.^3).*dh1;
        Q3=3*Re*(2/5*h1.^3-1/10*b*h1.^2+5*b^2*h1+3*b^3).*(h1.^3).*dh1;
        L11=-c1*I+diag(3*h1.^2+6*b*h1);
        L12=3*M*diag((4*h1.^3+21*b*h1.^2).*dh1);
        L13=3*Re*diag((12/5*h1.^5-1/2*b*h1.^4+20*b^2*h1.^3+9*b^3*h1.^2).*dh1);
        L14=diag((3*h1.^2+6*b*h1).*(S/3*dh3-cot(theta)*dh1));
        L21=S/3*k1^3*diag((h1.^3+3*b*h1.^2))*D3;
        L22=-k1*cot(theta)*diag((h1.^3+3*b*h1.^2))*D1;
        L23=3*Re*k1*diag((2/5*h1.^3-1/10*b*h1.^2+5*b^2*h1+3*b^3).*(h1.^3))*D1;
        L24=3*M*k1*diag(h1.^4+7*b*h1.^3)*D1;
        Lk1=3*M*(h1+7*b).*(h1.^3).*(D1*h1);
        Lk2=3*Re*(2/5*h1.^3-1/10*b*h1.^2+5*b^2*h+3*b^3).*(h1.^3).*(D1*h1);
        Lk3=-cot(theta)*(h1+3*b).*(h1.^2).*(D1*h1);
        Lk4=S*k1^2*(h1+3*b).*(h1.^2).*(D3*h1);
        F=-c1*h1+Q1+Q2+Q3-q1;
        G=(h1-h0)'*dh+(c1-c0)*dc+(q1-q0)*dq+(k1-k0)*dk-ds;
        % closed flow system
        Mass=sum(h1)/N-1;
        A=[L11+L12+L13+L14+L21+L22+L23+L24,-h1,-O(:,1),Lk1+Lk2+Lk3+Lk4;dh',dc,dq,dk;1/N*O(1,:),0,0,0;dh1',0,0,0];
        bu=[F;G;Mass;0];
        dhdcdqdk=-A\bu;
        h1=h1+dhdcdqdk(1:N);
        c1=c1+dhdcdqdk(N+1);
        q1=q1+dhdcdqdk(N+2);
        k1=k1+dhdcdqdk(N+3);
        change=norm(dhdcdqdk,inf);
    end
    change=1;
    h0=h1;
    c0=c1;
    q0=q1;
    k0=k1
    ff=[0*h1;1;0;0];
    dF=A\ff;
    dF=dF/norm(dF,2);
    dh=dF(1:N);
    dc=dF(N+1);
    dq=dF(N+2);
    dk=dF(N+3);
    h1=h0+dh*ds;
    c1=c0+dc*ds;
    q1=q0+dq*ds;
    k1=k0+dk*ds;
    hh(jj)=max(h0)-min(h0);
    hhmax(jj)=max(h0);
    cc(jj)=c0;
    kk(jj)=k0;
    jj
end
disp(['k1=',num2str(k1),';k0=',num2str(k0),';c1=',num2str(c1),';c0=',num2str(c0),';'])
figure(1),plot(cc,2*pi./kk,'r'),hold on;
figure(2),plot(kk,hh),hold on; plot(kk,hhmax)
figure(3),plot(x/k0,h0,'r')


