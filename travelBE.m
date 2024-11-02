function [h,c,q]=travelBE(h0,c0,q0,km,N,S,Re,M,b,theta)
[A,B1,D,x]=dftw(N);
D1=real(D);
D3=real(D^3);
I=eye(N);
O=ones(N);
change=1;
while change>1e-8
    dh1=km*D1*h0;
    dh3=km^3*D3*h0;
    Q1=h0.^2.*(h0+3*b)+(h0.^2).*(h0+3*b).*(S/3*dh3-cot(theta)*dh1);
    Q2=3*M*(h0+7*b).*(h0.^3).*dh1;
    Q3=3*Re*(2/5*h0.^3+12/5*b*h0.^2+5*b^2*h0+3*b^3).*(h0.^3).*dh1;
    L11=-c0*I+diag(3*h0.^2+6*b*h0);
    L12=3*M*diag((4*h0.^3+21*b*h0.^2).*dh1);
    L13=3*Re*diag((12/5*h0.^5+12*b*h0.^4+20*b^2*h0.^3+9*b^3*h0.^2).*dh1);
    L14=diag((3*h0.^2+6*b*h0).*(S/3*dh3-cot(theta)*dh1));
    L21=S/3*km^3*diag((h0.^3+3*b*h0.^2))*D3;
    L22=-km*cot(theta)*diag((h0.^3+3*b*h0.^2))*D1;
    L23=3*Re*km*diag((2/5*h0.^3+12/5*b*h0.^2+5*b^2*h0+3*b^3).*(h0.^3))*D1;
    L24=3*M*km*diag(h0.^4+7*b*h0.^3)*D1;
    F=-c0*h0+Q1+Q2+Q3-q0;
   %% traveling wave
   %closed system
   G=sum(h0)/N-1;
   A=[L11+L12+L13+L14+L21+L22+L23+L24,-h0,-O(:,1);1/N*O(1,:),0,0;dh1',0,0];
    bu=-[F;G;0];
    dhdcdq=A\bu;
    h0=h0+dhdcdq(1:N);
    c0=c0+dhdcdq(N+1);
    q0=q0+dhdcdq(N+2);
    change=norm(dhdcdq(1:N),inf)
end
h=h0;
c=c0;
q=q0;
