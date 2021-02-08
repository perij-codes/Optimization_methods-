function b_L = tryR(Ig1,Ig2,L1,d1,d2,m1,m2,th1,th2,thdot1,thdot2)
%TRYR
%    B_L = TRYR(IG1,IG2,L1,D1,D2,M1,M2,TH1,TH2,THDOT1,THDOT2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    07-Feb-2021 20:22:43

t2 = cos(th1);
t3 = sin(th1);
t4 = sin(th2);
t5 = L1.^2;
t6 = d1.^2;
t7 = d2.^2;
t8 = th1.*2.0;
t9 = th2.*2.0;
t10 = thdot1.^2;
t11 = thdot2.^2;
t12 = -th2;
t14 = -thdot2;
t13 = -t9;
t15 = t12+th1;
t16 = t14+thdot1;
t20 = (Ig1.*t10)./2.0;
t21 = (Ig2.*t11)./2.0;
t23 = (m2.*t5.*t10)./2.0;
t24 = (m1.*t6.*t10)./2.0;
t25 = (m2.*t7.*t11)./2.0;
t26 = L1.*m2.*t3.*(9.81e+2./1.0e+2);
t27 = d1.*m1.*t3.*(9.81e+2./1.0e+2);
t28 = d2.*m2.*t4.*(9.81e+2./1.0e+2);
t17 = cos(t15);
t18 = sin(t15);
t19 = t8+t13;
t30 = t26+t27+t28;
t22 = sin(t19);
t29 = L1.*d2.*m2.*t17.*thdot1.*thdot2;
t31 = sign(t30);
t32 = t20+t21+t23+t24+t25+t29;
t33 = dirac(t32);
t34 = sign(t32);
b_L = [L1.*m2.*t2.*t31.*(-9.81e+2./1.0e+2)-d1.*m1.*t2.*t31.*(9.81e+2./1.0e+2)-L1.*d2.*m2.*t11.*t18.*t34+L1.*d2.*m2.*t16.*t33.*thdot1.*thdot2.*(Ig1.*t18.*thdot1.*2.0+m1.*t6.*t18.*thdot1.*2.0+m2.*t5.*t18.*thdot1.*2.0+L1.*d2.*m2.*t22.*thdot2);d2.*m2.*t31.*cos(th2).*(-9.81e+2./1.0e+2)+L1.*d2.*m2.*t10.*t18.*t34+L1.*d2.*m2.*t16.*t33.*thdot1.*thdot2.*(Ig2.*t18.*thdot2.*2.0+m2.*t7.*t18.*thdot2.*2.0+L1.*d2.*m2.*t22.*thdot1)];