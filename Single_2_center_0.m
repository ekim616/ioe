n3 = 2.2119;
n1 = 2.1386;
n = sqrt (n3*n1);
delta_ngr = 0.08;

period = 21;
lc = (179 + 3/4)*period;
lpa = (20)*period;
lpb = (20 + 2/4)*period;
L = lc*8 + lpa*4 + lpb*3;

lambda0 = period*(n3 - n1);
lambdaj = lambda0 - 0.0016*10;

npoints = 2000;

S = pi/(2*3.02*10^4);

for int = 1 : npoints

lambda = 1.37 + 0.0001*int;
delta = 2*pi*(n3 - n1)*(1/lambda - 1/lambda0);
deltaj = 2*pi*delta_ngr*(1/lambdaj - 1/lambda0);
d = delta/2;
b1 = 2*pi/lambda*n1;
b2 = 2*pi/lambda*n3;

K1 = S*cos (deltaj*1/16*L);
K1a = S*cos (deltaj*1/16*L)*(1 + 0.5*cos (2*pi*(1/16 - 0.5)));
ac1 = cos (sqrt (K1^2 + d^2)*lc) + i*d/sqrt (K1^2 + d^2)*sin (sqrt (K1^2 + d^2)*lc);
bc1 = -i*K1/sqrt (K1^2 + d^2)*sin (sqrt (K1^2 + d^2)*lc);
ac1a = cos (sqrt (K1a^2 + d^2)*lc) + i*d/sqrt (K1a^2 + d^2)*sin (sqrt (K1a^2 + d^2)*lc);
bc1a = -i*K1a/sqrt (K1a^2 + d^2)*sin (sqrt (K1a^2 + d^2)*lc);
MC1 = [ac1*exp(-i*(b1 + d)*lc) bc1*exp(-i*(b1 + d)*lc); -conj(bc1)*exp(-i*(b2 - d)*lc) conj(ac1)*exp(-i*(b2 - d)*lc)];
MC1a = [ac1a*exp(-i*(b1 + d)*lc) bc1a*exp(-i*(b1 + d)*lc); -conj(bc1a)*exp(-i*(b2 - d)*lc) conj(ac1a)*exp(-i*(b2 - d)*lc)];

MP2 = [exp(-i*b1*lpa) 0; 0 exp(-i*b2*lpa)];
K2 = S*sin (deltaj*3/16*L);
K2a = S*sin (deltaj*3/16*L)*(1 + 0.5*cos (2*pi*(3/16 - 0.5)));
ac2 = cos (sqrt (K2^2 + d^2)*lc) + i*d/sqrt (K2^2 + d^2)*sin (sqrt (K2^2 + d^2)*lc);
bc2 = -i*K2/sqrt (K2^2 + d^2)*sin (sqrt (K2^2 + d^2)*lc);
ac2a = cos (sqrt (K2a^2 + d^2)*lc) + i*d/sqrt (K2a^2 + d^2)*sin (sqrt (K2a^2 + d^2)*lc);
bc2a = -i*K2a/sqrt (K2a^2 + d^2)*sin (sqrt (K2a^2 + d^2)*lc);
MC2 = [ac2*exp(-i*(b1 + d)*lc) bc2*exp(-i*(b1 + d)*lc); -conj(bc2)*exp(-i*(b2 - d)*lc) conj(ac2)*exp(-i*(b2 - d)*lc)];
MC2a = [ac2a*exp(-i*(b1 + d)*lc) bc2a*exp(-i*(b1 + d)*lc); -conj(bc2a)*exp(-i*(b2 - d)*lc) conj(ac2a)*exp(-i*(b2 - d)*lc)];
T2 = MC2*MP2;
T2a = MC2a*MP2;

MP3 = [exp(-i*b1*lpb) 0; 0 exp(-i*b2*lpb)];
K3 = S*cos (deltaj*5/16*L);
K3a = S*cos (deltaj*5/16*L)*(1 + 0.5*cos (2*pi*(5/16 - 0.5)));
ac3 = cos (sqrt (K3^2 + d^2)*lc) + i*d/sqrt (K3^2 + d^2)*sin (sqrt (K3^2 + d^2)*lc);
bc3 = -i*K3/sqrt (K3^2 + d^2)*sin (sqrt (K3^2 + d^2)*lc);
ac3a = cos (sqrt (K3a^2 + d^2)*lc) + i*d/sqrt (K3a^2 + d^2)*sin (sqrt (K3a^2 + d^2)*lc);
bc3a = -i*K3a/sqrt (K3a^2 + d^2)*sin (sqrt (K3a^2 + d^2)*lc);
MC3 = [ac3*exp(-i*(b1 + d)*lc) bc3*exp(-i*(b1 + d)*lc); -conj(bc3)*exp(-i*(b2 - d)*lc) conj(ac3)*exp(-i*(b2 - d)*lc)];
MC3a = [ac3a*exp(-i*(b1 + d)*lc) bc3a*exp(-i*(b1 + d)*lc); -conj(bc3a)*exp(-i*(b2 - d)*lc) conj(ac3a)*exp(-i*(b2 - d)*lc)];
T3 = MC3*MP3;
T3a = MC3a*MP3;

MP4 = [exp(-i*b1*lpa) 0; 0 exp(-i*b2*lpa)];
K4 = S*sin (deltaj*7/16*L);
K4a = S*sin (deltaj*7/16*L)*(1 + 0.5*cos (2*pi*(7/16 - 0.5)));
ac4 = cos (sqrt (K4^2 + d^2)*lc) + i*d/sqrt (K4^2 + d^2)*sin (sqrt (K4^2 + d^2)*lc);
bc4 = -i*K4/sqrt (K4^2 + d^2)*sin (sqrt (K4^2 + d^2)*lc);
ac4a = cos (sqrt (K4a^2 + d^2)*lc) + i*d/sqrt (K4a^2 + d^2)*sin (sqrt (K4a^2 + d^2)*lc);
bc4a = -i*K4a/sqrt (K4a^2 + d^2)*sin (sqrt (K4a^2 + d^2)*lc);
MC4 = [ac4*exp(-i*(b1 + d)*lc) bc4*exp(-i*(b1 + d)*lc); -conj(bc4)*exp(-i*(b2 - d)*lc) conj(ac4)*exp(-i*(b2 - d)*lc)];
MC4a = [ac4a*exp(-i*(b1 + d)*lc) bc4a*exp(-i*(b1 + d)*lc); -conj(bc4a)*exp(-i*(b2 - d)*lc) conj(ac4a)*exp(-i*(b2 - d)*lc)];
T4 = MC4*MP4;
T4a = MC4a*MP4;

MP5 = [exp(-i*b1*lpb) 0; 0 exp(-i*b2*lpb)];
K5 = S*cos (deltaj*9/16*L);
K5a = S*cos (deltaj*9/16*L)*(1 + 0.5*cos (2*pi*(9/16 - 0.5)));
ac5 = cos (sqrt (K5^2 + d^2)*lc) + i*d/sqrt (K5^2 + d^2)*sin (sqrt (K5^2 + d^2)*lc);
bc5 = -i*K5/sqrt (K5^2 + d^2)*sin (sqrt (K5^2 + d^2)*lc);
ac5a = cos (sqrt (K5a^2 + d^2)*lc) + i*d/sqrt (K5a^2 + d^2)*sin (sqrt (K5a^2 + d^2)*lc);
bc5a = -i*K5a/sqrt (K5a^2 + d^2)*sin (sqrt (K5a^2 + d^2)*lc);
MC5 = [ac5*exp(-i*(b1 + d)*lc) bc5*exp(-i*(b1 + d)*lc); -conj(bc5)*exp(-i*(b2 - d)*lc) conj(ac5)*exp(-i*(b2 - d)*lc)];
MC5a = [ac5a*exp(-i*(b1 + d)*lc) bc5a*exp(-i*(b1 + d)*lc); -conj(bc5a)*exp(-i*(b2 - d)*lc) conj(ac5a)*exp(-i*(b2 - d)*lc)];
T5 = MC5*MP5;
T5a = MC5a*MP5;

MP6 = [exp(-i*b1*lpa) 0; 0 exp(-i*b2*lpa)];
K6 = S*sin (deltaj*11/16*L);
K6a = S*sin (deltaj*11/16*L)*(1 + 0.5*cos (2*pi*(11/16 - 0.5)));
ac6 = cos (sqrt (K6^2 + d^2)*lc) + i*d/sqrt (K6^2 + d^2)*sin (sqrt (K6^2 + d^2)*lc);
bc6 = -i*K6/sqrt (K6^2 + d^2)*sin (sqrt (K6^2 + d^2)*lc);
ac6a = cos (sqrt (K6a^2 + d^2)*lc) + i*d/sqrt (K6a^2 + d^2)*sin (sqrt (K6a^2 + d^2)*lc);
bc6a = -i*K6a/sqrt (K6a^2 + d^2)*sin (sqrt (K6a^2 + d^2)*lc);
MC6 = [ac6*exp(-i*(b1 + d)*lc) bc6*exp(-i*(b1 + d)*lc); -conj(bc6)*exp(-i*(b2 - d)*lc) conj(ac6)*exp(-i*(b2 - d)*lc)];
MC6a = [ac6a*exp(-i*(b1 + d)*lc) bc6a*exp(-i*(b1 + d)*lc); -conj(bc6a)*exp(-i*(b2 - d)*lc) conj(ac6a)*exp(-i*(b2 - d)*lc)];
T6 = MC6*MP6;
T6a = MC6a*MP6;

MP7 = [exp(-i*b1*lpb) 0; 0 exp(-i*b2*lpb)];
K7 = S*cos (deltaj*13/16*L);
K7a = S*cos (deltaj*13/16*L)*(1 + 0.5*cos (2*pi*(13/16 - 0.5)));
ac7 = cos (sqrt (K7^2 + d^2)*lc) + i*d/sqrt (K7^2 + d^2)*sin (sqrt (K7^2 + d^2)*lc);
bc7 = -i*K7/sqrt (K7^2 + d^2)*sin (sqrt (K7^2 + d^2)*lc);
ac7a = cos (sqrt (K7a^2 + d^2)*lc) + i*d/sqrt (K7a^2 + d^2)*sin (sqrt (K7a^2 + d^2)*lc);
bc7a = -i*K7a/sqrt (K7a^2 + d^2)*sin (sqrt (K7a^2 + d^2)*lc);
MC7 = [ac7*exp(-i*(b1 + d)*lc) bc7*exp(-i*(b1 + d)*lc); -conj(bc7)*exp(-i*(b2 - d)*lc) conj(ac7)*exp(-i*(b2 - d)*lc)];
MC7a = [ac7a*exp(-i*(b1 + d)*lc) bc7a*exp(-i*(b1 + d)*lc); -conj(bc7a)*exp(-i*(b2 - d)*lc) conj(ac7a)*exp(-i*(b2 - d)*lc)];
T7 = MC7*MP7;
T7a = MC7a*MP7;

MP8 = [exp(-i*b1*lpa) 0; 0 exp(-i*b2*lpa)];
K8 = S*sin (deltaj*15/16*L);
K8a = S*sin (deltaj*15/16*L)*(1 + 0.5*cos (2*pi*(15/16 - 0.5)));
ac8 = cos (sqrt (K8^2 + d^2)*lc) + i*d/sqrt (K8^2 + d^2)*sin (sqrt (K8^2 + d^2)*lc);
bc8 = -i*K8/sqrt (K8^2 + d^2)*sin (sqrt (K8^2 + d^2)*lc);
ac8a = cos (sqrt (K8a^2 + d^2)*lc) + i*d/sqrt (K8a^2 + d^2)*sin (sqrt (K8a^2 + d^2)*lc);
bc8a = -i*K8a/sqrt (K8a^2 + d^2)*sin (sqrt (K8a^2 + d^2)*lc);
MC8 = [ac8*exp(-i*(b1 + d)*lc) bc8*exp(-i*(b1 + d)*lc); -conj(bc8)*exp(-i*(b2 - d)*lc) conj(ac8)*exp(-i*(b2 - d)*lc)];
MC8a = [ac8a*exp(-i*(b1 + d)*lc) bc8a*exp(-i*(b1 + d)*lc); -conj(bc8a)*exp(-i*(b2 - d)*lc) conj(ac8a)*exp(-i*(b2 - d)*lc)];
T8 = MC8*MP8;
T8a = MC8a*MP8;

C = T8*T7*T6*T5*T4*T3*T2*MC1;
PCE (int) = 10*log10 ((abs (C (1, 2)))^2);
Ca = T8a*T7a*T6a*T5a*T4a*T3a*T2a*MC1a;
PCEa (int) = 10*log10 ((abs (Ca (1, 2)))^2);

dvec (int) = delta;
end;

plot (dvec*L, PCE, ' : k');
hold on
plot (dvec*L, PCEa, ' k');
hold on
plot (dvec*L, -3, ' k', dvec*L, -10, ' k', dvec*L, -20, ' k');
hold off
set (gcf, ' Color', ' white');
xlabel ('\DeltaL');
ylabel (' PCE (dB)');
xlim ([-30, 30]);
ylim ([-30, 0]);
