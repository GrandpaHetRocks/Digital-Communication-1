T=10;
t2=-5:0.05:5;
t=-5:1:5; %sampled to get 10 discrete values
mu=0; %mean
sig=0.5^0.5; %standard deviation
n=sig*randn(1,11)+mu; 
%gaussian noise randomly generated at 11 points
A=2; %amplitude
s0_cont=A*cos(pi*t2/T);  %continuous s0 (unsampled)
s1_cont=A*cos(2*pi*t2/T); %continuous s1 (unsampled)
s0=A*cos(pi*t/T);  %sampled s0
s1=A*cos(2*pi*t/T); %sampled s1
r0=s0+n;
r1=s1+n;
subplot(2,1,1)
hold on
stem(t,r0,"linewidth",1.5)
plot(t,r0,"linewidth",1.5)
plot(t2,s0_cont,"r--")
legend("discrete r0","r0","s0",'location', 'westoutside');
title("signal s0 after adding noise at 11 random points")
subplot(2,1,2)
hold on
plot(t,r1,"linewidth",1.5)
stem(t,r1,"linewidth",1.5)
plot(t2,s1_cont,"b--")
% plot(t1,n)
legend("discrete r1","r1","s1",'location', 'westoutside');
title("signal s1 after adding noise at 11 random points")

%randn is a matlab function that generates gaussian noise with sigma=1 and
%mean=0. I've tweaked it a bit to get the desired variance and mean.

%using awgn routine in matlab variance=0.5 signal is sampled
figure
snr1=snr(s0,n); %find snr of signal to gaussian noise(var=0.5)
snr2=snr(s1,n); %find snr of signal to gaussian noise(var=0.5)
r0_new1=awgn(s0,snr1);
r1_new1=awgn(s1,snr2);
plot(t,r0_new1,"linewidth",1.5)
hold on
plot(t,r1_new1,"linewidth",1.5)
plot(t2,s0_cont,"--")
plot(t2,s1_cont,"--")
legend("r0","r1","s0","s1",'location', 'westoutside');
title("using awgn variance=0.5 sampled signal")

%using awgn routine in matlab variance=2 signal is sampled
sig=2^0.5;
n=sig*randn(1,11)+mu; %gaussian noise with variance=2
figure
snr1=snr(s0,n); %find snr of signal to gaussian noise(var=2)
snr2=snr(s1,n); %find snr of signal to gaussian noise(var=2)
r0_new=awgn(s0,snr1); 
r1_new=awgn(s1,snr2);
plot(t,r0_new,"linewidth",1.5)
hold on
plot(t,r1_new,"linewidth",1.5)
plot(t2,s0_cont,"--")
plot(t2,s1_cont,"--")
legend("r0","r1","s0","s1",'location', 'westoutside');
title("using awgn variance=2 sampled signal")

%using awgn routine in matlab variance=0.5 signal is not sampled
sig=0.5^0.5;
n=sig*randn(1,length(t2))+mu;  %gaussian noise with variance=0.5
figure
snr1=snr(s0_cont,n); %find snr of signal to gaussian noise(var=0.5)
snr2=snr(s1_cont,n); %find snr of signal to gaussian noise(var=0.5)
r0_new2=awgn(s0_cont,snr1); 
r1_new2=awgn(s1_cont,snr2);
plot(t2,r0_new2,"linewidth",1.5)
hold on
plot(t2,r1_new2,"linewidth",1.5)
plot(t2,s0_cont,"--")
plot(t2,s1_cont,"--")
legend("r0","r1","s0","s1",'location', 'westoutside');
title("using awgn variance=0.5 on the unsampled signal")

%using awgn routine in matlab variance=2 signal is not sampled
sig=2^0.5;
n=sig*randn(1,length(t2))+mu;  %gaussian noise with variance=2
figure
snr1=snr(s0_cont,n); %find snr of signal to gaussian noise(var=2)
snr2=snr(s1_cont,n); %find snr of signal to gaussian noise(var=2)
r0_new3=awgn(s0_cont,snr1); 
r1_new3=awgn(s1_cont,snr2);
plot(t2,r0_new3,"linewidth",1.5)
hold on
plot(t2,r1_new3,"linewidth",1.5)
plot(t2,s0_cont,"--")
plot(t2,s1_cont,"--")
legend("r0","r1","s0","s1",'location', 'westoutside');
title("using awgn variance=2 on the unsampled signal")

% we can see that on increasing variance, the distortion/noise increases
% it wasn't mentioned in the question whether to use awgn on the sampled
% signal or the unsampled signal so i plotted both.

eyediagram(r0_new,100)
title("eyediagram for r0 variance=2 sampled")
eyediagram(r1_new,100)
title("eyediagram for r1 variance=2 sampled")
eyediagram(r0_new1,100)
title("eyediagram for r0 variance=0.5 sampled")
eyediagram(r1_new1,100)
title("eyediagram for r1 variance=0.5 sampled")
eyediagram(r0_new2,100)
title("eyediagram for r0 variance=0.5 unsampled")
eyediagram(r1_new2,100)
title("eyediagram for r1 variance=0.5 unsampled")
eyediagram(r0_new3,100)
title("eyediagram for r0 variance=2 unsampled")
eyediagram(r1_new3,100)
title("eyediagram for r1 variance=2 unsampled")

%similar to noise, the eye diagram distorts more as the variance increases

figure
h0=A*cos(pi*(t-T)/T);
h1=A*cos(2*pi*(t-T)/T);
r0_conv=conv(h0,r0_new1);
r1_conv=conv(h1,r1_new1);
stem(r0_conv)
hold on
stem(r1_conv)
legend("r0'","r1'",'location', 'westoutside');
title("r' for variance=0.5 signal is sampled")

figure
r0_conv=conv(h0,r0_new);
r1_conv=conv(h1,r1_new);
stem(r0_conv)
hold on
stem(r1_conv)
legend("r0'","r1'",'location', 'westoutside');
title("r' for variance=2 signal is sampled")

h0=A*cos(pi*(t2-T)/T);
h1=A*cos(2*pi*(t2-T)/T);
figure
r0_conv=conv(h0,r0_new2);
r1_conv=conv(h1,r1_new2);
stem(r0_conv)
hold on
stem(r1_conv)
legend("r0'","r1'",'location', 'westoutside');
title("r' for variance=0.5 signal is unsampled")

figure
r0_conv=conv(h0,r0_new3);
r1_conv=conv(h1,r1_new3);
stem(r0_conv)
hold on
stem(r1_conv)
legend("r0'","r1'",'location', 'westoutside');
title("r' for variance=2 signal is unsampled")

% here we know that increasing varaince will increase the noise, hence on
% convolution the amplitude increases as the variance increases. we can see
% that the plot with variance=2 has higher amplitude than one with
% variance=0.5.

