image=imread("./assign1.jpg");
x = reshape((dec2bin(image,8)-'0').',1,[]);
stairs(x(1:24))
title("Time domain representaion of image (24 bits @ 1bit/s)") %1bit/sec
figure
[m,n]=size(image);
imshow(image,[0,255])
title("Input Image")
choice=input('bpsk or qpsk or qam ','s')  %bpsk or qpsk or qam


reversed=x;

if(choice=="bpsk")
bpsk=pskmod(reversed,2,pi);
scatterplot(bpsk)
title("Constellation Plot: BPSK")
t=0:0.001:1-0.001; %1 symbol/sec
time_domain=[];
for i=bpsk(1:4)
    s_i=real(i);
    s_q=imag(i);
    wave=s_i*(cos(2*pi*10*t))-s_q*sin(2*pi*10*t); 
    time_domain=[time_domain wave];
end
figure
plot([0:0.001:4-0.001],time_domain)
title("Time Domain: BPSK")

points=1024;
Fs=1000;
Y=fft(time_domain,points); 

figure

f1=(-points/2:points/2-1)*Fs/points;
Y=fftshift(Y);
stem(f1,abs(Y))
title("Frequency Domain: BPSK")

end

if(choice=="qpsk")
signal_qpsk=[];
for i=1:2:length(reversed) %taking 2 bits at a time for modulation
    signal_qpsk=[signal_qpsk 2*reversed(i)+reversed(i+1)];
end


qpsk=pskmod(signal_qpsk,4,pi/2);

t=0:0.001:1-0.001;
time_domain=[];
for i=qpsk(1:4) %4 symbols ie 8 bits modulated
    s_i=real(i);
    s_q=imag(i);
    wave=s_i*(cos(2*pi*10*t))-s_q*sin(2*pi*10*t);
    time_domain=[time_domain wave];
end
figure
plot([0:0.001:4-0.001],time_domain)
title("Time Domain QPSK")

points=1024;
Fs=1000;
Y=fft(time_domain,points); 

figure

f1=(-points/2:points/2-1)*Fs/points;
Y=fftshift(Y);
stem(f1,abs(Y))
title("Frequency Domain: QPSK")

scatterplot(qpsk)
title("Constellation Plot: QPSK")
end

if(choice=="qam")
signal_qam=[];
for i=1:4:length(reversed) %taking 4 bits for modulation at a time
    signal_qam=[signal_qam 8*reversed(i)+4*reversed(i+1)+2*reversed(i+2)+reversed(i+3)];
end

qam=qammod(signal_qam,16);

t=0:0.001:1-0.001;
time_domain=[];
for i=qam(1:4)
    s_i=real(i);
    s_q=imag(i);
    wave=s_i*(cos(2*pi*10*t))-s_q*sin(2*pi*10*t);
    time_domain=[time_domain wave];
end
figure
plot([0:0.001:4-0.001],time_domain)
title("Time Domain: QAM")

points=1024;
Fs=1000;
Y=fft(time_domain,points); 

figure

f1=(-points/2:points/2-1)*Fs/points;
Y=fftshift(Y);
stem(f1,abs(Y))
title("Frequency Domain: QAM")
ylim([0 500])

scatterplot(qam)
title("Constellation Plot: 16 QAM")
end

%demodulating

if(choice=="bpsk")
high_snr_bpsk=awgn(bpsk,40);
low_snr_bpsk=awgn(bpsk,1);

bpsk_demod=pskdemod(high_snr_bpsk,2,pi);
stream=[];
for i=bpsk_demod
    stream=[stream de2bi(i,1,'left-msb')];
end

figure
% stream=bit(stream);
s = num2cell(reshape(stream,8,[])',2);
b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
out=reshape(b,m,n);  
imshow(out,[0,255])
title("Reconstructed from BPSK (HIGH SNR)");


bpsk_demod=pskdemod(low_snr_bpsk,2,pi);
stream=[];
for i=bpsk_demod
    stream=[stream de2bi(i,1,'left-msb')];
end

figure
% stream=bit(stream);
s = num2cell(reshape(stream,8,[])',2);
b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
out=reshape(b,m,n);  
imshow(out,[0,255])
title("Reconstructed from BPSK (LOW SNR)");
end


if(choice=="qpsk")
high_snr_qpsk=awgn(qpsk,40);
low_snr_qpsk=awgn(qpsk,1);


qpsk_demod=pskdemod(high_snr_qpsk,4,pi/2);
stream=[];
for i=qpsk_demod
    stream=[stream de2bi(i,2,'left-msb')];
end

figure
% stream=bit(stream);
s = num2cell(reshape(stream,8,[])',2);
b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
out=reshape(b,m,n);  
imshow(out,[0,255])
title("Reconstructed from QPSK (HIGH SNR)");

qpsk_demod=pskdemod(low_snr_qpsk,4,pi/2);
stream=[];
for i=qpsk_demod
    stream=[stream de2bi(i,2,'left-msb')];
end

figure
% stream=bit(stream);
s = num2cell(reshape(stream,8,[])',2);
b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
out=reshape(b,m,n);  
imshow(out,[0,255])
title("Reconstructed from QPSK (LOW SNR)");
end

if(choice=="qam")
high_snr_qam=awgn(qam,40);
low_snr_qam=awgn(qam,1);

qam_demod=qamdemod(high_snr_qam,16);
stream=[];
for i=qam_demod
    stream=[stream de2bi(i,4,'left-msb')];
end

figure
% stream=bit(stream);
s = num2cell(reshape(stream,8,[])',2);
b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
out=reshape(b,m,n);  
imshow(out,[0,255])
title("Reconstructed from 16 QAM (HIGH SNR)");


qam_demod=qamdemod(low_snr_qam,16);
stream=[];
for i=qam_demod
    stream=[stream de2bi(i,4,'left-msb')];
end

figure
% stream=bit(stream);
s = num2cell(reshape(stream,8,[])',2);
b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
out=reshape(b,m,n);  
imshow(out,[0,255])
title("Reconstructed from 16 QAM (LOW SNR)");
end

%image reconstruction code picked from the internet








