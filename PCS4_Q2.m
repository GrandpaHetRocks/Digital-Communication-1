%raised cosine
tb=2;  
t=-10:1/(100*tb):10; %time vector

r=0.5; %roll off factor=0.5
p=(sin(pi*t/tb)./(pi*t/tb)).*(cos((pi*t/tb)*r)./(1-(4*r*r.*t.*t/(tb*tb)))); %p(t)
plot(t,p)

hold on

r=0.75; %roll off factor=0.75
p=(sin(pi*t/tb)./(pi*t/tb)).*(cos((pi*t/tb)*r)./(1-(4*r*r.*t.*t/(tb*tb)))); %p(t)
plot(t,p)

r=0; %roll off factor=0 [sinc function]
p=(sin(pi*t/tb)./(pi*t/tb)).*(cos((pi*t/tb)*r)./(1-(4*r*r.*t.*t/(tb*tb)))); %p(t)
plot(t,p)
legend("r=0.5","r=0.75","r=0->sinc(pi*t/tb)")
title("time domain")


%frequency response of raised cosine
figure
f=-5/tb:1/(100*tb):5/tb;  %frequency range

r=0.5;
P=[]; %frequency response to be calculated
for i=1:length(f) %calculating frequency response at each point
    if(abs(f(i))<=((1-r)/(2*tb)))
        P=[P 1];
    else if(abs(f(i))<=((1+r)/(2*tb)) && abs(f(i))>((1-r)/(2*tb)))
            P=[P (0.5*(1+cos(pi*tb/r*(abs(f(i))-((1-r)/(2*tb))))))];  
    else
        P=[P 0];
        end
    end
end
hold on
plot(f,P)


r=0.75; 
P=[]; %frequency response to be calculated
for i=1:length(f) %calculating frequency response at each point
    if(abs(f(i))<=((1-r)/(2*tb)))
        P=[P 1];
    else if(abs(f(i))<=((1+r)/(2*tb)) && abs(f(i))>((1-r)/(2*tb)))
            P=[P (0.5*(1+cos(pi*tb/r*(abs(f(i))-((1-r)/(2*tb))))))];  
    else
        P=[P 0];
        end
    end
end

plot(f,P)

r=0;   %at r=0 raised cosine is basically a sinc function
P=[];  %frequency response to be calculated
for i=1:length(f)  %calculating frequency response at each point
    if(abs(f(i))<=((1-r)/(2*tb)))
        P=[P 1];
    else if(abs(f(i))<=((1+r)/(2*tb)) && abs(f(i))>((1-r)/(2*tb)))
            P=[P (0.5*(1+cos(pi*tb/r*(abs(f(i))-((1-r)/(2*tb))))))];  
    else
        P=[P 0];
        end
    end
end
hold on
plot(f,P)
legend("r=0.5","r=0.75","r=0->sinc(pi*t/tb)")
title("frequency domain")

%%we can see that the sinc function has a lesser bandwith than the raised
%%cosine. we know that bandwidth= ((1+r)/2)*(1/tb). so it is directly
%%proportional to r. more the r, more is the bandwidth.

%%i chose to not go with the fft function in matlab because of many 
%%failed attempts to get the plot right after using the same. most of
%%the times i was getting NaN as an output, and on setting a specified
%%range, the output came out wrong. so i manually implemented the frequency
%%domain representation of the signals.
            


