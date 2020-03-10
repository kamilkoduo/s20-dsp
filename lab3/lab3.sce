funcprot(0)
clear all;
x = [1,2,2,1];
h = [1,1,1,1];

disp(x,'x:')
disp(h,'h:')

function y=conv_direct(x,h)
    m = length(x);
    n = length(h);
    for i = 1: n + m -1
        conv_sum = 0;
        for j = 1:i
            if ( ( (i-j+1)<=n ) & (j<=m) )
                conv_sum = conv_sum + x(j)*h(i-j+1);
            end;
            y(i) = conv_sum;
        end;
    end;
endfunction;

function y=conv_freqdom(x,h)
    m = length(x);
    n = length(h);
    
    N = n + m-1;
    x = [x zeros(1,N-m)];
    h = [h zeros(1,N-n)];
    f1 = fft(x)
    f2 = fft(h)
    f3 = f1 .* f2 ;
    f4 = ifft(f3)
    y = f4
endfunction;


tic()
y = conv_direct(x,h)
t = toc()
mprintf('Convolution x*h using Direct Formula Method took time = %f ms',1000*t)
disp(y')

tic()
y = convol(x,h)
t=toc()
mprintf('Convolution x*h using Built-In ""convol"" took time = %f ms',1000*t)
disp(y)

tic()
y = conv_freqdom(x,h)
t=toc()
mprintf('Convolution x*h using Frequency Domain multiplication took time = %f ms',1000*t)
disp(y)

