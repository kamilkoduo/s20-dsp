clf()

/*= Task 1 =*/

function X = my_dft(x)
    N = length(x)
    X = zeros(0:N-1)
    for k = 1:N
        for n = 1:N
            X(k) = X(k) + x(n) * exp(-2 * %i * %pi * (k-1) * (n-1) / N)
        end
    end
endfunction

function X = my_fft(x)
  N = length(x)
  
  if (N <= 1) then
    X = x
  elseif (modulo(N,2)>0) then
        error('We assert signal length to be a power of 2 | N!=2^M')
  else
    X_even = my_fft(x(1:2:$)) 
    X_odd = my_fft(x(2:2:$))
    
    k = [0:1:N/2-1]
    factor = exp(k .* (-2*%i*%pi/N))
    factor = factor .* X_odd
   
    X = [X_even + factor, X_even - factor]
  end
endfunction


function [t,signal] = gen_signal(Func,f,A,fs,T)
    dt = 1/fs
    t = 0:dt:T-dt
    signal = Func(t .* (2 * %pi * f)) .*A
endfunction


/*= Task 2 =*/

sin_fs = 128

[t1,sin1] = gen_signal(sin,4,1,sin_fs,1)
[t2,sin2] = gen_signal(sin,5,1,sin_fs,1)
[t3,sin3] = gen_signal(sin,11,1,sin_fs,1)

figure(0)
plot(t1,sin1,"b-o-")
plot(t2,sin2,"g-o-")
plot(t3,sin3,"r-o-")
xlabel("Time, t (s)", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title("sin signals sampled at 128 Hz rate", 'fontsize', 4)
legend(['4 Hz';'5 Hz';'11 Hz']);



/*= Task 2. With leakege [sin1 + sin2 + sin3] =*/
sin_sum = sin1+sin2+sin3

figure(1)
plot(t1,sin_sum,"b-o-")
xlabel("Time, t (s)", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title("4 Hz, 5 Hz, 11 Hz sin signals combined, sampled at 128 Hz rate", 'fontsize', 4)


sin_sum_spectrum = abs(fft(sin_sum))
sin_N = length(sin_sum_spectrum)
sin_freqs = (0:sin_N-1) .* (sin_fs/sin_N)

figure(2)
plot(sin_freqs, sin_sum_spectrum,"b-o-")
xlabel("Frequency, f (Hz)", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title("Spectrum of 4 Hz, 5 Hz, 11 Hz sin signals combined, Fs=128 Hz", 'fontsize', 4)


/*= Task 3 =*/

sin_sum_i = resize_matrix(sin_sum, 1, 512) 
dt = 1/sin_fs
t_i = 0:dt:512*dt-dt

figure(3)
plot(t_i,sin_sum_i,"b-o-")
xlabel("Time, t (s)", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title("4 Hz, 5 Hz, 11 Hz sin signals combined, sampled at 128 Hz rate", 'fontsize', 4)


sin_sum_i_spectrum = abs(fft(sin_sum_i))
sin_N_i = length(sin_sum_i_spectrum)
sin_freqs_i = (0:sin_N_i-1) .* (sin_fs/sin_N_i)

figure(4)
plot(sin_freqs, sin_sum_spectrum,"b-o-")
plot(sin_freqs_i, sin_sum_i_spectrum,"r-o-")
xlabel("Frequency, f (Hz)", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title("Spectrum of 4 Hz, 5 Hz, 11 Hz sin signals combined, Fs=128 Hz", 'fontsize', 4)
legend(['original';'interpolated'])


/*= Task 4=*/


fs = 1024
T = 0.5

[t,sig1] = gen_signal(cos,190,0.5,fs,T)
[t,sig2] = gen_signal(cos,10,2,fs,T)

figure(5)
plot(t,sig1,"b-o-")
xlabel("Time, t (s)", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title("190 Hz sin signal", 'fontsize', 4)

figure(6)
plot(t,sig2,"-o-")
xlabel("Time, t (s)", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title("10 Hz sin signal", 'fontsize', 4)
