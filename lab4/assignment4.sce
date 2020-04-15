clear()

/* --- Load data --- :
- signal_with_noise
- sampling frequency Fs
- signal_filtered
*/
DATA = 'data/';
load(DATA + "signal_with_noise_and_filtered.sod")

disp(Fs)
//playsnd(signal_with_noise,Fs)
//playsnd(signal_filtered,Fs)

PLOTS = 'plots/'

/* --- Add Library function --- */
exec('cshift.sci');


/* --- Signal with noise (TIME domain) --- */
s = signal_with_noise
// Take the first channel
s = s(1, :)

// Plot 
figure(0)
plot(s) 
xlabel("Time, n", 'fontsize', 2)
ylabel("Amplitude", 'fontsize', 2)
title("Signal with noise in time domain", 'fontsize', 3)
xs2png(gcf(), PLOTS + "signal_with_noise_time.png")
close() // close immediately after saving


/* --- Sinal with noice (FREQ domain - spectrum) --- */
s_len = length(s)
frequencies = (0:s_len-1)/s_len * Fs; // calculate frequences
disp(s_len)
mprintf('Frequencies length %d',length(frequencies))
figure(1)
plot2d("nl", frequencies, abs(fft(s)), color("blue")) // plot
xlabel("Frequency component, n", 'fontsize', 2)
ylabel("Freq Amplitude", 'fontsize', 2)
title("Signal with noise in frequency domain", 'fontsize', 3)
xs2png(gcf(), PLOTS + "signal_with_noise_freq.png")
close() // close immediately after saving


/* 
--- Ideal Low Pass filter ---

N:          integer length of FIR filter 
cutoff:     fraction of Fs, at which frequencies are stopped
stop_value: the value for frequencies in the stop band (after cutoff frequency)
return:     frequency representation of an ideal
            low pass FIR filter of length N+1 if N is even or N if N is odd
*/
function H = ideal_lowpass(N, cutoff, stop_value)
    N = (N - modulo(N,2)) / 2
    cutoff = floor(2 * N * cutoff)
    H = ones(1, N) * stop_value
    H(1,1:cutoff) = 1.
    H = [1. H flipdim(H, 2)] // need to make N odd // <---- line 14
endfunction

с
/* 
--- Ideal High Pass filter ---

N:          integer length of FIR filter 
cutoff:     fraction of Fs, at which frequencies are stopped
stop_value: the value for frequencies in the stop band (before cutoff frequency)
return:     frequency representation of an ideal
            low pass FIR filter of length N+1 if N is even or N if N is odd
*/
function H = ideal_highpass(N, cutoff, stop_value)
    N = (N - modulo(N,2)) / 2
    cutoff = floor(2 * N * cutoff)
    H = ones(1, N) * stop_value
    H(1,cutoff:N) = 1.
    H = [0. H flipdim(H, 2)] // need to make N odd // <---- line 14
endfunction


// Plot ideal lowpass freq response
// calculate lowpass
 
H_l = ideal_lowpass(2000, 0.15, 0.); // Filter will have length 257
hl_len = length(H_l)
frequencies = (0:hl_len-1)/hl_len * Fs; // calculate frequencies
figure(2)
plot2d("nn", frequencies, H_l, color("blue")) // plot
xlabel("Frequency, Hz", 'fontsize', 2)
ylabel("Freq amplitude", 'fontsize', 2)
title("Frequency response of ideal low-pass filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "ideal_lowpass_freq.png")
close() // close immediately after saving


// Plot ideal highpass freq response
// calculate highpass
 
H_h = ideal_highpass(2000, 0.001, 0.); // Filter will have length 257
hh_len = length(H_h)
frequencies = (0:hh_len-1)/hh_len * Fs; // calculate frequencies
figure(3)
plot2d("nn", frequencies, H_l, color("blue")) // plot
xlabel("Frequency, Hz", 'fontsize', 2)
ylabel("Freq amplitude", 'fontsize', 2)
title("Frequency response of ideal high-pass filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "ideal_highpass_freq.png")
close() // close immediately after saving


/* 
--- Compute lowpass impulse response --- 
     
    project into temporal domain
    imaginary part should be close to 0
*/
h_l = real(ifft(H_l))
figure(4)
plot2d('nn', 0:length(h_l)-1, h_l, color("blue"))
xlabel("Time, n", 'fontsize', 2)
ylabel("Amplitude", 'fontsize', 2)
title("Impulse response of ideal low-pass filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "ideal_lowpass_time.png")
close()

/* 
--- Compute highpass impulse response --- 

    project into temporal domain
    imaginary part should be close to 0
*/
h_h = real(ifft(H_h))
figure(5)
plot2d('nn', 0:length(h_h)-1, h_h, color("blue"))
xlabel("Time, n", 'fontsize', 2)
ylabel("Amplitude", 'fontsize', 2)
title("Impulse response of ideal high-pass filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "ideal_lowpass_time.png")
close()



/*
--- Shift ---
*/
h_l = cshift(h_l, [0 (hl_len - modulo(hl_len,2)) / 2])
h_h = cshift(h_h, [0 (hh_len - modulo(hh_len,2)) / 2])

/*
--- Window Function on lowpass---
*/

h_l = h_l .* window('kr', length(h_l), 8)
figure(6)
plot2d('nn', 0:length(h_l)-1, h_l, color("blue"))
xlabel("Time, n", 'fontsize', 2)
ylabel("Amplitude", 'fontsize', 2)
title("Windowed impulse response of ideal low-pass filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "windowed_ideal_lowpass_time.png")
//close()

/*
--- Window Function on highpass---
*/
h_h = h_h .* window('kr', length(h_h), 8)
figure(7)
plot2d('nn', 0:length(h_h)-1, h_h, color("blue"))
xlabel("Time, n", 'fontsize', 2)
ylabel("Amplitude", 'fontsize', 2)
title("Windowed impulse response of ideal high-pass filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "windowed_ideal_highpass_time.png")
//close()


figure(8)
plot2d("nl", frequencies,  abs(fft(h_l)), color("blue"))
xlabel("Time, n", 'fontsize', 2)
ylabel("Amplitude", 'fontsize', 2)
title("Frequency response of the final lowpass FIR filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "freq_resp_FIR_Filter.png")
//close()

figure(9)
plot2d("nl", frequencies,  abs(fft(h_h)), color("blue"))
xlabel("Time, n", 'fontsize', 2)
ylabel("Amplitude", 'fontsize', 2)
title("Frequency response of the final highpass FIR filter", 'fontsize', 3)
xs2png(gcf(), PLOTS + "freq_resp_FIR_Filter.png")
//close()


con = convol(s, h_l)
con = convol(con, h_h)

c_len = length(con)
c_frequencies = (0:c_len-1)/c_len * Fs;
s_frequencies = (0:s_len-1)/s_len * Fs;


savewave(PLOTS + "filtered", con, Fs)
figure(10)
clf()
plot2d("nl", s_frequencies, abs(fft(s)), color("blue"))
plot2d("nl", c_frequencies, fft(con), color("red"))
xlabel("Frequency, Hz", 'fontsize', 2)
ylabel("Frquency amplitude", 'fontsize', 2)
title("Result of filtering", 'fontsize', 3)
legend(['Original';'Filtered'])
xs2png(gcf(), PLOTS + "filtering_result.png")

