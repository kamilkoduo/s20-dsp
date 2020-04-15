clf()

DATA = 'data/'

[irc, F_irc, _] = wavread(DATA+ '1a_marble_hall.wav')
irc = irc(1, :)

[sig, F_sig, _] = wavread(DATA+'7cef8230.wav')
sig = sig(1, :)


irc_inv = 1./fft(irc)
irc_inv(1, find(isinf(irc_inv))) = 0.
h = ifft(irc_inv)

h_len = length(h)
h = cshift(h, [0 (h_len-modulo(h_len, 2))/2])
h = h .* window('kr', length(h), 8)


con = convol(sig, irc)//clear with irc
filtered = convol(con, h) //echo with reversed echo

con = con ./60

figure(0)
plot2d("nn", 1:length(con), con, color("blue"))

xlabel("Time, t", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title('Convolution result', 'fontsize', 4)

figure(1)
plot2d("nn", 1:length(filtered), filtered, color("red"))

xlabel("Time, t", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title('Filter result', 'fontsize', 4)

//playsnd(filtered,F_sig)
savewave(DATA + 'filered2.wav ',filtered,F_sig)

kronecker = convol(irc,h)

figure(2)
plot2d("nn", 1:length(kronecker), kronecker, color("red"))
title('Kronecker', 'fontsize', 4)
xlabel("Time, t", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)

figure(3)
plot2d("nn", 1:length(h), h, color("blue"))

xlabel("Time, t", 'fontsize', 3)
ylabel("Amplitude", 'fontsize', 3)
title('Filter for cancelling echo', 'fontsize', 4)
/**/
