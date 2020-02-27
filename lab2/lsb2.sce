function plot_sig(fig,data)
    f = figure(fig)
    clf
    plot(data, '--o')
    gca.data_bounds = [0,-0.5;size(data)(1),0.5]
    xlabel('Samples')
    ylabel('Amplitude')
endfunction

function run()

b = chdir('.')
exec('ADC.sce')

// configs
n = 5
fs = 50000
quant_levels = -0.5:0.0005:0.1

// recorded data
data = ADC(n, quant_levels, fs)

plot_sig(1,data)

// shifted data
data = data + 0.1

plot_sig(2,data)

// sin noise configs
sin_freq = 165
sin_ampl = 0.1
step_size = sin_freq*(2*%pi)/fs;
samples = [1:size(data)(1)]*step_size;

// sin signal subtraction
sin_sig = sin_ampl*sin(samples)
data = data - sin_sig'

plot_sig(3,data)

playsnd(data, fs)

endfunction

run()
