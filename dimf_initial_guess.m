function [ampl,phase] = dimf_initial_guess( signal )
  N = length(signal);
  dft = fft(signal);
  spec = abs(dft);
  [ampl,freq] = max(spec(2:floor(end/2)));
  ampl = 2 * ampl / N;
  phase = (0:(length(signal)-1)) / N * 2 * pi * freq + angle(dft(freq+1));  
end
