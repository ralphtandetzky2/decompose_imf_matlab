function retval = cost( 

% signal = input signal
% N = The number of complex parameters to fit. 
%     This must be an integer >= 3. 
% best_fit_samples = The samples of the first IMF fitted to the signal. 
%     The length of this vector is equal to the length of the signal vector.
% best_fit_ampl = The amplitude vector of the best_fit_samples.
% best_fit_phase = The phase vector of the best_fit_samples.
% By definition best_fit_samples equals best_fit_ampl * cos(best_fit_phase). 
function [best_fit_samples,best_fit_ampl,best_fit_phase] 
  = decompose_imf( signal, N )
  
  % We optimize gamma_coeffs which is a complex valued vector of length N. 
  % The resulting samples, amplitude and phase function can be calculated
  % in the following manner:
  %
  %   t        = ( (1:length(signal)) - 0.5 ) / length(signal)
  %   gamma    = sum_{i=1}^N gamma_coeffs(i) * b_spline( (N-2)*t + i - 3 )
  %   phase    = imag( gamma )
  %   log_ampl = real( gamma )
  %   ampl     = exp( log_ampl )
  %   samples  = ampl * cos( phase )
  %
  % The cost of this fit is given by
  %
  %   cost = sumsq( samples - signal ) + punishing_term( gamma_coeffs )
  %
  % Here punishing_term(...) is zero, if the boundary condition for 
  % the function gamma hold true. 
  % Otherwise, it is a term that will be quiet large compared to the 
  % square sum to the left. 
  % Therefore, invalid solutions will be extinguished very quickly 
  % during optimization. 
  %
  % Note, that adding a complex constant named offset to all values 
  % gamma_coeffs changes gamma in the same way. 
  % The phase is changed by adding the imaginary part of that offset, 
  % the log_ampl is changed by adding the real part of the offset, 
  % hence the amplitude is changed by multiplying with exp(real(offset)). 
  % Hence samples is changed to 
  %
  %   A*ampl*cos(phase) + B*ampl*sin(phase)
  %
  % where
  %
  %   A =   exp(real(offset)) * cos(imag(offset))
  %   B = - exp(real(offset)) * sin(imag(offset))
  %
  % With a little linear algebra, real A and B can be chosen such that
  %
  %   sumsq( A*ampl*cos(phase) + B*ampl*sin(phase) - signal )
  %
  % is minimized. The value of offset can be chosen to yield any values 
  % of real A and B by putting
  %
  %   offset = sqrt(A*A+B*B) + i*arg(A+iB)
  %
  % Note also, that adding such an offset does not change the punishing_term. 
  % Hence, an optimal offset can be easily found and added in every 
  % optimization step. 
  
  % make an initial guess
  % distribute swarm around initial guess
  % optimize swarm with differential evolution
  % return the amplitude and phase functions of the best fit
endfunction