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
  %   gamma(gamma_coeffs) = sum_{i=1}^N gamma_coeffs(i) * b_spline( (N-2)*t + i - 3 )
  %   phase(gamma) = imag( gamma )
  %   log_ampl(gamma) = real( gamma )
  %   ampl(gamma) = exp( log_ampl(gamma) )
  %   samples(gamma) = ampl(ampl) * cos( phase(gamma) )
  %
  % The cost of this fit is given by
  %
  %   cost(gamma_coeffs) = sumsq( samples(gamma(gamma_coeffs) - signal ) + punishing_term( gamma_coeffs )
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
  
  % 1. Find an initial guess
  % 1.1. Make an initial guess of ampl' and phase'
  % 1.2. Calculate a matching vector gamma'
  % 1.3. Find gamma_coeffs, such that gamma(gamma_coeffs) and
  %      gamma'(ampl',phase') have minimal l2 distance. 
  % 2. Distribute a swarm around the found gamma_coeffs. 
  % 3. Subtract a constant from each member of the swarm, such that 
  %    each swarm member has zero average. (This is for better numerical stability).
  % 4. Find optimal values (modulo a constant offset) for gamma_coeffs
  % 4.1. Let offset(gamma_coeffs,signal) be the function that calculates 
  %      the optimal offset. 
  % 4.2. Let cost_off(gamma_coeffs) be the function that calculates 
  %      the cost of gamma_coeffs changed by an optimal offset. 
  %      cost_off(gamma_coeffs) = 
  %          cost(gamma_coeffs+offset(gamma_coeffs))
  % 4.3. best_fit_gamma_coeffs_off = DE( swarm, cost_off )
  % 5. Calculate the results
  % 5.1. best_fit_gamma_coeffs = best_fit_gamma_coeffs_off + offset(best_fit_gamma_coeffs_off)
  % 5.2. best_fit_
  
  % make an initial guess in terms of amplitude and phase function
  % find a vector gamma_coeffs yields similiar amplitude and phase
  % distribute swarm around initial guess
  % optimize swarm with differential evolution
  % return the amplitude and phase functions of the best fit
endfunction