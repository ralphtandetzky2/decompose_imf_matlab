function [best_fit_samples,best_fit_ampl,best_fit_phase,best_fit_gamma_coeffs] ...
  = decompose_imf( signal, N, swarm_size, std_dev_log_ampl, std_dev_phase, ...
                   n_iters, differential_weight, crossover_probability )
% DECOMPOSE_IMF Find an intrinsic mode function approximating a signal.
%
% signal                = input signal
% N                     = The number of complex parameters to fit. 
%                         This must be an integer in the interval 
%                         [3,length(signal)]. 
% swarm_size            = The size of the swarm to be used for 
%                         differential evolution.
% std_dev_log_ampl      = The standard deviation of the logarithm of 
%                         the amplitide during swarm generation.
% std_dev_phase         = The standard deviation of the phase during 
%                         swarm generation.
% n_iters               = number of iterations during global optimization.
% differential_weight   = differential weight used in 
%                         differential evolution algorithm.
% crossover_probability = crossover probability used in 
%                         differential evolution algorithm.
%
% OUTPUT PARAMETERS
%
% best_fit_samples      = The samples of the first IMF fitted to the 
%                         signal. The length of this vector is equal 
%                         to the length of the signal vector.
% best_fit_ampl         = The amplitude vector of the best_fit_samples.
% best_fit_phase        = The phase vector of the best_fit_samples.
%
% By definition best_fit_samples equals best_fit_ampl*cos(best_fit_phase). 
%
% See also DIFFERENTIAL_EVOLUTION
  
  % We optimize gamma_coeffs which is a complex valued vector of length N. 
  % The resulting samples, amplitude and phase function can be calculated
  % in the following manner:
  %
  %   t                   = ( (1:length(signal)) - 0.5 ) / length(signal)
  %   gamma(gamma_coeffs) = sum_{i=1}^N gamma_coeffs(i) * b_spline( (N-2)*t - i + 3 )
  %   phase(gamma)        = imag( gamma )
  %   log_ampl(gamma)     = real( gamma )
  %   ampl(gamma)         = exp( log_ampl(gamma) )
  %   samples(gamma)      = ampl(ampl) * cos( phase(gamma) )
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
  % Hence samples are changed to 
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
  %   offset = log(A-iB)
  %
  % Note also, that adding such an offset does not change the punishing_term. 
  % Hence, an optimal offset can be easily found and added in every 
  % optimization step. 
  %
  % OUTLINE OF THE ALGORITHM
  %
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
  % 5.2. best_fit_gamma   = gamma  (best_fit_gamma_coeffs)
  % 5.3. best_fit_samples = samples(best_fit_gamma)
  % 5.4. best_fit_ampl    = ampl   (best_fit_gamma)
  % 5.5. best_fit_phase   = phase  (best_fit_gamma)
  
  % IMPLEMENTATION

  if ( N < 3 )
    error('The number of complex parameters to be fitted must be at least 3.');
  end
  
  % 1. Find an initial guess
  % 1.1. Make an initial guess of ampl' and phase'
  
  if ( N == 3 )
    [ampl,phase] = dimf_initial_guess(signal);
  else
    [~,ampl,phase] = ...
        decompose_imf( signal, floor(N/2)+1, ceil(swarm_size/2), ...
            std_dev_log_ampl*2, std_dev_phase*2, ...
            n_iters/2, differential_weight, crossover_probability );
  end
  
  % 1.2. Calculate a matching vector gamma'
  
  gamma = log(ampl) + 1i*phase;
  
  % 1.3. Find gamma_coeffs, such that gamma(gamma_coeffs) and
  %      gamma'(ampl',phase') have minimal l2 distance. 

  t = ( (1:length(signal)) - 0.5 ) / length(signal);
  gamma_base = zeros(N,length(t));
  for iter = 1:N
    gamma_base(iter,:) = dimf_b_spline( (N-2)*t - iter + 3 );
  end;
  gamma_coeffs = gamma / gamma_base;
  
  % 2. Distribute a swarm around the found gamma_coeffs. 
  
  swarm = ones(swarm_size,1) * gamma_coeffs;
  swarm = swarm  + normrnd(0,std_dev_log_ampl,size(swarm)) ...
              + 1i*normrnd(0,std_dev_phase   ,size(swarm));
  
  % 3. Subtract a constant from each member of the swarm, such that 
  %    each swarm member has zero average. (This is for better numerical stability).
  
  swarm = swarm - mean(swarm,2) * ones(1,size(swarm,2));
  
  % 4. Find optimal values (modulo a constant offset) for gamma_coeffs
  % 4.1. Let offset(gamma_coeffs,signal) be the function that calculates 
  %      the optimal offset. 
  % 4.2. Let cost_off(gamma_coeffs) be the function that calculates 
  %      the cost of gamma_coeffs changed by an optimal offset. 
  %      cost_off(gamma_coeffs) = 
  %          cost(gamma_coeffs+offset(gamma_coeffs))
  
  function retval = dimf_cost_off_functor( x )
    retval = dimf_cost_off( x, signal );
  end
  
  % 4.3. best_fit_gamma_coeffs_off = DE( swarm, cost_off )
  
  best_fit_gamma_coeffs_off = differential_evolution( ...
    swarm, @dimf_cost_off_functor, n_iters, ...
    differential_weight, crossover_probability );
  
  % 5. Calculate the results
  % 5.1. best_fit_gamma_coeffs = best_fit_gamma_coeffs_off + offset(best_fit_gamma_coeffs_off)
  
  best_fit_gamma_coeffs = best_fit_gamma_coeffs_off + ...
      dimf_offset(dimf_gamma(best_fit_gamma_coeffs_off,t),signal);
  
  % 5.2. best_fit_gamma   = gamma  (best_fit_gamma_coeffs)
  
  best_fit_gamma = dimf_gamma( best_fit_gamma_coeffs, t );
  
  % 5.3. best_fit_samples = samples(best_fit_gamma)
  
  best_fit_samples = dimf_samples( best_fit_gamma );
  
  % 5.4. best_fit_ampl    = ampl   (best_fit_gamma)
  
  best_fit_ampl = exp( real( best_fit_gamma ) );
  
  % 5.5. best_fit_phase   = phase  (best_fit_gamma)
  
  best_fit_phase = imag( best_fit_gamma );

end
