/**

SIR Covid model
negative binomial observation density

**/
model tn_covid {

  input N
  
  noise e 
  
  obs y

  state S
  state I
  state R
  state x
  state Z

  param k
  param gamma
  param sigma // Noise driver
  param I0
  param R0
  param x0
  param p_rep // reporting rate
  param p_rep0 // reporting rate
  param p_rep1 // reporting rate
  param p_phi // overdispersion in reporting
  
  // prior over parameters
  sub parameter {
    gamma ~ truncated_gaussian(5.058, 1.519, lower=2.228, upper=11.8) // gamma is the period
    k ~ gamma(shape=1.058, scale=2.174)
    sigma ~ uniform(0, 1)  // std dev of noise process
    x0 ~ uniform(-5, 1)  
    I0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.5, 0.15, lower = 0, upper = 1)
    p_rep0 ~ uniform(0.0, 0.5)
    p_rep1 ~ uniform(0.3, 1)
    p_phi ~ uniform(0, 0.5)
  }

  // prior over states
  sub initial {
    S <- N
    R <- R0 * S
    S <- S - R
    I <- exp(I0 + log(S))
    S <- S - I
    x <- x0
    Z <- 0
  }

  sub transition(delta = 1) {
    p_rep <- (t_now > 28 ? p_rep1 : p_rep0)
    Z <- ((t_now) % 7 == 0 ? 0 : Z)
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1e-3, rtoler = 1e-8) {
      dx/dt = sigma * e // geometric Brownian motion
      dS/dt = -(exp(x) * S * I)/N
      dI/dt = (exp(x) * S * I)/N - I/gamma
      dR/dt = I/gamma
      dZ/dt = (exp(x) * S * I)/N
    }
  }

  sub observation {
    y ~ truncated_gaussian(p_rep * Z, sqrt(max(p_rep * (1 - p_rep) * Z + (p_rep * Z * p_phi)**2, 1)), lower=0)
  }

  sub proposal_parameter {
    k ~ gaussian(k, 1.0)
    sigma ~ gaussian(sigma, 1.0)
    gamma ~ gaussian(gamma, 1.0)
    x0 ~ gaussian(x0, 1.0)   
    I0 ~ gaussian(I0, 1.0)
    R0 ~ gaussian(R0, 1.0)
    p_rep ~ truncated_gaussian(p_rep, 0.1, lower = 0, upper = 1.0)
    p_rep0 ~ truncated_gaussian(p_rep0, 0.1, lower = 0, upper = 0.5)
    p_rep1 ~ truncated_gaussian(p_rep1, 0.1, lower = 0.3, upper = 1.0)
    p_phi ~ truncated_gaussian(p_phi, 1, lower = 0, upper = 0.5)   
  }
}
