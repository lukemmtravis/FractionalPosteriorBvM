library(ggplot2)
library(tidyverse)
library(gdata)
library(gridExtra)
library(MASS)
library(latex2exp)
library(gdata)
set.seed(123)
setwd('/Users/lmt15/Documents/phd/Fractional Posterior BvM/codes')
# install.packages("sbde_0.0-14.tar.gz", repos=NULL)
library(sbde)

basis_functions = function(K, L=100){
  # K is the (even) number of basis functions to compute 
  x_range = seq(0,1,length.out = L)
  args = seq(1,K/2,1)*2*pi
  args = as.matrix(x_range)%*%t(as.matrix(args))
  # L is the size of the range discretisation
  cosines = sqrt(2)*cos(args)
  sines = sqrt(2)*sin(args)
  ones = as.matrix(rep(1,L))
  X = cbind(ones, t(interleave(t(cosines), t(sines))))
  return(X)
}

basis_functions_x = function(K, xs){
  # K is the (even) number of basis functions to compute 
  args = seq(1,K/2,1)*2*pi
  args = as.matrix(xs)%*%t(as.matrix(args))
  # L is the size of the range discretisation
  cosines = sqrt(2)*cos(args)
  sines = sqrt(2)*sin(args)
  ones = as.matrix(rep(1,length(xs)))
  X = cbind(ones, t(interleave(t(cosines), t(sines))))
  return(X)
}

phi_bars = function(K,xs){
  temp = basis_functions_x(K,xs)
  bars = apply(temp,2,mean)
  return(bars)
}

compute_fn_from_coefs = function(coefs, L=100){
  # Computes the function on [0,1] on a grid of size L
  K = length(coefs) - 1
  basis_fn_vals = basis_functions(K,L)
  fn_vals = basis_fn_vals%*%coefs
  return(fn_vals)
}

compute_fn_from_coefs_x = function(coefs, xs){
  K = length(coefs) - 1
  basis_fn_vals = basis_functions_x(K,xs)
  fn_vals = basis_fn_vals%*%coefs
  return(fn_vals)
}

sample_from_density = function(f_0, num_samples){
  L = length(f_0)
  x_range = seq(0, 1, length.out=L)
  sample = sample(x_range,size=num_samples,prob=f_0,replace=TRUE)
  return(sample)
}

n=10000 #number of training points
alpha = 1
K=401 #truncation point
L=1000

#true smoothness of log(f)
beta = 1/2
#coeffs of log(f)
g_ks = (1:K)^(-1/2 - beta)
#log(f)
g = compute_fn_from_coefs(g_ks, L=L)
#f
f = exp(g)/sum(exp(g))

#smoothness of a
mu = 1
a_ks = (1:K)^(-1/2 - mu)
a = compute_fn_from_coefs(a_ks, L=L)

gamma = 1

posterior_sample = function(beta,mu,gamma,n=10000,L=1000,alpha=1,num_samples=1000,rescale=FALSE){
  
  #compute true density and sample from it
  g_ks = (1:K)^(-1/2 - beta)
  g = compute_fn_from_coefs(g_ks, L=L)
  # f = exp(g)/sum(exp(g))
  f = exp(g)/sum(exp(g))*L
  xn = sample_from_density(f_0 = f, num_samples = n)
  
  #compute functional
  a_ks = (1:K)^(-1/2 - mu)
  a = compute_fn_from_coefs(a_ks, L=L)
  
  inverse_lengthscale = n^(1/(1+2*gamma))
  # prox = 1/inverse_lengthscale
  prox = exp(-inverse_lengthscale^2/100)
  
  fit = sbde(xn, nsamp=num_samples, fbase="unif", par="pmean", blocking="gp",
             temp=alpha, prox.range = prox)
  
  pred = predict(fit, burn=0.5, nmc=num_samples, yRange=c(0,1),yLength=L)
  
  psi_samples = apply(pred$fsamp, 2, function(t) mean(t*a))
  
  #rescale
  if(rescale){
    #recenter
    psi_0 = sum(f*a)
    psi_hat = mean(compute_fn_from_coefs_x(a_ks,xn))
    psi_samples = psi_samples - mean(psi_samples)
    #compute quantities for rescaling
    norm_a_sq = mean(a*a*f)
    norm_eff_inf_sq = norm_a_sq - (mean(a*f))^2
    sig = sqrt(norm_eff_inf_sq)/2
    #rescale
    psi_samples = sqrt(n*alpha)*psi_samples/sig
  }
  
  return(psi_samples)
}

alpha=1/4
#compute samples for different beta
betas = c(0.1, 1)
mu = 1
gamma = 1
rescaled_samples = sapply(betas, function(b) posterior_sample(beta=b,mu,gamma, alpha = alpha,rescale=TRUE,num_samples=10000))

beta_plot_dat = data.frame(beta_0_1 = rescaled_samples[,1], beta_1 = rescaled_samples[,2]) %>%  
  gather(key='beta',value='psi_f')


density_x = seq(-6,6,length.out=100)
ref_density = dnorm(density_x)
# ref_dat = data.frame(density_x=density_x,ref_density = ref_density,beta=rep('beta',100))
ref_dat = data.frame(density_x=density_x,ref_density = ref_density)

beta_plot = ggplot(ref_dat, aes(x=density_x, y = ref_density)) + geom_line() +
  geom_histogram(data=beta_plot_dat, aes(x=psi_f, y=..density.., fill=beta), position='identity', alpha=0.5, bins=50) +
  xlab(TeX(r'($\frac{\sqrt{n\alpha_n}(\psi(f) - \hat{\psi})}{||\tilde{\psi} ||_L})')) +
  scale_fill_discrete(name = TeX(r'($\beta)'), labels = c("0.1", "1")) +
  ylab("Density") + 
  # ggtitle("Rescaled samples from the Fractional Posterior (alpha = 1/4)") +
  xlim(-6,6) +
  ylim(0,0.5) +
  theme(text = element_text(size = 13))

beta_plot
#compute samples for different gamma
beta = 1
mu = 1
gammas = c(10,1)
rescaled_samples_gamma = sapply(gammas, function(g) posterior_sample(beta,mu,gamma=g, alpha=alpha, rescale=TRUE,num_samples=10000))

#note unconventional naming of dataframe columns to get colors to match 'right' and 'wrong'
gamma_plot_dat = data.frame(gamma_1 = rescaled_samples_gamma[,1], gamma_2 = rescaled_samples_gamma[,2]) %>% 
  gather(key='gamma',value='psi_f')


density_x = seq(-6,6,length.out=100)
ref_density = dnorm(density_x)
# ref_dat = data.frame(density_x=density_x,ref_density = ref_density,beta=rep('beta',100))
ref_dat = data.frame(density_x=density_x,ref_density = ref_density)

gamma_plot = ggplot(ref_dat, aes(x=density_x, y = ref_density)) + geom_line() +
  geom_histogram(data=gamma_plot_dat, aes(x=psi_f, y=..density.., fill=gamma), position='identity', alpha=0.5, bins=50) +
  xlab(TeX(r'($\frac{\sqrt{n\alpha_n}(\psi(f) - \hat{\psi})}{||\tilde{\psi} ||_L})')) +
  scale_fill_discrete(name = TeX(r'($\gamma)'), labels = c("10", "1")) +
  ylab(" ") + 
  xlim(-6,6) +
  ylim(0,0.5) +
  # ggtitle('') +
  theme(text = element_text(size = 13))

combined = grid.arrange(beta_plot, gamma_plot, nrow = 1)
ggsave(plot = combined, '/Users/lmt15/Documents/phd/Thesis/Writing/Figures/gp_bvm_condition_illustration2.pdf', width = 10, height = 5, units='in')
#NOTE: Save as 8x14 pdf

### To compare credible sets with corrected credible sets
beta = 1
mu = 1
gamma = 1

L = 1000
K = 401

#compute truth
g_ks = (1:K)^(-1/2 - beta)
#log(f)
g = compute_fn_from_coefs(g_ks, L=L)
#f
f = exp(g)/sum(exp(g))*L

#compute functional
a_ks = (1:K)^(-1/2 - mu)
a = compute_fn_from_coefs(a_ks, L=L)

psi_0 = mean(a*f)

alpha = n^{-1/4}
xn = seq(0,1,length.out=1000)
sample_dat = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=alpha, 
                        n=10000, L = 1000, num_samples = 10000, 
                        rescale=FALSE)
sample_dat = sample_dat - mean(sample_dat) + psi_0

alpha_lower_bound = quantile(sample_dat, 0.025)
alpha_upper_bound = quantile(sample_dat, 0.975)

corrected_lower_bound = psi_0 + sqrt(alpha)*(alpha_lower_bound - psi_0)
corrected_upper_bound = psi_0 + sqrt(alpha)*(alpha_upper_bound - psi_0)

sample_dat_2 = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=1, 
                              n=10000, L = 1000, num_samples = 10000, 
                              rescale=FALSE)
sample_dat_2 = sample_dat_2 - mean(sample_dat_2) + psi_0

lower_bound_2 = quantile(sample_dat_2, 0.025)
upper_bound_2 = quantile(sample_dat_2, 0.975)

shading = 0.5

plot_1 = ggplot(data.frame(psi_f = sample_dat), aes(x=psi_f,y=..density..)) + 
  geom_histogram(bins=50, alpha=0.75) + 
  geom_vline(aes(xintercept=alpha_lower_bound), color='red', alpha=shading) + 
  geom_vline(aes(xintercept=alpha_upper_bound), color='red', alpha=shading) +
  geom_vline(aes(xintercept=corrected_lower_bound), color='blue', alpha=shading) +
  geom_vline(aes(xintercept=corrected_upper_bound), color='blue', alpha=shading) +
  scale_color_discrete(name='Set Type') +
  xlim(psi_0 - 0.03, psi_0 + 0.03) +
  ylim(0,150) +
  ggtitle(TeX(r'(Sample from the $\alpha_n$-posterior ($\alpha_n = n^{-1/4}$))')) +
  theme(axis.text = element_text(size = 13)) +
  xlab(TeX(r'($\psi(f)$)')) +
  ylab('Density')

plot_2 = ggplot(data.frame(psi_f = sample_dat_2), aes(x=psi_f,y=..density..)) + 
  geom_histogram(bins=85, alpha=0.75) + 
  geom_vline(aes(xintercept=lower_bound_2), color='#198038', alpha=0.75) + 
  geom_vline(aes(xintercept=upper_bound_2), color='#198038', alpha=0.75) +
  scale_color_discrete(name='Set Type') +
  xlim(psi_0 - 0.03, psi_0 + 0.03) +
  ylim(0,150) +
  ggtitle('Sample from the full posterior') +
  theme(axis.text = element_text(size = 13)) +
  xlab(TeX(r'($\psi(f)$)')) +
  ylab('Density')

grid.arrange(plot_1, plot_2, nrow=2)
#Note: save as 8x10


### To compare corrected credible sets with different alpha decays
beta = 1
mu = 1
gamma = 1

L = 1000
K = 401

#compute truth
g_ks = (1:K)^(-1/2 - beta)
#log(f)
g = compute_fn_from_coefs(g_ks, L=L)
f = exp(g)/sum(exp(g))*L

#compute functional
a_ks = (1:K)^(-1/2 - mu)
a = compute_fn_from_coefs(a_ks, L=L)

psi_0 = mean(a*f)

xn = seq(0,1,length.out=1000)

sample_dat_full = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=1, 
                                n=n, L = 1000, num_samples = 10000, 
                                rescale=FALSE)
centering = mean(sample_dat_full)

alpha_slow = n^{-1/4}*sqrt(log(n))
sample_dat_alpha_slow = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=alpha_slow,
                              n=n, L = 1000, num_samples = 10000,
                              rescale=FALSE)
alpha_fast = n^{-1/4}/sqrt(log(n))
sample_dat_alpha_fast = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=alpha_fast,
                                         n=n, L = 1000, num_samples = 10000,
                                         rescale=FALSE)

# alpha_equal = n^{-1/4}
# sample_dat_alpha_equal = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=alpha_equal, 
#                                          n=n, L = 1000, num_samples = 10000, 
#                                          rescale=FALSE)
{
{
sample_dat_full = sample_dat_full + (psi_0 - centering)
sample_dat_alpha_slow = sample_dat_alpha_slow + (psi_0 - centering)
sample_dat_alpha_fast = sample_dat_alpha_fast + (psi_0 - centering)
sample_dat_alpha_equal = sample_dat_alpha_equal + (psi_0 - centering)

lower_bound_full = quantile(sample_dat_full, 0.025)
upper_bound_full = quantile(sample_dat_full, 0.975)

lower_bound_alpha_slow = quantile(sample_dat_alpha_slow, 0.025)
upper_bound_alpha_slow = quantile(sample_dat_alpha_slow, 0.975)

lower_bound_alpha_fast = quantile(sample_dat_alpha_fast, 0.025)
upper_bound_alpha_fast = quantile(sample_dat_alpha_fast, 0.975)

lower_bound_alpha_equal = quantile(sample_dat_alpha_equal, 0.025)
upper_bound_alpha_equal = quantile(sample_dat_alpha_equal, 0.975)

corrected_lower_bound_slow =  mean(sample_dat_alpha_slow) + sqrt(alpha_slow)*(lower_bound_alpha_slow - mean(sample_dat_alpha_slow))
corrected_upper_bound_slow =  mean(sample_dat_alpha_slow) + sqrt(alpha_slow)*(upper_bound_alpha_slow - mean(sample_dat_alpha_slow))

corrected_lower_bound_fast =  mean(sample_dat_alpha_fast) + sqrt(alpha_fast)*(lower_bound_alpha_fast - mean(sample_dat_alpha_fast))
corrected_upper_bound_fast =  mean(sample_dat_alpha_fast) + sqrt(alpha_fast)*(upper_bound_alpha_fast - mean(sample_dat_alpha_fast))

corrected_lower_bound_equal =  mean(sample_dat_alpha_equal) + sqrt(alpha_equal)*(lower_bound_alpha_equal - mean(sample_dat_alpha_equal))
corrected_upper_bound_equal =  mean(sample_dat_alpha_equal) + sqrt(alpha_equal)*(upper_bound_alpha_equal - mean(sample_dat_alpha_equal))

cred_regions = data.frame(alpha = c('Slow', 'Fast', 'Full'), 
                          lower_bounds = c(corrected_lower_bound_slow, corrected_lower_bound_fast, lower_bound_full),
                          upper_bounds = c(corrected_upper_bound_slow, corrected_upper_bound_fast, upper_bound_full))

# sample_dat_2 = sample_dat_2 - mean(sample_dat_2) + psi_0

sample_dat_alpha = data.frame(c1 = sample_dat_alpha_slow, 
                              c2 = sample_dat_alpha_fast,
                              c3 = sample_dat_full)
colnames(sample_dat_alpha) = c('Slow', 'Fast', 'Full')
sample_dat_alpha = sample_dat_alpha %>% gather(key = 'alpha', value = 'sample')

uncorrected_cred_regions  = data.frame(alpha = c('Slow', 'Fast', 'Full'), 
                                      lower_bounds = c(lower_bound_alpha_slow, lower_bound_alpha_fast, lower_bound_full),
                                      upper_bounds = c(upper_bound_alpha_slow, upper_bound_alpha_fast, upper_bound_full))

cred_regions['type'] = rep('Corrected', 3)
uncorrected_cred_regions['type'] = rep('Normal', 3)


shading = 0.5

sample_dat_alpha$alpha = factor(sample_dat_alpha$alpha, levels=c('Fast', 'Slow', 'Full'))
cred_regions$alpha = factor(cred_regions$alpha, levels=c('Fast', 'Slow', 'Full'))
uncorrected_cred_regions$alpha = factor(cred_regions$alpha, levels=c('Fast', 'Slow', 'Full'))
}
#Plot for reviewers
ggplot(sample_dat_alpha, aes(x=sample,y=..density.., fill = alpha)) + 
  geom_histogram(bins=50, alpha=0.75) + 
  geom_vline(data = uncorrected_cred_regions, aes(xintercept=lower_bounds, color = alpha), linetype = 'dashed') +
  geom_vline(data = uncorrected_cred_regions, aes(xintercept=upper_bounds, color = alpha), linetype = 'dashed') +
  geom_vline(data = cred_regions, aes(xintercept=lower_bounds, color = alpha)) +
  geom_vline(data = cred_regions, aes(xintercept=upper_bounds, color = alpha)) +
  scale_color_discrete(name=TeX(r'($\alpha_n$)'), labels=c(TeX(r'($n^{-1/4}/sqrt(ln(n))$)'), TeX(r'($n^{-1/4}\cdot sqrt(ln(n))$)'), '1')) +
  scale_fill_discrete(name=TeX(r'($\alpha_n$)'), labels=c(TeX(r'($n^{-1/4}/sqrt(ln(n))$)'), TeX(r'($n^{-1/4}\cdot sqrt(ln(n))$)'), '1')) +
  xlim(1.2, 1.27) +
  facet_wrap(~alpha, nrow = 3) +
  ggtitle(TeX(r'(Sample from the $\alpha_n$-posterior)')) +
  theme(axis.text = element_text(size = 13)) +
  xlab(TeX(r'($\psi(f)$)')) +
  ylab('Density') 

#Plot for Figure 2
alpha_plot = ggplot(data.frame(psi = sample_dat_alpha_equal), aes(x = psi, y = ..density..)) +
  geom_histogram(bins = 75, alpha = 0.75) +
  geom_vline(aes(xintercept = lower_bound_alpha_equal), color = 'red', alpha = 0.75) +
  geom_vline(aes(xintercept = upper_bound_alpha_equal), color = 'red', alpha = 0.75) + 
  geom_vline(aes(xintercept = corrected_lower_bound_equal), color = 'blue', alpha = 0.75) +
  geom_vline(aes(xintercept = corrected_upper_bound_equal), color = 'blue', alpha = 0.75) +
  ylab('Density') +
  xlab(TeX(r'($\psi(f)$)')) +
  xlim(1.22, 1.28) +
  ylim(0, 150) +
  ggtitle(TeX(r'(Sample from the $\alpha_n$-posterior ($\alpha_n = n^{-1/4}$))'))

full_plot = ggplot(data.frame(psi = sample_dat_full), aes(x = psi, y = ..density..)) +
  geom_histogram(bins = 75, alpha = 0.75) +
  geom_vline(aes(xintercept = lower_bound_full), color='#198038', alpha = 0.75) +
  geom_vline(aes(xintercept = upper_bound_full), color='#198038', alpha = 0.75) + 
  ylab('Density') +
  xlab(TeX(r'($\psi(f)$)')) +
  xlim(1.22, 1.28) + 
  ylim(0, 150) +
  ggtitle(TeX(r'(Sample from the full posterior)')) 


grid.arrange(alpha_plot, full_plot,nrow=2)
}



mu = mean(sample_dat_full)
sig = sd(sample_dat_full)

mu_alpha_slow = mean(sample_dat_alpha_slow)
mu_alpha_fast = mean(sample_dat_alpha_fast)

mu_alpha_slow
pnorm((mu + 1.96*sig - mu_alpha_slow)/sig) - pnorm((mu - 1.96*sig - mu_alpha_slow)/sig)


sim_cred_set = function(n, beta, mu, gamma){
  L = 1000
  K = 401
  #compute truth
  g_ks = (1:K)^(-1/2 - beta)
  #log(f)
  g = compute_fn_from_coefs(g_ks, L=L)
  #f
  f = exp(g)/sum(exp(g))*L
  
  #compute functional
  a_ks = (1:K)^(-1/2 - mu)
  a = compute_fn_from_coefs(a_ks, L=L)
  
  psi_0 = mean(a*f)
  
  xn = seq(0,1,length.out=1000)
  
  sample_dat_full = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=1, 
                                     n=10000, L = 1000, num_samples = 10000, 
                                     rescale=FALSE)
  centering = mean(sample_dat_full)
  
  alpha_slow = n^{-1/4}*log(n)
  sample_dat_alpha_slow = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=alpha_slow, 
                                           n=10000, L = 1000, num_samples = 10000, 
                                           rescale=FALSE)
  alpha_fast = n^{-1/4}/log(n)
  sample_dat_alpha_fast = posterior_sample(beta=beta, mu=mu, gamma=gamma, alpha=alpha_fast, 
                                           n=10000, L = 1000, num_samples = 10000, 
                                           rescale=FALSE)
  
  noise = rnorm(1, mean = 0, sd = sd(sample_dat_full))
  sample_dat_full = sample_dat_full + (psi_0 - centering + noise)
  sample_dat_alpha_slow = sample_dat_alpha_slow + (psi_0 - centering + noise)
  sample_dat_alpha_fast = sample_dat_alpha_fast + (psi_0 - centering + noise)
  
  lower_bound_full = quantile(sample_dat_full, 0.025)
  upper_bound_full = quantile(sample_dat_full, 0.975)
  
  lower_bound_alpha_slow = quantile(sample_dat_alpha_slow, 0.025)
  upper_bound_alpha_slow = quantile(sample_dat_alpha_slow, 0.975)
  
  lower_bound_alpha_fast = quantile(sample_dat_alpha_fast, 0.025)
  upper_bound_alpha_fast = quantile(sample_dat_alpha_fast, 0.975)
  
  corrected_lower_bound_slow =  mean(sample_dat_alpha_slow) + sqrt(alpha_slow)*(lower_bound_alpha_slow - mean(sample_dat_alpha_slow))
  corrected_upper_bound_slow =  mean(sample_dat_alpha_slow) + sqrt(alpha_slow)*(upper_bound_alpha_slow - mean(sample_dat_alpha_slow))
  
  corrected_lower_bound_fast =  mean(sample_dat_alpha_fast) + sqrt(alpha_fast)*(lower_bound_alpha_fast - mean(sample_dat_alpha_fast))
  corrected_upper_bound_fast =  mean(sample_dat_alpha_fast) + sqrt(alpha_fast)*(upper_bound_alpha_fast - mean(sample_dat_alpha_fast))
  
  full_posterior_stats = cred_stats(lower_bound_full, upper_bound_full, psi_0)
  alpha_posterior_stats_slow = cred_stats(lower_bound_alpha_slow, upper_bound_alpha_slow, psi_0)
  alpha_posterior_stats_fast = cred_stats(lower_bound_alpha_fast, upper_bound_alpha_fast, psi_0)
  ss_posterior_stats_slow = cred_stats(corrected_lower_bound_slow, corrected_upper_bound_slow, psi_0)
  ss_posterior_stats_fast = cred_stats(corrected_lower_bound_fast, corrected_upper_bound_fast, psi_0)
  return(rbind(full_posterior_stats, 
               alpha_posterior_stats_slow, 
               alpha_posterior_stats_fast,
               ss_posterior_stats_slow,
               ss_posterior_stats_fast))
}

cred_stats = function(lower, upper, actual){
  return(c((lower <= actual & actual <= upper), upper - lower, (lower + upper)/2-actual) )
}

t1 = Sys.time()
temp = replicate(50, sim_cred_set(1000, 1, 1, 1))
t2 = Sys.time()
cat('Time Elapsed: ', difftime(t2, t1, units = 'mins'), '\n')
full_post = apply(temp[1,,], 1, mean)
alpha_slow = apply(temp[2 ,,], 1, mean)
alpha_fast = apply(temp[3 ,,], 1, mean)
ss_slow = apply(temp[4 ,,], 1, mean)
ss_fast = apply(temp[5 ,,], 1, mean)


#slow
bias = mean(sample_dat_alpha_slow) - psi_0
len = upper_bound_alpha_slow - lower_bound_alpha_slow
sig = len/1.96
cat('Original, slow alpha')
cat('Bias: ', bias, '\n')
cat('Len: ', len, '\n')
cat('Prob: ', pnorm((1.96*sig - bias)/sig)  - pnorm((-1.96*sig - bias)/sig), '\n')

len = corrected_upper_bound_slow - corrected_lower_bound_slow
sig = len/1.96
cat('Corrected, slow alpha')
cat('Bias: ', bias, '\n')
cat('Len: ', len, '\n')
cat('Prob: ', pnorm((1.96*sig - bias)/sig)  - pnorm((-1.96*sig - bias)/sig), '\n')

#fast
bias = mean(sample_dat_alpha_fast) - psi_0
len = upper_bound_alpha_fast - lower_bound_alpha_fast
sig = len/1.96
cat('Original, fast alpha')
cat('Bias: ', bias, '\n')
cat('Len: ', len, '\n')
cat('Prob: ', pnorm((1.96*sig - bias)/sig)  - pnorm((-1.96*sig - bias)/sig), '\n')

len = corrected_upper_bound_fast - corrected_lower_bound_fast
sig = len/1.96
cat('Corrected, fast alpha')
cat('Bias: ', bias, '\n')
cat('Len: ', len, '\n')
cat('Prob: ', pnorm((1.96*sig - bias)/sig)  - pnorm((-1.96*sig - bias)/sig), '\n')

#full
bias = mean(sample_dat_full) - psi_0
len = upper_bound_full - lower_bound_full
sig = len/1.96
cat('Original, fast alpha')
cat('Bias: ', bias, '\n')
cat('Len: ', len, '\n')
cat('Prob: ', pnorm((1.96*sig - bias)/sig)  - pnorm((-1.96*sig - bias)/sig), '\n')


n = 1e4
alpha_slow = n^{-1/3.5}
alpha_fast = n^{-1/4.5}
cat('Slow: ', alpha_slow, '\n')
cat('Fast: ', alpha_fast, '\n')

