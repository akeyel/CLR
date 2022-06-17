#### Simple Power analysis for logistic regression

#### Overview ####

# Delta Emergence analysis had non-significant results.
# However, the same magnitude Odds Ratios were significant for Omicron.
# Goal is to demonstrate that non-significance for the Delta Emergence Analysis
# was a function of low statistical power, and that a real effect is the
# most likely outcome.

# Note that simple logistic regression for a binary variable was used here
# to simplify the power calculations. A failed attempt was made to use the
# samplesizeCMH package for conditional logistic regression.

#### Design ####
# Calculate power for an OR of 2 for the Delta Emergence Sample size
# Calculate power for changing sample sizes 50 - 350
# Calculate power for changing OR's: 2,3,4 for changing sample sizes

# Use WebPower for logistic regression
#install.packages('WebPower')
library(WebPower)

#### Define a function for converting from outcome probabilities to Odds Ratios
calculate.odds.ratio = function(p1, p2){
  odds1 = p1 / (1-p1)
  odds2 = p2 / (1-p2)
  OR = odds1/odds2
  return(OR)
}

#### Identify probability pairs that produce an OR of ~2
# These were adjusted by trial and error
calculate.odds.ratio(0.1, 0.0526) #2.001
calculate.odds.ratio(0.2, 0.111) # 2.002
calculate.odds.ratio(0.4, 025)   # exactly 2
# etc.

# Look at probability in increments of 0.1
p0.vec = seq(0.1,0.9, 0.1)
p1.vec = c(0.0526, 0.111, 0.176, 0.25, 0.333,0.428, 0.528, 0.666, 0.818)

# Calculate power for OR 2 across a range of probability pairs
for (i in 1:length(p1.vec)){
  p0 = p0.vec[i]
  p1 = p1.vec[i]
  power = wp.logistic(p0 = p0, p1 = p1, n = 110, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
  print(power$power)
}
# 0.1 corresponds to 0.15 minimum
# 0.7 corresponds to 0.45 maximum

# Examine effect of power for 0.7, 0.528 pair, which had maximum power
p0 = p0.vec[7]
p1 = p1.vec[7]
for (n in c(50,100,150,200,250,300,350)){
  power = wp.logistic(p0 = p0, p1 = p1, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
  print(power$power)
}
# 79% power for sample size 250, so 80% power is >250 (but not much greater)

# Identify where 80% power is reached
p0 = p0.vec[7]
p1 = p1.vec[7]
for (n in seq(250,260)){
  print(n)
  power = wp.logistic(p0 = p0, p1 = p1, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
  print(power$power)
}

# Examine power different ORs at p = 0.70
# OR = 4
p0 = 0.7
p1 = 0.368
for (n in c(50,100,150,200,250,300,350)){
  power = wp.logistic(p0 = p0, p1 = p1, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
  print(power$power)
}
power = wp.logistic(p0 = p0, p1 = p1, n = 110, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
print(power$power)


# OR = 3
p0 = 0.7
p1 = 0.437
for (n in c(50,100,150,200,250,300,350)){
  power = wp.logistic(p0 = p0, p1 = p1, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
  print(power$power)
}
power = wp.logistic(p0 = p0, p1 = p1, n = 110, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
print(power$power)

# With 24 samples, how big would OR have to be?
# OR = 4
p0 = 0.7
p1 = 0.368
n = 24 # 12 pairs
power = wp.logistic(p0 = p0, p1 = p1, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
# 0.20 power to detect 

calculate.odds.ratio(0.7,0.1) #21
power = wp.logistic(p0 = 0.7, p1 = 0.1, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
# 0.75 for an OR of 21

calculate.odds.ratio(0.7,0.05) #44.3
power = wp.logistic(p0 = 0.7, p1 = 0.05, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
# 0.73 for an OR of 44.3 # Lower than for 21, which is odd, but expected due to the shift away from the center.

calculate.odds.ratio(0.7,0.01)
power = wp.logistic(p0 = 0.7, p1 = 0.01, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
# 0.45. Statistical power is even lower. 

power = wp.logistic(p0 = 0.77, p1 = 0.13, n = n, alpha = 0.05, alternative = 'two.sided', family = 'Bernoulli')
power$power
calculate.odds.ratio(0.77,0.13)
