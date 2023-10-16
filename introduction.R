##intermediate variables 

Y = C + G                         # Output
Tax = Theta * Y                   # Tax 
Yd = Y - Tax                      # Disposable income
C_target= a1*Yd+a2*S              # Consumption targeted

##time derivatives
S = Yd - C                        # Savings
C = beta_C*(C_target-C)
##initial values
S = 0
C = 20
##parameters

Theta = 0.2                         # tax rate 
a1 = 0.6                            # Marginal propensity to consume out of income
a2 = 0.4                            # Marginal propensity to consume out of wealth
G = 20                              # Exogenously defined government spending
beta_C = 2


#_______________________________________________________________________________________________________________________________#
##time
begin = 0                      # Note that time can be used as a variable "t"
end = 84
by = 0.05