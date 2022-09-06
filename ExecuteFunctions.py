from EconPy import *

''' Iteration technique '''

var = 0.00001
# Fill in condition 
condition = 500*PA(var*100, 25) # HERE 
# Fill in target value
target = 8000 # HERE
while abs(condition - target) / condition > 0.0001:
    var += 0.00001
    # Copy and paste condition
    condition = 500*PA(var*100, 25) # HERE
print(var*100)


