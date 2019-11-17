from numpy import sin

def f122a(x):
    return 1/(1 + x**4)

def f122b(x):
    return (sin(x) - 5)/6

error = 1
guess1 = -0.123
while error > 1e-8:
  temp_guess = f122a(guess1)
  error = abs(temp_guess - guess1)
  guess1 = temp_guess

print("The solution to x^5 + x = 1 is ",guess1)

error = 1
guess1 = -0.123
while error > 1e-8:
  temp_guess = f122b(guess1)
  error = abs(temp_guess - guess1)
  guess1 = temp_guess

print("The solution to sin(x) = 6x + 5 is ",guess1)


