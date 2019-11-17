from numpy import sqrt,exp,log,sin,arcsinh,pi
from matplotlib.pyplot import plot, show,yscale,xscale
#For integrals in Computer Problem 1, calculate the approximation error of the composite
#Simpson's Rule for h = (b - a)/2, h/2, h/4, .../c/2^8, and plot. Make a log-log plot.
#What is the slope of the plot and does it agree with theory?

#Define Simpson's Method
def Simpson(f,a,b):
    h = (b-a)/(2)
    return h * (f(b) + 4*f(a+h) + f(a)) / 3

#Define Composite of any method
def Comp(f,a,b,nPanels):
    summation = 0
    for i in range(nPanels):
        summation += Simpson(f,a+i*(b-a)/nPanels,a+(i+1)*(b-a)/nPanels)
    return summation

#Define our functions
def fa(x):
    return x/sqrt(x**2 + 9)

def fb(x):
    return x**3 / (x**2 + 1)

def fc(x):
    return x * exp(x)

def fd(x):
    return x**2 * log(x)

def fe(x):
    return x**2 * sin(x)

def ff(x):
    return x**3 / sqrt(x**4 - 1)

def fg(x):
    return 1/sqrt(x**2 + 4)

def fh(x):
    return x/sqrt(x**4 + 1)

#Put them, solutions, and bounds into lists
functions = [fa,fb,fc,fd,fe,ff,fg,fh]
solutions = [2,(1-log(2))/2,1,-26/9+9*log(3),pi**2-4,-sqrt(5)*(sqrt(3)-4)/2,arcsinh(sqrt(3)),arcsinh(1)/2]
a = [0,0,0,1,0,2,0,0]
b = [4,1,1,3,pi,3,2*sqrt(3),1]
for j in range(len(functions)): #iterate over all functions
    error = []
    for i in range(8): #iterate over specified number of panels
        error.append(abs(solutions[j]-Comp(functions[j],a[j],b[j],2**i)))
    plot(error)
    yscale('log')
    show()
