import sys
sys.path.append('../')
import methods as mine

problem = input("Test problem C5-P")

#C5-P6 a.k.a. Computer Problem 5.2.1ac
#Use the composite Trapezoid Rule wiht m = 16 and 32 panels to approximate the definite
#integral. Compare with the correct integral and report the two errors.
if problem == '6':
    print("\nPart a")
    def fa(x):
        from numpy import sqrt
        return x / sqrt(x**2 + 9)
    
    a = 0
    b = 4
    m = 16
    integratorA = mine.Integration(a,b,fa,m)
    approx16 = integratorA.Trap()
    approx32 = integratorA.Trap(32)
    actual = 2
    print("The approximation for 16 panels is ",approx16," with an error of ",abs(approx16 - actual))
    print("The approximation for 32 panels is ",approx32," with an error of ",abs(approx32 - actual))

    print("\nPart b")

    def fc(x):
        from numpy import exp
        return x * exp(x)

    a = 0
    b = 1
    m = 16
    integratorC = mine.Integration(a,b,fc,m)
    approx16 = integratorC.Trap()
    approx32 = integratorC.Trap(32)
    actual = 1
    print("The approximation for 16 panels is ",approx16," with an error of ",abs(approx16 - actual))
    print("The approximation for 32 panels is ",approx32," with an error of ",abs(approx32 - actual))

#C5-P7 a.k.a. Computer Problem 5.2.2de
#Apply the composite Simpson's Rule to the integrals in Computer Problem 1. Use m = 16
#and 32, and report errors.
if problem == '7':
    print("\nPart d")
    def fd(x):
        from numpy import log
        return x**2 * log(x)

    from numpy import log
    a = 1
    b = 3
    m = 16
    integratorD = mine.Integration(a,b,fd,m)
    approx16 = integratorD.Simpson()
    approx32 = integratorD.Simpson(32)
    actual = 9 * log(3) - 26 / 9
    print("The approximation for 16 panels is ",approx16," with an error of ",abs(approx16 - actual))
    print("The approximation for 32 panels is ",approx32," with an error of ",abs(approx32 - actual))
   
    print("\nPart e")
    def fe(x):
        from numpy import sin
        return x**2 * sin(x)

    from numpy import pi
    a = 0
    b = pi
    m = 16
    integratorE = mine.Integration(a,b,fe,m)
    approx16 = integratorE.Simpson()
    approx32 = integratorE.Simpson(32)
    actual = pi**2 - 4
    print("The approximation for 16 panels is ",approx16," with an error of ",abs(approx16 - actual))
    print("The approximation for 32 panels is ",approx32," with an error of ",abs(approx32 - actual))

#C5-P8 a.k.a. Computer Problem 5.2.9bf
#For the integrals in Computer Problem 1, calculate the approximation error of the
#composite Trapezoid Rule for h = b - a,h/2,h/4,...,h/2^8, and plot. Make a log-log
#plot. What is the slope of the plot, and does it agree with theory?
if problem == '8':
    print("\nPart b")
    def fb(x):
        return x**3 / (x**2 + 1)

    from numpy import log
    a = 0
    b = 1
    integratorB = mine.Integration(a,b,fb)
    actual = (1 - log(2)) / 2
    error = []
    for i in range(25):
        error.append(abs(actual - integratorB.Trap(2**i)))
    
    print("This agrees perfectly with theory, a constant slope")
    from matplotlib.pyplot import plot, show, yscale
    plot(error)
    yscale('log')
    show()

    print("\nPart f")
    def ff(x):
        from numpy import sqrt
        return x**3 / sqrt(x**4 - 1)

    from numpy import sqrt
    a = 2
    b = 3
    integratorF = mine.Integration(a,b,ff)
    actual = -sqrt(5) * (sqrt(3) - 4) / 2
    error = []
    for i in range(25):
        error.append(abs(actual - integratorF.Trap(2*i)))

    print("This one is very strange. It definitely follows a pattern, but the error jumps around a lot")
    plot(error)
    yscale('log')
    show()
