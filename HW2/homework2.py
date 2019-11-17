import methods as mine

problem = input("Test problem C1-P")

#C1-P9
#Consider the function
#f(x) = e^(sin^3(x)) + x^6 - 2x^4 - x^3 - 1
#on the interval [-2,2]. Plot the function on the interval, and
#find all three roots to six correct decimal places. Determine
#which roots converge quadratically, and find the multiplicity
#of the roots that converge linearly.
if problem == '9':
    from numpy import exp, sin, linspace, cos
    from matplotlib.pyplot import plot, show
    def f(x):
        return exp(sin(x)**3) + x**6 - 2*x**4 - x**3 - 1

    def fp(x):
        return -3*x**2 - 8*x**3 + 6*x**5 + 3*exp(sin(x)**3)*cos(x)*sin(x)**2

    x = linspace(-2,2,60,endpoint = True)
    y = f(x)
    plot(x,y)
    show()

    tol = 0.5e-6
    error = 1
    x_guess = [-1.2,0,1.5]
    for i in range(3):
        while error > tol:
            new_guess = x_guess[i] - f(x_guess[i]) / fp(x_guess[i])
            error = abs(new_guess - x_guess[i])
            x_guess[i] = new_guess
        x_guess[i] = round(x_guess[i],6)
        print('The root is',x_guess[i])
        hasMult = fp(x_guess[i])
        if hasMult == 0:
            print('Root',x_guess[i],'converges linearly and has multiplicity 4') #Calculated in Mathematica
        else:
            print('Root',x_guess[i],'converges quadratically')

#C1-P11
#Compare the results obtained using the Bisection Method,
#Newton's Method, and the Secant Method to solve the
#equation ln(x) + x^2 = 3.
    #Solve the problem using each of the three methods.
    #Plot the errors on the same graph using a log-scaled
    #y-axis. Explanation on paperwork.
if problem == '11':

    from numpy import log, abs
    from matplotlib.pyplot import plot, show, legend, yscale

    def f(x):
        return log(x) + x**2 - 3
    def fp(x):
        return 1 / x + 2 * x
    
    solver = mine.RootFinding(f,0.5e-6,2,x1 = 1,fp=fp,plot = True)
    BisData, BisRoot = solver.BisectionMethod()
    NewtData, NewtRoot = solver.NewtonsMethod()
    SecData, SecRoot = solver.SecantMethod()

    print('Bisection Method solves with a solution of',round(BisRoot,6))
    print('Newtons Method solves with a solution of',round(NewtRoot,6))
    print('Secant Method solves with a solution of',round(SecRoot,6))

    plot(abs(BisData - BisRoot), label = "Bisection")
    plot(abs(NewtData - NewtRoot), label = "Newton")
    plot(abs(SecData - SecRoot), label = "Secant")
    yscale("log")
    legend()
    show()

else:
    print("This problem is not here! I have problems C1-P9 and C1.P11.")
