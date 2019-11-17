import sys
sys.path.append('../')
import methods as mine

problem = input("Test problem C3-P")

#C3-P7
#Rebuild Program 3.3 to implement the Chebyshev interpolation polynomial with four
#nodes on the interval [0,pi/2]. Then plot the polynomial and the sine function on
#the interval [-2,2]. To demonstrate that you function works, interpolate sin(x) on
#the interval [-pi,pi] using the nodes -pi,-pi/2,0,pi/2,pi. Plot your interpolating
#polynomial Use the numpy functions polyfit and polyval to plot an interpolating
#polynomial from Python on the same graph. Label the graph so it is clear which curve
#is yours and which on came from Python. Also, plot sin(x) on the same graph.
if problem == '7':
    from numpy import pi, sin, arange, array, polyfit, polyval
    from matplotlib.pyplot import plot, show, legend

    def f(x):
        return sin(x)

    #Part 1: What the book says
    fit = mine.PolynomialFitting(num_nodes = 4, left_lim = 0, right_lim = pi / 2, step = 0.01,function = f)
    interpPoly, points, poly = fit.ChebyPoly()
    plot(fit.x_array, interpPoly)
    plot(fit.x_array, f(fit.x_array))
    show()

    #Part 2: What Canvas adds
    x_array = arange(-pi, pi+0.1, 0.1)
    xdata = array([-pi,-pi/2,0,pi/2,pi])
    ydata = f(xdata)
    fit2 = mine.PolynomialFitting(x_array = x_array, xdata = xdata, ydata = ydata)
    interpPoly2 = fit2.LagrPoly()
    plot(x_array, interpPoly2, label = "My Sine")

    polyCoeff = polyfit(xdata, ydata, len(xdata) - 1)
    plot(x_array, polyval(polyCoeff,x_array), label = "Numpy Sine")

    plot(x_array, f(x_array),label = "Actual Sine")
    legend()
    show()
