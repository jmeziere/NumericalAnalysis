# -*- coding: utf-8 -*-
#Problem 2

#Use both np.linspace() and np.arange() to create a list of floating point num-
#bers starting at 1.0, ending at 4.0, equally spaced and separated by 0.2.  
#In other words, the output should be 1.0, 1.2, 1.4, ... , 3.8, 4.0

from numpy import linspace, arange

print(linspace(1.0,4.0,16))
print(arange(1.0,4.1,0.2))

#Problem 3

#Create  an  array  consisting  of  the  floats  1.0,2.0,3.0,4.0,5.0.
#Create  a  second  array that contains the square root of these numbers.  
#Use a for loop to subtract the arrays component-wise, square each difference, 
#and add-up the total.  In other words, compute the sum of the squared differences 
#using a for loop.  Note:  this can be done without using a for loop.

from numpy import sqrt, sum

#hard way
p2_array = arange(1.0,5.1,1.0)
p2_array_2 = sqrt(p2_array)
p2_sum = 0
for i in range(len(p2_array_2)):
    p2_sum += (p2_array_2[i] - p2_array[i])**2

#easy way
p2_sum = sum((sqrt(p2_array) - p2_array)**2)

#Problem 4

#Start with x = 1. Use a while loop to divide x by 2 until x < 1e-4.
#Display the list 1.0,0.5,0.25, ....  Also, determine k, the number of divisions 
#needed so that the (k−1)th division produces x > 10e−4 and the kth division produces
#x < 1e−4.

x = 1
counter = -1
while x >= 1e-4:
    print(x)
    x /= 2
    print(x)
    counter += 1
print('k is '+ str(counter))

#Problem 5

#Create a function whose input is a number x and whose output is e^−x cos(x), or 
#f(x) = e^−x cos(x).  Evaluate f(x) at the points 0, 0.1, 0.2, ... , 1.0

from numpy import exp, cos

def f(x):
    return exp(-x) * cos(x)

print(f(arange(0,1.1,0.1)))

#Problem 6

#Plot the function h(x) = e^x cos^2(x) − 2 on the interval − 0.5 to 5.5 and 
#visually estimate the roots of h(x) on that interval.

from matplotlib.pyplot import plot, show

xdata = linspace(0.5,5.5,100,endpoint = True)
ydata = exp(xdata) * cos(xdata)**2 - 2

plot(xdata,ydata)
show()
