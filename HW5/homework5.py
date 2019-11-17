import sys
sys.path.append('../')
import methods as mine

problem = input("Test problem C5-P")

#C5-P11 a.k.a. Exercise 5.41acd
#Use Adaptive Trapezoid Quadrature to approximate the definite integral within 0.5x10^-8.
#Report the answer with eight correct decimal places and the number of subintervals required.
if problem == '11':
    #Part a
    #Int_0^4 x/sqrt(x^2+9) dx
    print("\nPart a")
    def fa(x):
        from numpy import sqrt
        return x/sqrt(x**2+9)

    integraterA = mine.Integration(0,4,fa)
    valueA = integraterA.AdaptQuad('Trap',0.5e-8)
    print("The approximate integral is",valueA,"with error",abs(2-valueA))
    print("The number of subintervals is",integraterA.index)

    #Part c
    #Int_0^1 xe^x
    print("\nPart c")
    def fc(x):
        from numpy import exp
        return x*exp(x)


    integraterC = mine.Integration(0,1,fc)
    valueC = integraterC.AdaptQuad('Trap',0.5e-8)
    print("The approximate integral is",valueC,"with error",abs(1-valueC))
    print("The number of subintervals is",integraterC.index)

    #Part d
    #Int_1^3 x^2ln(x)
    print("\nPart d")
    def fd(x):
        from numpy import log
        return x**2 * log(x)

    from numpy import log
    integraterD = mine.Integration(1,3,fd)
    valueD = integraterD.AdaptQuad('Trap',0.5e-8)
    print("The approximate integral is",valueD,"with error",abs(9*log(3)-26/9-valueD))
    print("The number of subintervals is",integraterD.index)

#C5-P12 a.k.a. 5.4.2
#Modify the MATLAB code for Adaptive Trapezoid Rule Quadrature to use Simpson's Rule
#instead, applying the criterion (5.42) with the 15 replaced by 10. Approximate the
#integral in Example 5.12 within 0.005, and compare with Figure 5.5(b). How many
#subintervals were required?
if problem == '12':
    def f(x):
        from numpy import sin, exp
        return 1 + sin(exp(3*x))

    integrater = mine.Integration(-1,1,f)
    value = integrater.AdaptQuad('Simpson',0.005)
    print("The approximate integral is",value,"with error",abs(2.500809110336424-value))
    print("The number of subintervals is",integrater.index)

#C5-P13 a.k.a. 5.4.3acd
#Carry out the steps of Computer Problem 1 for adaptive Simpson's Rule, developed in 
#Computer problem 2
if problem =='13':
    #Part a
    #Int_0^4 x/sqrt(x^2+9) dx
    print("\nPart a")
    def fa(x):
        from numpy import sqrt
        return x/sqrt(x**2+9)

    integraterA = mine.Integration(0,4,fa)
    valueA = integraterA.AdaptQuad('Simpson',0.5e-8)
    print("The approximate integral is",valueA,"with error",abs(2-valueA))
    print("The number of subintervals is",integraterA.index)

    #Part c
    #Int_0^1 xe^x
    print("\nPart c")
    def fc(x):
        from numpy import exp
        return x*exp(x)
    
    
    integraterC = mine.Integration(0,1,fc)
    valueC = integraterC.AdaptQuad('Simpson',0.5e-8)
    print("The approximate integral is",valueC,"with error",abs(1-valueC))
    print("The number of subintervals is",integraterC.index)
    
    #Part d
    #Int_1^3 x^2ln(x)
    print("\nPart d")
    def fd(x):
        from numpy import log
        return x**2 * log(x)
    
    from numpy import log
    integraterD = mine.Integration(1,3,fd)
    valueD = integraterD.AdaptQuad('Simpson',0.5e-8)
    print("The approximate integral is",valueD,"with error",abs(9*log(3)-26/9-valueD))
    print("The number of subintervals is",integraterD.index)
