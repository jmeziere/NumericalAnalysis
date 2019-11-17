class RootFinding:
    def __init__(self,f,tol,x0,x1 = 0,initial_error = 1,fp=0,plot = False):
        self.f = f
        self.tol = tol
        self.x0 = x0
        self.x1 = x1
        self.plot = plot
        self.initial_error = initial_error
        self.guesses = [x0]
        self.fp = fp

    def BisectionMethod(self):
        guesses = self.guesses.copy()
        a = self.x0
        b = self.x1
        error = abs(a - b)
        isANeg = False
        if self.f(a) == 0:
            return a
        elif self.f(b) == 0:
            return b
        elif self.f(a) < 0:
            isANeg = True
            if self.f(b) < 0:
                raise ValueError('Both f(x0) and f(x1) are negative')
        else:
            if self.f(b) > 0:
                raise ValueError('Both f(x0) and f(x1) are positive')
        while error > self.tol:
            error /= 2
            c = (a + b) / 2
            if self.plot:
                guesses.append(c)
            if self.f(c) == 0:
                return c
            elif self.f(c) < 0:
                if isANeg:
                    a = c
                else:
                    b = c
            else:
                if isANeg:
                    b = c
                else:
                    a = c
        if self.plot:
            from numpy import array
            return (array(guesses), c)
        return c

    def NewtonsMethod(self):
        if self.fp == 0:
            raise ValueError('We need a function')
        guesses = self.guesses.copy()
        guess = self.x0
        error = self.initial_error
        if self.f(guess) == 0:
            return guess
        while error > self.tol:
            new_guess = guess - self.f(guess) / self.fp(guess)
            if self.plot:
                guesses.append(new_guess)
            error = abs(new_guess - guess)
            guess = new_guess
        if self.plot:
            from numpy import array
            return (array(guesses), guess)
        return guess

    def SecantMethod(self):
        guesses = self.guesses.copy()
        a = self.x0
        b = self.x1
        if self.f(a) == 0: #We have found a root!
            return a
        elif self.f(b) == 0: #We have found a root!
            return b
        error = self.initial_error
        while error > self.tol: #Get under our tolerance
            c = b - self.f(b) * (b - a) / (self.f(b) - self.f(a)) #Our new x-guess
            if self.plot:
                guesses.append(c)
            error = abs(b - c) #Update the error
            a = b
            b = c
        if self.plot:
            from numpy import array
            return (array(guesses), b)
        return b

    def FixedPointMethod(self):
        guesses = self.guesses.copy()
        guess = self.x0
        error = self.initial_error
        if guess == self.f(guess):
            return guess
        while error > self.tol:
            new_guess = self.f(guess)
            if self.plot:
                guesses.append(new_guess)
            error = abs(new_guess - guess)
            guess = new_guess
        if self.plot:
            from numpy import array
            return (array(guesses), guess)
        return guess

class PolynomialFitting:
    def __init__(self,xdata = [],ydata = [],x_array = [],num_nodes = 0,left_lim = 0,right_lim = 0, step = 0, function = 0):
        from numpy import array, arange
        self.xdata = array(xdata).astype(float)
        self.ydata = array(ydata).astype(float)
        if len(x_array) == 0:
            self.x_array = arange(left_lim, right_lim + step, step)
        else:
            self.x_array = array(x_array).astype(float)
        self.num_nodes = num_nodes
        self.left_lim = left_lim
        self.right_lim = right_lim
        self.function = function

    def NewtDifDiv(self):
        from numpy import zeros

        row = zeros((len(self.xdata),3))
        row[:,0] = self.ydata
        row[:,1] = self.xdata
        row[:,2] = self.xdata
        coeff = [self.ydata[0]]

        while len(row) > 1:
            new_row = zeros((len(row) - 1,3))
            for i in range(len(new_row)):
                new_row[i,0] = (row[i+1,0] - row[i,0]) / (row[i+1,2] - row[i,1])
                new_row[i,1] = row[i,1]
                new_row[i,2] = row[i+1,2]
            coeff.append(new_row[0,0])
            row = new_row

        value = 0
        for i in range(len(coeff)):
            new_value = coeff[i]
            for j in range(i):
                new_value *= (self.x_array - self.xdata[j])
            value += new_value 
        return value
   
    def LagrPoly(self):
        value = 0
        for i in range(len(self.ydata)):
            new_value = self.ydata[i]
            for j in range(len(self.xdata)):
                if i != j:
                    new_value *= (self.x_array - self.xdata[j])
                    new_value /= (self.xdata[i] - self.xdata[j])
            value += new_value
        return value

    def ChebyPoly(self):
        from numpy import pi, cos, array
        points = []
        value = 1
        for i in range(1,self.num_nodes):
            points.append((self.right_lim + self.left_lim)/2 + (self.right_lim-self.left_lim) * cos((2*i-1)*pi/(2*self.num_nodes)) / 2)
            value *= (self.x_array - points[-1])
        
        tempx = self.xdata
        tempy = self.ydata
        
        self.xdata = array(points)
        self.ydata = self.function(points)
        interpPoly = self.LagrPoly()

        self.xdata = tempx
        self.ydata = tempy
        
        return (interpPoly,points,value)

class Integration:
    def __init__(self,a,b,f,n = 0):
        self.f = f
        self.a = a
        self.b = b
        self.n = n
        self.int_type = {'Trap': self.Trap, 
                         'Simpson': self.Simpson, 
                         'Midpoint': self.Midpoint,
                         'GaussQuad': self.GaussQuad}
    
    def Trap(self,a,b):
        h = b - a
        value = self.f(a) + self.f(b)
        value *= h/2
        return value

    def Simpson(self,a,b):
        h = (b - a)/2
        value = self.f(a) + self.f(b)
        value += 4*self.f(a + h)
        value *= h/3
        return value

    def Midpoint(self,a,b):
        h = b - a
        value = h*self.f(a + h/2)
        return value

    def GaussQuad(self,a,b):
        from numpy import sqrt, array

        roots = [[-1/sqrt(3),1/sqrt(3)],[-sqrt(3/5),0,sqrt(3/5)],[-sqrt((15+2*sqrt(30))/35),-sqrt((15-2*sqrt(30))/35), sqrt((15-2*sqrt(30))/35), sqrt((15+2*sqrt(30))/35)]]

        coeff = [[1,1],[5/9,8/9,5/9],[(90-5*sqrt(30))/180,(90+5*sqrt(30))/180,(90+5*sqrt(30))/180,(90-5*sqrt(30))/180]]
        new_roots = ((b - a)*array(roots[self.n-2]) + b + a)/2
        return (b-a)*sum(array(coeff[self.n - 2]) * self.f(new_roots))/2

    def Integrate(self,int_type):
        return(self.int_type[int_type](self.a,self.b))

    def Composite(self, int_type, m):
        if int_type == 'GaussQuad':
            ValueError('Composite integration not defined for Gaussian Quadrature method')

        from numpy import linspace

        points = linspace(self.a, self.b, m)
        result = 0
        for i in range(len(points)-1):
            result += self.int_type[int_type](points[i], points[i+1])
        return result

    def AdaptQuad(self, int_type, tol, a = None, b = None, result = 0):
        if int_type == 'GaussQuad' or int_type == 'Midpoint':
            ValueError('Adaptive Quadrature not defined for Gaussian Quadrature or Midpoint method')

        tol_mult = {'Trap': 3, 'Simpson': 10}
        if a == None:
            self.points = [self.a, self.b]
            self.index = 1
            a = self.a
            b = self.b
            result = self.int_type[int_type](a, b)

        left = self.int_type[int_type](a, (a+b)/2)
        right = self.int_type[int_type]((a+b)/2, b)
        self.points.append((a+b)/2)
        self.points[-1], self.points[-2] = self.points[-2], self.points[-1]
        self.index += 1
        if abs(result - (left+right)) < tol_mult[int_type] * tol:
            return left+right
        else:
            return self.AdaptQuad(int_type, tol/2, a, (a+b)/2, left) + self.AdaptQuad(int_type, tol/2, (a+b)/2, b, right)
'''
class SysEquats:
    def __init__(self, H, b):
        self.H = H
        self.b = b

    def GaussElim(self):
'''

