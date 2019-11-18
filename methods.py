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
        if self.f(a) == 0:
            return a
        elif self.f(b) == 0:
            return b
        error = self.initial_error
        while error > self.tol:
            c = b - self.f(b) * (b - a) / (self.f(b) - self.f(a))
            if self.plot:
                guesses.append(c)
            error = abs(b - c)
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

class SysEquats:
    def __init__(self, H, b, tol = 1e-6, initGuess = 0):
        from numpy import array,float64,zeros_like
        self.H = array(H,float64)
        self.b = array(b,float64)
        self.A = self.H
        self.tol = tol
        if initGuess == 0:
            self.guess = zeros_like(b)
        else:
            self.guess = array(initGuess,float64)

    def reduceColumn(self,index,A):
        from collections import deque
        multipliers = deque()
        for i in range(index+1,len(self.H)):
            multipliers.append(A[i,index]/A[index,index])
            A[i] -= multipliers[-1] * A[index]
        return (multipliers,A)

    def swapColumns(self,index,A):
        swappers = sorted(range(len(A)-index),key = lambda k: A[index:,index][k])
        swappers.reverse()
        A[index:] = A[index:][swappers]
        return (swappers,A)

    def solveSystem(self,A,b):
        from numpy import flipud, fliplr, array
        flip_again = False
        
        if A[-1,0] == 0:
            A = flipud(A)
            b = flipud(b)
            flip_again = True
        else:
            A = fliplr(A)
        
        solution = []
        for i in range(len(A)):
            right_side = b[i]
            for j in range(i):
                right_side -= solution[j]*A[i,-j-1]
            solution.append(right_side/A[i,-i-1])
        if flip_again:
            return flipud(array(solution))
        return array(solution)


    def GaussElim(self):
        gauss_A = self.H #no sense keeping gauss_A since it only works for this one b
        b = self.b

        for i in range(len(self.H)-1):
            if gauss_A[i,i] == 0:
                raise ValueError("I can't have a zero pivot point! Now I'm broken :(")
            (multipliers,gauss_A) = self.reduceColumn(i,gauss_A)
            for j in range(i+1,len(self.H)):
                b[j] -= multipliers.popleft() * b[i]

        return self.solveSystem(gauss_A,b)

    def LUFact(self,b = None):
        if b != None:
            return self.solveSystem(self.U,self.solveSystem(self.L,b))

        from numpy import eye, copy
        self.U = copy(self.H)
        self.L = eye(len(self.H))
        b = self.b
        
        for i in range(len(self.H)-1):
            if self.U[i,i] == 0:
                raise ValueError("I can't have a zero pivot point! Now I'm borken :(")
            (multipliers,self.U) = self.reduceColumn(i,self.U)
            for j in range(i+1,len(self.H)):
                self.L[j,i] = multipliers.popleft()

        return self.solveSystem(self.U,self.solveSystem(self.L,b))
    
    def PALUFact(self,b = None):
        if b != None:
            return self.solveSystem(self.U2,self.solveSystem(self.L2,matmul(self.P,b)))

        from numpy import eye, matmul, copy
        self.U2 = copy(self.H)
        self.L2 = eye(len(self.H))
        self.P = eye(len(self.H))
        b = self.b
        
        for i in range(len(self.H)-1):
            (swappers,self.U2) = self.swapColumns(i,self.U2)
            self.P[i:] = self.P[i:][swappers]
            self.L2[i:,:i] = self.L2[i:,:i][swappers]
            (multipliers,self.U2) = self.reduceColumn(i,self.U2)
            for j in range(i+1,len(self.H)):
                self.L2[j,i] = multipliers.popleft()
        
        return self.solveSystem(self.U2,self.solveSystem(self.L2,matmul(self.P,b)))

    def Jacobi(self, guess = 0, solution = 0):
        from numpy import triu, tril, diagonal, matmul, array, multiply, float64
        for i in range(len(self.H)):
            if sum(self.H[i]) > 2*self.H[i,i]:
                print("Matrix Is NOT Diagonally Dominant! May not converge.")
                break;

        if type(guess) != int:
            self.guess = array(guess,float64)
        U = triu(self.H,k = 1)
        L = tril(self.H,-1)
        DInv = 1/diagonal(self.H)
        if type(solution) == int:
            error = max(abs(matmul(self.H,self.guess) - self.b))
        else:
            error = max(abs(solution - self.guess))
        errors = [error]
        while error > self.tol:
            self.guess = multiply(DInv,self.b - matmul(L+U,self.guess))
            if type(solution) == int:
                error = max(abs(matmul(self.H,self.guess) - self.b))
            else:
                error = max(abs(solution - self.guess))
            errors.append(error)

        return (self.guess,array(errors))

    def GaussSeid(self, guess = 0, solution = 0):
        from numpy import triu, tril, diag, matmul, array, float64
        for i in range(len(self.H)):
            if sum(self.H[i]) > 2*self.H[i,i]:
                print("Matrix Is NOT Diagonally Dominant! May not converge.")
                break;

        if type(guess) != int:
            self.guess = array(guess,float64)
        if type(solution) == int:
            error = max(abs(matmul(self.H,self.b) - self.guess))
        else:
            error = max(abs(guess - solution))
        errors = [error]
        while error > self.tol:
            for i in range(len(self.H)):
                self.guess[i] = self.b[i]
                for j in range(len(self.H)):
                    if i != j:
                        self.guess[i] -= self.H[i,j] * self.guess[j]
                c = self.guess[i] / self.H[i,i]
                self.guess[i] = c
            if type(solution) == int:
                error = max(abs(matmul(self.H, self.guess) - self.b))
            else:
                error = max(abs(guess - solution))
            errors.append(error)

        return (self.guess,array(errors))

class LeastSquares:
    def __init__(self,A,b):
        from numpy import array, float64
        self.A = array(A,float64)
        self.b = array(b,float64)
        
    def LS(self):
        from numpy import transpose, matmul, sqrt
        a_mat = matmul(transpose(self.A),self.A)
        b_mat = matmul(transpose(self.A),self.b)

        linearSolver = SysEquats(a_mat,b_mat)
        solution = linearSolver.GaussElim()

        r = self.b - matmul(self.A,solution)
        rmse = sqrt(sum(r**2))

        return (solution,rmse)

    def GramSchmidt(self, matrix, method = 'classic'):
        from numpy.linalg import norm
        from numpy import zeros_like, zeros, dot, float64

        Q = zeros_like(matrix,float64)
        R = zeros((len(matrix[0]),len(matrix[0])),float64)
        for i in range(len(Q[0])):
            q = matrix[:,i]
            for j in range(i):
                if method == 'classic' or j == 0:
                    R[j,i] = dot(Q[:,j],matrix[:,i])
                elif method == 'modern':
                    R[j,i] = dot(Q[:,j],q)
                q -= R[j,i]*Q[:,j]
            R[i,i] = norm(q)
            if R[i,i] < 1e-5:
                return (zeros(0),i)
            Q[:,i] = q/R[i,i]
        return (R,Q)
                
    def QRFact(self, gsmethod = 'classic',maxIts = 10):
        from numpy import transpose, matmul, append, expand_dims
        from numpy.random import uniform
        A = self.A
        for i in range(len(A[0]),len(A)):
            new_vec = uniform(-50,50,len(A))
            new_vec = expand_dims(new_vec,axis = 1)
            A = append(A,new_vec,1)
        Itnum = 0
        repeat = True
        Q = None
        while repeat == True and Itnum < maxIts:
            if Q != None:
                A[i] = uniform(-50,50,len(A))
            (R,Q) = self.GramSchmidt(A,gsmethod)
            if R.size == 0:
                replace_vector = Q
            else:
                repeat = False
            Itnum += 1
        new_b = matmul(transpose(Q),self.b)
        solver = SysEquats(self.A,self.b)
        solution = solver.solveSystem(R[:len(self.A[0]),:len(self.A[0])],new_b[:len(self.A[0])])
        return (solution[:len(self.A)],abs(new_b[len(self.A[0]):]))


