'''
This file holds one class and an implementation for that class to solve the Stewart ...
problem. The PDF attached describes the problem in more detail. In general, the following is the way the code is used. An implementation of the class Stewart takes in 9 parameters.
These are:
    x1 - x position of strut 2
    x2 - x position of strut 3
    y2 - y position of strut 3
    L1 - length of side 1 of the triangle
    L2 - length of side 2 of the triangle
    L3 - length of side 3 of the triangle
    p1 - length of strut 1
    p3 - length of strut 3
    gamma - angle of the triangle, sort of
Conspicuously missing are the parameters x1, y1, x2, and p2. For x1, y1, and x2, we assume these to be 0, and we vary p2 as part of the problem.

The Solve function runs the main algorithm. First, it calculates the function f(theta) over the domain and then calculates the zeros of the function, and finally plots the triangle that corresponds to that solution.

The FindIntervals function does adds in the functionality of calculating the intervals over which p2 is certain values, returning a list that contains the intervals for when there are 2,4, and 6 solutions.
'''
class Stewart:
    # The init class is where I define all of my variables, including the domain from
    # -pi to pi.
    def __init__(self, x1,x2,y2,gamma,p1,p3,L1,L2,L3):
        from numpy import linspace
        
        self.x1 = x1
        self.x2 = x2
        self.y2 = y2
        self.gamma = gamma
        self.p1 = p1
        self.p3 = p3
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.theta = linspace(-pi,pi,300,endpoint = True)
    
    # This returns the function f(theta) evaluated on the input.
    def f(self, theta, iszero = False):
        from numpy import cos, sin

        A2 = self.L3 * cos(theta) - self.x1
        B2 = self.L3 * sin(theta)
        A3 = self.L2 * cos(theta + self.gamma) - self.x2
        B3 = self.L2 * sin(theta + self.gamma) - self.y2

        N1 = B3*(self.p2**2-self.p1**2-A2**2-B2**2)-B2*(self.p3**2-self.p1**2-A3**2-B3**2)
        N2 = -A3*(self.p2**2-self.p1**2-A2**2-B2**2)+A2*(self.p3**2-self.p1**2-A3**2-B3**2)
        D = 2 * (A2 * B3 - B2 * A3)
        if iszero:
            x = N1/D
            y = N2/D
            return x,y,theta
        return N1**2 + N2**2 - self.p1**2 * D**2
    
    # This plots the function f(theta)
    def fPlot(self):
        from matplotlib.pyplot import plot, show, title
        plot(self.theta, self.fTheta)
        title("f(theta) for p = "+str(self.p2))
        show()

    # This calculates an estimate for the location of the zeros in the array
    def getEstLocZeros(self):
        from numpy import where
        return where(self.fTheta[:-1] * self.fTheta[1:] < 0)[0]

    # This finds the roots according to the secant method, which is defined in my
    # methods.py file
    def getXY(self,zero_loc):
        from methods import RootFinding
        
        # Note, if your grid isn't fine enough, you will miss zeros
        findMyRoots = RootFinding(self.f,0.5e-8,self.theta[zero_loc],self.theta[zero_loc+1])
        return self.f(findMyRoots.SecantMethod(),True);#Secant method is commented in
                                                       #methods.py

    # This plots the solution to the Stewart platform problem based off of the zeros
    # we've found.
    def plotTriangles(self, vals_array, num_zeros,sol_num):
        from matplotlib.pyplot import plot, show, title
        from numpy import sin, cos

        x = vals_array[0]
        y = vals_array[1]
        theta = vals_array[2]

        x_vals = [x, x + self.L3 * cos(theta), x + self.L2 * cos(theta + self.gamma), x]
        y_vals = [y, y + self.L3 * sin(theta), y + self.L2 * sin(theta + self.gamma), y]

        plot(x_vals,y_vals,'r-')
        plot([0,self.x1,self.x2],[0,0,self.y2],'bo')
        plot([0,x],[0,y],'b-')
        plot([self.x1,x_vals[1]],[0,y_vals[1]],'b-')
        plot([self.x2,x_vals[2]],[self.y2,y_vals[2]],'b-')
        title("For p2 = "+str(self.p2))
        show()

    # The runs the solver
    def Solve(self, p2, isPlot = True):
        self.p2 = p2

        #Find the function f(theta)
        self.fTheta = self.f(self.theta)
        if isPlot:
            #Plot the function f(theta)
            self.fPlot()
        
        #Get a good estimate for the location of the zeros in the array
        zero_locs = self.getEstLocZeros()
        if isPlot:
            for i in zero_locs:
                #Plot the triangles after finding the zeros
                self.plotTriangles(self.getXY(i),len(zero_locs),i)
        else:
            #Gives us back the number of zeros if we don't want to plot
            return len(zero_locs)

    #Gets the intervals for the number of solutions that correspond to certain values
    #of p2. I used the bisection method
    def FindIntervals(self, lower_bound, max_zeros, upper_bound):
        min_int_len = 7
        lowest = self.Solve(lower_bound, False)
        middle = self.Solve(max_zeros, False)
        highest = self.Solve(upper_bound, False)

        num_zeros = [lowest, middle,highest]
        intervals = [lower_bound,max_zeros,upper_bound]
        #We're going to go until the smallest interval is less that our tolerance.
        #This will happen when our solutions are odd, and so there is a solution only at
        #a point
        while min_int_len > 0.5e-8:
            min_len = 1
            for i in range(len(intervals)-1):
                #If the number of zeros on either side of us is the same, we are inside
                #the interval, so we continue and don't save that number
                if num_zeros[i] == num_zeros[i+1]:
                    continue
                midpoint = (intervals[i+1] + intervals[i])/2
                #Find the number of zeros for the midpoint
                midpoint_num_zeros = self.Solve(midpoint, False)
                #Here we check to see if the midpoint has the same number of zeros
                #as the left
                if midpoint_num_zeros == num_zeros[i]:
                    #We check to see if one more farther on the left side has the same
                    #number of zeros. If so, this replaces the left point as the midpoint
                    if i > 0 and midpoint_num_zeros == num_zeros[i-1]:
                        num_zeros[i] = midpoint_num_zeros
                        intervals[i] = midpoint
                    #if not, we know that this is the right side of the interval that 
                    #corresponds to the left point
                    else:
                        num_zeros.insert(i+1,midpoint_num_zeros)
                        intervals.insert(i+1,midpoint)
                #Similar logic holds here
                elif midpoint_num_zeros == num_zeros[i+1]:
                    if i < len(intervals) - 2 and midpoint_num_zeros == num_zeros[i+2]:
                        num_zeros[i+1] = midpoint_num_zeros
                        intervals[i+1] = midpoint
                    else:
                        num_zeros.insert(i+1,midpoint_num_zeros)
                        intervals.insert(i+1,midpoint)
                #If it's not equal to either side, then we know that this is a new
                #interval.
                else:
                    num_zeros.insert(i+1,midpoint_num_zeros)
                    intervals.insert(i+1,midpoint)
                
                try_min = intervals[i+1] - intervals[i]
                if try_min < min_len:
                    min_len = try_min
            min_int_len = min_len
        #clean up the intervals and num_zeros a bit
        new_intervals = [[intervals[2*i],intervals[2*i+1]] for i in range(len(intervals) // 2)]
        num_zeros = [num_zeros[2*i] for i in range(len(intervals) // 2)]
        return new_intervals, num_zeros
        
from numpy import linspace, pi, sqrt

#Part 1-3
L1 = 2
L2 = sqrt(2)
L3 = sqrt(2)
gamma = pi / 2
p1 = sqrt(5)
p2 = sqrt(5)
p3 = sqrt(5)
x1 = 4
x2 = 0
y2 = 4
LetsSolveStewart = Stewart(x1,x2,y2,gamma,p1,p3,L1,L2,L3)
LetsSolveStewart.Solve(p2) #Part 2 & 3
print(LetsSolveStewart.f(pi/4)) #Part 1
print(LetsSolveStewart.f(-pi/4)) #Part 1

#Part 4-7
L1 = 3
L2 = 3 * sqrt(2)
L3 = 3
gamma = pi / 4
p1 = 5
p2 = 5
p3 = 3
x1 = 5
x2 = 0
y2 = 6
NowSolveAgain = Stewart(x1,x2,y2,gamma,p1,p3,L1,L2,L3)
NowSolveAgain.Solve(p2) #Part 4
NowSolveAgain.Solve(7) #Part 5
intervals, num_zeros = NowSolveAgain.FindIntervals(0,7,10) #Part 6 & 7
[print("On the interval",intervals[i],"there are",num_zeros[i],"solutions.") for i in range(len(intervals))]
