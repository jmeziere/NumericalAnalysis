import sys
sys.path.append('../')
import methods as mine

chapter = input("Chapter:")
problem = input("Problem:")

if chapter == '2':
    if problem == '12':
#C2-P12 a.k.a. Exercise 2.5.2
#Use the Jacobi method to solve the sparse system within 3 correct decimal places (forward
#error in the infinity norm) for n = 100. The correct solution is [1,-1,...,1,-1]. Report
#the number of steps needed and the backward error. The system has A matrix with 3 on the
#diagonals and 1 on the off-diagonals. The system has b vector of [1,0,0,...,0,0,-1].
#Note: By the homework we are asked to solve using both Jacobi and Gauss-Seidel,
#and compare results.
        from numpy import eye, zeros, arange

        A = 2*eye(100)+eye(100,k=1)+eye(100,k=-1)
        b = zeros(100)
        b[0] = 1
        b[-1] = -1
        solution = zeros(100)
        solution[(arange(100)%2==0)] = 1
        solution[(arange(100)%2==1)] = -1
        solver = mine.SysEquats(A,b,tol = 0.5e-3)
        solutionJaco,errorsJaco = solver.Jacobi()
        solutionGau,errorsGau = solver.GaussSeid()
        print("The solution found by the Jacobi method is",solutionJaco)
        print("The solution found by the Gauss-Seidel method is",solutionGau)
        for i in range(max(len(errorsJaco),len(errorsGau))):
            print("Step:",i,",Jacobi Error:","Jacobi has finished" if i > len(errorsJaco) else errorsJaco[i],",Gauss-Seidel Error:","Gauss-Seidel has finished" if i > len(errorsGau) else errorsGau[i])
        
if chapter == '4':
    if problem == '2':
#C4-P2 a.k.a. Exercise 4.1.5
#A company test-markets a new soft dring in 22 cities of approximately equal size. The 
#selling price (in dollars) and the number sold per week in the cities are given in the
#defined dataset.
        price = [0.59,0.8,0.95,0.45,0.79,0.99,0.9,0.65,0.79,0.69,0.79,0.49,1.09,0.95,0.79,0.65,0.45,0.6,0.89,0.79,0.99,0.85]
        sales = [3980,2200,1850,6100,2100,1700,2000,4200,2440,3300,2300,6000,1190,1960,2760,4330,6960,4160,1990,2860,1920,2160]

#Part a
#First, the company wants to find the "demand curve": how many it will sell at each
#potential price. Let P denote price and S denote sales per week. Find the line
#S = c1 + c2P that best fits the data from the table in the sense of least squares. Find
#the normal equations and the coefficients c1 and c2 of the least wquares line. Plot the
#least squares line along with the data, and calculate the root mean square error.
        from numpy import array, full_like, transpose, linspace

        price = array(price)
        sales = array(sales)
        ones_vec = full_like(price,1)
        price = transpose(array([ones_vec,price]))
        print(price)

        lsFinder = mine.LeastSquares(price,sales)
        solution,rmse = lsFinder.LS()
        
        print("The RMSE is",rmse)

        from matplotlib.pyplot import plot, show
        x = linspace(0.4,1.2,2,endpoint = True)
        plot(x,solution[0] + solution[1]*x)
        plot(price[:,1],sales,'bo')
        show()

#Part b
#After studying the results of the test marketing, the company will set a single selling
#price P throughout the country. Give an manufacturing cost of $0.23 per unit, the total
#profit (per city, per week) is S(P - 0.23) dollars. Use the results of the previous
#least squares approximation to find the selling price for which the company's profits
#will be maximized.
        
        print("We know that S =",solution[0],"+",solution[1],"P. Thus, the total profit is (",solution[0],"+",solution[1],"P)(P-0.23). Taking the derivative of this, setting the resulting equation equal to zero, and solving for P leads us to find that the maximum profit will be at a price of $0.687 per unit.")
