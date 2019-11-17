import sys
sys.path.append('../')
import methods as mine

problem = input("Test problem C2-P")

#C2-P2 a.k.a. Exercise 2.1.2ac
#Let H denote the n x n Hilbert matrix, whose (i,j) entry is 1/(i+j-1). Solve Hx = b,
#where b is the vector of all ones
if problem == '2':
    def solution(n):
        from numpy import arange, meshgrid, full
        indeces = arange(1,n+1,1)
        i,j = meshgrid(indeces,indeces)
        H = 1/(i+j-1)
        b = full(n,1)
        solver = mine.SysEquats(H,b)
        return(solver.GaussElim())

    #Part a
    #n = 2
    print("\nPart a")
   
    print("The solution vector for n = 2 is\n",solution(2))

    #Part c
    #n = 10
    print("\nPart c")

    print("The solution vector for n = 10 is\n",solution(10))
