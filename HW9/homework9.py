import sys
sys.path.append('../')
import methods as mine

chapter = input("Chapter:")
problem = input("Problem:")

if chapter == '4':
    if problem == '6':
#C4-P6 a.k.a. Ex. 4.3.4(a)
#write a program that implements classical Gram-Schmidt to find the full QR factorization.
#Check your work by comparing factorizations of the matrices in exercise 4.3.2.
        from numpy import array,float64
        matrix_a = [[2,3],[-2,-6],[1,0]]
        rand_b = [1,1,1]
        solverA = mine.LeastSquares(matrix_a,rand_b)
        (R,Q) = solverA.GramSchmidt(array(matrix_a,float64))
        print(R)
        print(Q)

        matrix_b = [[-4,-4],[-2,7],[4,-5]]
        solverB = mine.LeastSquares(matrix_b,rand_b)
        (R,Q) = solverB.GramSchmidt(array(matrix_b,float64))
        print(R)
        print(Q)
