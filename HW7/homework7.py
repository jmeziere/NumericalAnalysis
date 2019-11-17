import sys
sys.path.append('../')
import methods as mine

problem = input("Test problem C2-P")

#C2-P6 a.k.a. Exercise 2.2.1ab
#Use the code fragments for Gaussian elimination in the previous section to write a 
#Python script to take a matrix A as input and output L and U. No row exchanges are
#allowed-the program should be designed to shut down if it encounters a zero pivot.
#Check your program by factoring the matrices in Exercise 2.
if problem == '6':
    Aa = [[1,2],[3,4]]
    linEquatsA = mine.SysEquats(Aa,[0,0,0])
    linEquatsA.LUFact()
    print(linEquatsA.L)
    print(linEquatsA.U)

    Ab = [[1,3],[2,2]]
    linEquatsB = mine.SysEquats(Ab,[0,0,0])
    linEquatsB.LUFact()
    print(linEquatsB.L)
    print(linEquatsB.U)
