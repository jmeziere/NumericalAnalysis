import methods as mine

A = [[1,-4],[2,3],[2,2]]
b = [-3,15,9]
solver = mine.LeastSquares(A,b)
print(solver.QRFact())
