import methods as mine

A = [[1,-4],[2,3],[2,2]]
b = [1,1,1]
solver = mine.LeastSquares(A,b)
solver.GrammSchmidt()
