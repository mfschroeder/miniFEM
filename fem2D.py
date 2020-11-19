import numpy as np
from scipy.sparse import dok_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# right-hand side
f = -1.0

# quadrature points and weights
quad = simpson = [(0.0, 1 / 6), (0.5, 4 / 6), (1.0, 1 / 6)]  # simpson rule

# basis functions in 1D
hat0 = {"eval": lambda x: x, "nabla": lambda x: 1}
hat1 = {"eval": lambda x: 1 - x, "nabla": lambda x: -1}
hatFunction = [hat0, hat1]

# basis functions in 2D = tensor product of 1D basis functions
class Phi2D:
    def __init__(self, x0, y0):
        self.xBasis = hatFunction[x0]
        self.yBasis = hatFunction[y0]

    # function evaluation
    def eval(self, x, y):
        return self.xBasis["eval"](x) * self.yBasis["eval"](y)

    # derivative
    def nabla(self, x, y):
        return np.array(
            [
                self.xBasis["nabla"](x) * self.yBasis["eval"](y),
                self.xBasis["eval"](x) * self.yBasis["nabla"](y),
            ]
        )


# linear basis functions in 2D
basisFunction = {0: Phi2D(0, 0), 1: Phi2D(1, 0), 2: Phi2D(0, 1), 3: Phi2D(1, 1)}


class Cell:
    def __init__(self, origin, right, up, upRight):
        self.dofs = [origin, right, up, upRight]


class DoF:
    def __init__(self, x, y, ind=-1):
        self.x, self.y = x, y
        self.ind = ind


class Grid:
    def __init__(self, xMin=0.0, xMax=1.0, yMin=0.0, yMax=1.0, stepSize=0.5):
        self.xMin, self.xMax = xMin, xMax
        self.yMin, self.yMax = yMin, yMax
        self.h = stepSize
        self.dofs, self.cells = [], []

    def decompose(self):
        # create mesh and transform it into a list of x and y coordinates
        xRange = np.arange(self.xMin, self.xMax, self.h)
        xRange = np.append(xRange, [self.xMax])
        self.xNum = len(xRange)
        yRange = np.arange(self.yMin, self.yMax, self.h)
        yRange = np.append(yRange, [self.yMax])
        self.yNum = len(yRange)
        self.xCoord, self.yCoord = np.meshgrid(xRange, yRange)
        xList, yList = self.xCoord.ravel(), self.yCoord.ravel()

        # create DoFs
        for i, (x, y) in enumerate(zip(xList, yList)):
            self.dofs.append(DoF(x, y, ind=i))

        # create cells
        for i, dof in enumerate(self.dofs):
            if dof.x != self.xMax and dof.y != self.yMax:
                self.cells.append(
                    Cell(
                        dof,
                        self.dofs[i + 1],
                        self.dofs[i + self.xNum],
                        self.dofs[i + self.xNum + 1],
                    )
                )

    def assembleSystem(self):
        # system matrix
        A = dok_matrix((len(self.dofs), len(self.dofs)), dtype=np.float32)
        # system right hand side
        F = np.zeros(len(self.dofs), dtype=np.float32)

        for cell in self.cells:
            for x, y, quadWeight in [ (x, y, wX * wY) for (x, wX) in quad for (y, wY) in quad ]:
                for j, dof_j in enumerate(cell.dofs):
                    # assemble rhs
                    F[dof_j.ind] += (
                        quadWeight * f * basisFunction[j].eval(x, y) * self.h ** 2
                    )
                    for i, dof_i in enumerate(cell.dofs):
                        # assemble matrix
                        A[dof_i.ind, dof_j.ind] += quadWeight * basisFunction[i].nabla(x, y).dot(
                            basisFunction[j].nabla(x, y)
                        )

        # apply homogeneous Dirichlet boundary conditions
        for dof in self.dofs:
            if dof.x in [self.xMin, self.xMax] or dof.y in [self.yMin, self.yMax]:
                _, nonZeroColumns = A[dof.ind, :].nonzero()
                for j in nonZeroColumns:
                    A[dof.ind, j] = 0.0
                A[dof.ind, dof.ind] = 1.0
                F[dof.ind] = 0.0

        return A.tocsr(), F

    def plotSolution(self, solution):
        # 2D contour plot
        plt.contourf(
            self.xCoord,
            self.yCoord,
            solution.reshape(self.yNum, self.xNum),
            20,
            cmap="viridis",
        )
        plt.colorbar()
        plt.show()

    def printSolutionAtCenter(self,solution):
        #prints the solution at center of the domain
        xCenter = (self.xMax + self.xMin)/2
        yCenter = (self.yMax + self.yMin)/2
        for dof in self.dofs:
            if abs(dof.x - xCenter) < 1e-10 and abs(dof.y - yCenter) < 1e-10:
                print(f'The Solution at the center x={xCenter},y={yCenter} has the value {solution[dof.ind]}.')
                break
        else:
            print(f'The center point x={xCenter},y={yCenter} is not on the grid.')



if __name__ == "__main__":
    grid = Grid(stepSize=0.1)
    grid.decompose()
    A, F = grid.assembleSystem()
    solution = spsolve(A, F)  # U = A^{-1}F
    grid.printSolutionAtCenter(solution)
    grid.plotSolution(solution)
    
