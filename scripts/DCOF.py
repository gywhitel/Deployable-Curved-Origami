import numpy as np
from function import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as ptColor
from copy import deepcopy


PI = np.pi
INF = np.inf
S3 = np.sqrt(3)
S3_2 = S3 / 2
DEG = 180 / PI

class DCOF:
    '''
    The class for design of deployable curved origami flashers (DCOF)
    '''
    # self.vertex
    def __init__(self, R:float, A:float):
        '''
        Specify the radius R of the spherical surface  and the edge length A of the central hub.  \\
        The method is not only capable of designing flashers on a spherical surface, but also can be exploited for flashers on other shapes by user-specified surface equation.
        Parameters
        ---
        param : type
            introduction
        Return
        ---
        return instruction
        '''
        self.__R = R
        self.__A = A
        self.hub = []
        self.LValley = []
        self.ridge = []
        self.RValley = []

    def getRadius(self):
        return self.__R
    
    def setRadius(self, R:float):
        self.__R = R
    
    def getLength(self):
        return self.__A

    def setLength(self, A):
        self.__A = A

    def computeZ(self, p:np.array)->np.array:
        '''
        brief
        Parameters
        ---
        param : type
            introduction
        Return
        ---
        return instruction
        '''
        R = self.__R
        x, y = p[0], p[1]
        z = -np.sqrt(R*R - x*x - y*y)
        return np.array([x, y, z])

    def setSearchParameter(self, step:float, limit:float, angleTolerence:float, distanceTolerence:float):
        '''
        brief
        Parameters
        ---
        param : type
            introduction
        Return
        ---
        return instruction
        '''
        self.step = step
        self.limit = limit
        self.angleTolerence = angleTolerence
        self.distanceTolerence = distanceTolerence

    def searchVertex(self):
        '''
        Compute a ridge vertex that satisfies both surface equation constraints and folding constraints, using numerical method.
        Do 3D brutal force search in the neighborhood of initialGuess.
        Parameters
        ---
        initialGuess : np.array([x, y, z]) 3D point
            introduction
        lastVertex :  np.array([x, y, z]) 3D point
            last ridge vertex (the previous vertex we solved)
        LVVertex : np.array([x, y, z]) 3D point]
            the vertex on the left valley
        RVVertex : np.array([x, y, z]) 3D point
            the vertex on the right valley
        last_LV : np.array([x, y, z]) 3D point
            the last vertex on the left valley
        last_RV : np.array([x, y, z]) 3D point
            the last vertex on the right valley
        limit : float
            this determines the size of searching region

        Return
        vertex : np.array([x, y, z]) 3D point
        ---
        the vertex on the ridge
        '''

        # 3D brutal force search
        cost = INF
        limit = self.limit
        step = self.step
        angleTolerence, distanceTolerence = self.angleTolerence, self.distanceTolerence
        A, R = self.__A, self.__R
        TILT = angle(self.ridge[0], self.hub[0], self.RValley[0])

        sol = np.zeros(3)

        initialGuess, lastVertex = self.ridge[-1], self.ridge[-2]
        LVVertex, last_LV = self.LValley[-1], self.LValley[-2]
        RVVertex, last_RV = self.RValley[-1], self.RValley[-2]

        for x in np.arange(-limit, limit, step):
            for y in np.arange(-limit, limit, step):
                for z in np.arange(-limit, limit, step):
                    
                    vertex = initialGuess + np.array([x, y, z])     # to be solved
                # angular constraint 1
                    gamma1 = angle(vertex, lastVertex, RVVertex)
                    gamma2 = angle(vertex, lastVertex, LVVertex)
                    delta1 = angle(RVVertex, lastVertex, last_RV)
                    delta2 = angle(LVVertex, lastVertex, last_LV)
                    constraint1 = abs(gamma1 + delta1 - gamma2 - delta2)

                    if constraint1 * DEG > angleTolerence:
                        continue
                    else:
                    # angular constraint 2
                        theta1 = angle(lastVertex, vertex, RVVertex)
                        theta2 = angle(lastVertex, vertex, LVVertex)
                        constraint2 = abs(theta1 - theta2)
                        
                        if constraint2 * DEG > angleTolerence:
                            continue
                        else:
                            # length constraint
                            constraint3 = abs(norm(vertex - lastVertex) * np.cos(TILT) - A)
                            if constraint3 > distanceTolerence:
                                continue
                            else:
                            # points that get here satisfy the geometric constaints of folding
                            # Find the point with smallest fitting error from these feasible solutions
                                fitErr = abs(norm(vertex) - R)
                                if fitErr < cost:
                                    cost = fitErr
                                    
                                # if the fitting error is small enough, we consider the vertex is the solution
                                    if cost < 0.001:
                                        print("Found optimal solution:", vertex)
                                        return vertex
                                # if none of the feasible solution satisfies the fitting error tolerance, return the one with smallest fitting error
                                    sol = vertex
        if cost is INF:
            # print("Fail to find a solution, return the initial guess")
            # return None
            raise ValueError("Fail to find a solution. Suggest tuning searching parameters")

        print("Found a solution: ", sol, "Fitting error: ", cost)
        return sol

    def circularArray(self, angle:float):
        self.hub = rotateZ(self.hub, angle)
        self.LValley = rotateZ(self.LValley, angle)
        self.RValley = rotateZ(self.RValley, angle)
        self.ridge = rotateZ(self.ridge, angle)
    
class TriangleDCOF(DCOF):
    
    def __init__(self, R:float, A:float):
        DCOF.__init__(self, R, A)
    
    def computeHub(self):
        A = self.getLength()
        self.hub = np.array([
            self.computeZ([0, A / S3]),
            self.computeZ([A/2, -A/2/S3]),
            self.computeZ([-A/2, -A/2/S3]),
        ])

    def getHubOrigin(self):
        '''
        The origin is located at the centroid of the hub
        Parameters
        ---
        param : type
            introduction
        Return
        ---
        return instruction
        '''
        return (self.hub[0] + self.hub[1] + self.hub[2]) / 3

    def computeLayer(self):
        A = self.getLength()
        layer = len(self.LValley) + 1

        LValleyVertex = self.hub[0][:2] + np.array([A/2, S3 * A/2]) * layer
        LValleyVertex = self.computeZ(LValleyVertex)
        self.LValley.append(LValleyVertex)

        RValleyVertex = self.hub[0][:2] + np.array([A/2, -S3*A/2]) * layer
        RValleyVertex = self.computeZ(RValleyVertex)
        self.RValley.append(RValleyVertex)

        ridgeVertex = self.hub[0][:2] + np.array([2*A ,0]) * layer
        ridgeVertex = self.computeZ(ridgeVertex)    # initial guess
        self.ridge.append(ridgeVertex)
        if len(self.LValley) > 1:
            # The projected ridge vertex of the first layer yields little error, so there is no need to do search on the first layer
            self.ridge[-1] = self.searchVertex()
        print("Layer NO.", len(self.ridge), "is generated.")

class SquareDCOF(DCOF):
    
    def __init__(self, R:float, A:float):
        DCOF.__init__(self, R, A)
    
    def computeHub(self):
        A = self.getLength()
        self.hub = np.array([
            self.computeZ([-A/2, A/2]),
            self.computeZ([A/2, A/2]),
            self.computeZ([A/2, -A/2]),
            self.computeZ([-A/2, -A/2])
        ])

    def getHubOrigin(self):
        return (self.hub[0] + self.hub[1] + self.hub[2] + self.hub[3]) / 4

    def computeLayer(self):
        A = self.getLength()
        layer = len(self.LValley) + 1

        LValleyVertex = self.hub[0][:2] + np.array([0, A]) * layer
        LValleyVertex = self.computeZ(LValleyVertex)
        self.LValley.append(LValleyVertex)

        RValleyVertex = self.hub[0][:2] + np.array([A, 0]) * layer
        RValleyVertex = self.computeZ(RValleyVertex)
        self.RValley.append(RValleyVertex)

        ridgeVertex = self.hub[0][:2] + np.array([A, A]) * layer
        ridgeVertex = self.computeZ(ridgeVertex)    # initial guess
        self.ridge.append(ridgeVertex)
        if len(self.LValley) > 1:
            # The projected ridge vertex of the first layer yields little error, so there is no need to do search on the first layer
            self.ridge[-1] = self.searchVertex()
        print("Layer NO.", len(self.ridge), "is generated.")

class HexagonDCOF(DCOF):
    
    def __init__(self, R:float, A:float):
        DCOF.__init__(self, R, A)

    def computeHub(self):
        A = self.getLength()
        self.hub = np.array([
                    self.computeZ([-A * S3_2, A / 2]),
                    self.computeZ([0,  A]),
                    self.computeZ([A * S3_2,  A / 2]),
                    self.computeZ([A * S3_2, -A / 2]),
                    self.computeZ([0, -A]),
                    self.computeZ([-A * S3_2, -A / 2])
                    ])
    
    def getHubOrigin(self):
        return (self.hub[0] + self.hub[3]) / 2

    def computeLayer(self):
        '''
        Compute a layer of a sector of the plannar flasher
        Parameters
        ---
        param : type
            introduction
        Return
        ---
        return instruction
        '''
        A, R = self.getLength(), self.getRadius()
        layer = len(self.LValley) + 1
        LValleyVertex = self.hub[0][:2] + np.array([0, A]) * layer
        LValleyVertex = self.computeZ(LValleyVertex)
        self.LValley.append(LValleyVertex)

        RValleyVertex = self.hub[0][:2] + np.array([A * S3_2, A/2]) * layer
        RValleyVertex = self.computeZ(RValleyVertex)
        self.RValley.append(RValleyVertex)

        ridgeVertex = self.hub[0][:2] + np.array([A / S3, A]) * layer
        ridgeVertex = self.computeZ(ridgeVertex)    # initial guess
        self.ridge.append(ridgeVertex)
        if len(self.LValley) > 1:
            # The projected ridge vertex of the first layer yields little error, so there is no need to do search on the first layer
            self.ridge[-1] = self.searchVertex()
        print("Layer NO.", len(self.ridge), "is generated.")

 
        



class OctagonDCOF(DCOF):
    
    def __init__(self, R:float, A:float):
        DCOF.__init__(self, R, A)
    
    def computeHub(self):
        A = self.getLength()
        T8 = np.tan(PI/8)
        self.hub = np.array([
                    self.computeZ([-A/2/T8, A / 2]),
                    self.computeZ([-A/2,  A/2/T8]),
                    self.computeZ([A/2,  A/2/T8]),
                    self.computeZ([A/2/T8, A / 2]),
                    self.computeZ([A/2/T8, -A / 2]),
                    self.computeZ([A/2,  -A/2/T8]),
                    self.computeZ([-A/2,  -A/2/T8]),
                    self.computeZ([-A/2/T8, -A / 2]),
                    ])

    def getHubOrigin(self):
        return (self.hub[0] + self.hub[4]) / 2

    def computeLayer(self):
        A = self.getLength()
        layer = len(self.LValley) + 1

        LValleyVertex = self.hub[0][:2] + np.array([0, A]) * layer
        LValleyVertex = self.computeZ(LValleyVertex)
        self.LValley.append(LValleyVertex)

        S2 = np.sqrt(2)
        RValleyVertex = self.hub[0][:2] + np.array([A/S2, A/S2]) * layer
        RValleyVertex = self.computeZ(RValleyVertex)
        self.RValley.append(RValleyVertex)

        T3_8 = np.tan(3*PI/8)
        ridgeVertex = self.hub[0][:2] + np.array([A/T3_8, A]) * layer
        ridgeVertex = self.computeZ(ridgeVertex)    # initial guess
        self.ridge.append(ridgeVertex)
        if len(self.LValley) > 1:
            # The projected ridge vertex of the first layer yields little error, so there is no need to do search on the first layer
            self.ridge[-1] = self.searchVertex()
        print("Layer NO.", len(self.ridge), "is generated.")



class Visual:

    def __init__(self, fig:plt.figure):
        self.ax = Axes3D(fig)
    
    def drawSphere(self, dcof:DCOF):
        R = dcof.getRadius()
        xgrid = np.arange(-5, 6, 0.5)
        ygrid = np.arange(-5, 6, 0.5)
        x, y = np.meshgrid(xgrid, ygrid)
        z = -np.real(np.sqrt(R*R - x*x - y*y))
        # self.ax.plot_wireframe(x, y, z, color='lightcoral')
        self.ax.contour3D(x, y, z, 80, cmap='binary')
        self.ax.set_xlabel("X (m)")
        self.ax.set_ylabel("Y (m)")
        self.ax.set_zlabel("Z (m)")

    
    def drawHub(self, dcof:DCOF, origin:np.array):

        color = ptColor.rgb2hex(np.random.rand(3))
        hub = dcof.hub
        for i in range(-1, len(hub)-1):
            tri = Poly3DCollection([hub[i], hub[i+1], origin], facecolor=color)
            self.ax.add_collection3d(tri)


    def drawTriangle(self, p1, p2, p3, color):
        triangle = Poly3DCollection([p1, p2, p3])
        triangle.set_color(color)
        triangle.set_edgecolor('k')
        self.ax.add_collection3d(triangle)
    
    def drawLayer(self, dcof:DCOF, layer:int):
        color = ptColor.rgb2hex(np.random.rand(3))
        ridge = dcof.ridge
        LV, RV = dcof.LValley, dcof.RValley
        if layer is 1:
            origin = dcof.hub[0]
            self.drawTriangle(ridge[0], RV[0], origin, color)
            self.drawTriangle(ridge[0], LV[0], origin, color)
        else:
            current = layer-1  # 2nd layer corresponds to ridge[1]
            last = layer - 2
            self.drawTriangle(ridge[current], RV[current], ridge[last], color)  # upper right facet
            self.drawTriangle(ridge[last], RV[current], RV[last], color)    # lower right facet
            self.drawTriangle(ridge[current], LV[current], ridge[last], color)      
            self.drawTriangle(LV[current], LV[last], ridge[last], color)

    def centroSymmetry(self, dcof:DCOF):

        if isinstance(dcof, TriangleDCOF):
            self.rotation = -2*PI/3
            for i in range(1,3):
                dcof2 = deepcopy(dcof)
                dcof2.circularArray(self.rotation * i)
                for layer in range(len(dcof2.LValley)):
                    self.drawLayer(dcof2, layer+1)

        if isinstance(dcof, SquareDCOF):
            self.rotation = -PI/2
            for i in range(1,4):
                dcof2 = deepcopy(dcof)
                dcof2.circularArray(self.rotation * i)
                for layer in range(len(dcof2.LValley)):
                    self.drawLayer(dcof2, layer+1)

        if isinstance(dcof, HexagonDCOF):
            self.rotation = -PI / 3
            for i in range(1,6):
                dcof2 = deepcopy(dcof)
                dcof2.circularArray(self.rotation * i)
                for layer in range(len(dcof2.LValley)):
                    self.drawLayer(dcof2, layer+1)

        if isinstance(dcof, OctagonDCOF):
            self.rotation = -PI / 4
            for i in range(1,8):
                dcof2 = deepcopy(dcof)
                dcof2.circularArray(self.rotation * i)
                for layer in range(len(dcof2.LValley)):
                    self.drawLayer(dcof2, layer+1)
        
        

def main():
    figure = plt.figure()
    visual = Visual(figure)

    hdcof = OctagonDCOF(10, 1)
    hdcof.computeHub()
    
    visual.drawSphere(hdcof)
    hubOrigin = hdcof.getHubOrigin()
    visual.drawHub(hdcof, hubOrigin)

    hdcof.computeLayer()    # generate first layer
    visual.drawLayer(hdcof, 1)
    hdcof.setSearchParameter(0.01, 0.15, 5, 0.15)
    hdcof.computeLayer()    # generate second layer
    visual.drawLayer(hdcof ,2)
    hdcof.computeLayer()
    visual.drawLayer(hdcof, 3)
    # hdcof.computeLayer()
    # visual.drawLayer(hdcof, 4)
    visual.centroSymmetry(hdcof)

    plt.show()

if __name__ == '__main__':
    main()