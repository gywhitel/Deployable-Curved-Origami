'''
@ Author: Yinghao Gao

This script demonstrate the computation of a four-layer deployable origami array that can deploy on a spherical surface.
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as ptColor 
from function import *


# --------------------- define constants --------------------------
PI = np.pi
INF = np.inf
S3 = np.sqrt(3)
S3_2 = S3 / 2
DEG = 180 / PI

A = 1   # the edge length of the hub
R = 10  # the radius of the spherical surface

angleTolerence = 5          # slack of angular constraints
distanceTolerence = 0.15    # slack of length constraints
step = 0.01         # search step

# -----------------------------------------------------------------

# -------------------- projection function ------------------------------
def computeZ(p:np.array)->np.array:
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
    x, y = p[0], p[1]
    z = -np.sqrt(R*R - x*x - y*y)
    return np.array([x, y, z])

# -------------------- Compute vertices of first layer ----------------------
# compute 2D coordinates of vertices of plannar array
# First layer
p001 = np.array([-A * S3_2, A/2])
p002 = np.array([0,     A])
p111 = np.array([-A * S3_2, 3*A/2])
p121 = p111 + np.array([A/S3, 0])

p002 = [0, A]
p003 = [A * S3_2, A / 2]
p004 = [A * S3_2, -A / 2]
p005 = [0, -A]
p006 = [-A * S3_2, -A / 2]

# project onto the spherical surface
p001 = computeZ(p001)
p002 = computeZ(p002)
p111 = computeZ(p111)
p121 = computeZ(p121)

TILT = angle(p121, p001, p002)

p002 = computeZ(p002)
p003 = computeZ(p003)
p004 = computeZ(p004)
p005 = computeZ(p005)
p006 = computeZ(p006)



def searchVertex(initialGuess, lastVertex, LVVertex, RVVertex, last_LV, last_RV, limit:float=0.1):
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
    sol = np.zeros(3)
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
                                    return vertex
                            # if none of the feasible solution satisfies the fitting error tolerance, return the one with smallest fitting error
                                sol = vertex
    if cost is INF:
        print("Fail to find a solution")
        return None

    print("Found a solution: ", sol, "Fitting error: ", cost)
    return sol

class Visual:

    def __init__(self, fig:plt.figure):
        self.ax = Axes3D(fig)
    
    def drawSphere(self):
        xgrid = np.arange(-6, 7, 0.5)
        ygrid = np.arange(-6, 7, 0.5)
        x, y = np.meshgrid(xgrid, ygrid)
        z = -np.real(np.sqrt(R*R - x*x - y*y))
        # x[z>=0] = np.nan
        # y[z>=0] = np.nan
        # self.ax.plot_wireframe(x, y, z, color='lightcoral')
        self.ax.contour3D(x, y, z, 80, cmap='binary')
        # self.ax.set_alpha(0.1)
        self.ax.set_xlabel("X (m)")
        self.ax.set_ylabel("Y (m)")
        self.ax.set_zlabel("Z (m)")

    
    def drawHub(self):
        origin = (p001 + p004)/2

        color = ptColor.rgb2hex(np.random.rand(3))

        hexagon1 = Poly3DCollection([p001, p002, origin], facecolor=color)
        hexagon2 = Poly3DCollection([p002, p003, origin], facecolor=color)
        hexagon3 = Poly3DCollection([p003, p004, origin], facecolor=color)
        hexagon4 = Poly3DCollection([p004, p005, origin], facecolor=color)
        hexagon5 = Poly3DCollection([p005, p006, origin], facecolor=color)
        hexagon6 = Poly3DCollection([p006, p001, origin], facecolor=color)

        self.ax.add_collection3d(hexagon1)
        self.ax.add_collection3d(hexagon2)
        self.ax.add_collection3d(hexagon3)
        self.ax.add_collection3d(hexagon4)
        self.ax.add_collection3d(hexagon5)
        self.ax.add_collection3d(hexagon6)


    def drawTriangle(self, p1, p2, p3, color):
        triangle = Poly3DCollection([p1, p2, p3])
        triangle.set_color(color)
        triangle.set_edgecolor('k')
        self.ax.add_collection3d(triangle)


def main():
    fig = plt.figure()
    
    # Second layer
    p211 = p111 + p111 - p001
    p221 = p121 + p121 - p001
    p112 = p002 + p002 - p001

    p211 = computeZ(p211)
    p221 = computeZ(p221)
    p112 = computeZ(p112)
    visual = Visual(fig)
    visual.drawSphere()
    visual.drawHub()
    
    color1 = ptColor.rgb2hex(np.random.rand(3))
    visual.drawTriangle(p001, p002, p121, color1)
    visual.drawTriangle(p001, p111, p121, color1)

    # solve p221
    p221 = searchVertex(p221, p121, p211, p112, p111, p002, 0.2)
    color2 = ptColor.rgb2hex(np.random.rand(3))
    visual.drawTriangle(p112, p221, p121, color2)
    visual.drawTriangle(p112, p121, p002, color2)
    visual.drawTriangle(p221, p121, p211, color2)
    visual.drawTriangle(p211, p111, p121, color2)

    # baseVertex = p001
    ridge = [p121, p221]
    LValley = [p111, p211]
    Rvalley = [p002, p112]

    # project vertices of plannar array onto spherical surface
    for i in range(2,5):
        ridgeVertex = np.array(p121[:,2]) + i * np.array([A/S3, A])
        ridge.append(computeZ(ridgeVertex))

        LValleyVertex = np.array([-A * S3_2, 3*A/2 + i * A])
        LValley.append(computeZ(LValleyVertex))

        RValleyVertex = np.array([np.array([0, A]) + i*np.array(A * S3_2, A/2)])
        Rvalley.append(computeZ(RValleyVertex))

    for i in range(2,5):



    plt.show()
    



if __name__ == '__main__':
    main()
