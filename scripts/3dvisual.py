import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from dcof import DCOF

class Visual:

    def __init__(self, fig:plt.figure):
        self.ax = Axes3D(fig)
    
    def setDCOF(self, dcof:DCOF):
        self.dcof = DCOF

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
