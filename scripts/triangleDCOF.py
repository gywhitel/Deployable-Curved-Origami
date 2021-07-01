import DCOF
import matplotlib.pyplot as plt

if __name__ == '__main__':
    figure = plt.figure()
    visual = DCOF.Visual(figure)

    hdcof = DCOF.TriangleDCOF(10, 1)
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