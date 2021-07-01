# Design of Deployable Curved-Surface Rigid Origami Flashers

Author: Sen Wang, Yinghao Gao, Hailin Huang, Bing Li, Hongwei Guo, Rongqiang Liu

Affiliation: Harbin Institute of Technology

# Background

![](figures/flasher_CG.png)
Fig 1. Possible usage of deployable origami flashers for solar array on a spacecraft [1]

> [1]	S. A. Zirbel and R. J. Lang, "Accommodating Thickness in Origami-Based Deployable Arrays," Journal of Mechanical Design, vol. 135, 2013.


![](figures/CG.png)
Fig 2. Four computed models with different hubs 

This paper proposes an approach to designing rigid origami flashers that can be deployed onto curved-surface configurations. The method of designing planar origami flashers that can be wrapped around regular polygonal central hubs is presented. Based on the principle of parallel projection, planar origami flashers are projected onto target spherical surfaces to obtain the vertices on the boundary creases between sections of adjacent origami flashers. The geometric relationships of thin-panel curved origami flashers are established in terms of foldability, and other vertices in each section are calculated using numerical methods. Flexible and elastic hinges modify the thin-panel curved origami flashers into thick-panel rigid flashers. The thick rigid panels maintain the shape of flashers, and the flexible creases treated using thickness-accommodating methodologies are used to maintain the foldability of the thick panels. Additionally, we analyzed the parameters that affect the accuracy of the surface fitting. The feasibility of the proposed approach is verified using numerical simulations and physical prototypes. The novel designing of deployable curved-surface rigid origami flashers facilitates their potential applicability in solid surface antennas, surface reflectors, and other space engineering applications.

# Introduction to the Repo

<img src="figures/flow.jpg" width="400" height="500">  \
Fig 3. The flow of the design method

The core algorithm is all in [`DCOF.py`](scripts/DCOF.py) (Deployable Curve Origami Flasher).

Four demos of design DCOF with four hubs the are given:
- [three-layer DCOF with triangle hub](scripts/triangleDCOF.py)
- [four-layer DCOF with square hub](scripts/squareDCOF.py)
- [four-layer DCOF with hexagon hub](scripts/hexagonDCOF.py)
- [four-layer DCOF with octagon hub](scripts/octagonDCOF.py)


# Import computed vertices into CAD software
I am using AutoDesk Invertor, so I will take it for example to demonstrate how to import computed vertices into CAD softwares. 

Users can use `DCOF.exportVertices(filename)` to export computed vertices into an excel file (`.csv`).

![](figures/CAD.jpg)
Fig 4. Import the computed vertices into CAD softwares

![](figures/prototype.png)
Fig 5. Membrane prototype and thin-panel prototype
