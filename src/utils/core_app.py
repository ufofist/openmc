#!/usr/bin/env python

import sys

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import numpy as np

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)

        # Read data from source or leakage fraction file
        self.get_file_data()

        self.setWindowTitle('Core Map Tool')
        self.main_frame = QWidget() 
        self.setCentralWidget(self.main_frame)
       
        # Create the Figure, Canvas, and Axes
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.axes = self.fig.add_subplot(111)
        
        # Create the navigation toolbar, tied to the canvas
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # Other GUI controls
        self.label_axial_level = QLabel("Axial Level:")
        self.axial_level = QComboBox()
        self.axial_level.addItems([str(i+1) for i in range(self.nz)])
        self.connect(self.axial_level, SIGNAL('currentIndexChanged(int)'),
                     self.on_draw)

        self.label_stage = QLabel("Stage:")
        self.stage = QComboBox()
        self.stage.addItems([str(i+1) for i in range(self.nstage)])
        self.connect(self.stage, SIGNAL('currentIndexChanged(int)'),
                     self.on_draw)
        
        # Grid layout at bottom
        grid = QGridLayout()
        grid.addWidget(self.label_axial_level, 0, 0)
        grid.addWidget(self.axial_level, 0, 1)
        grid.addWidget(self.label_stage, 1, 0)
        grid.addWidget(self.stage, 1, 1)
        
        # Overall layout
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(grid)
        self.main_frame.setLayout(vbox)

        self.on_draw()

    def get_file_data(self):
        # Get filename from command line
        filename = sys.argv[-1]
        fh = open(filename,'r')

        # Read size of mesh
        words = fh.readline().split()
        self.nx = int(words[0])
        self.ny = int(words[1])
        self.nz = int(words[2])
        self.nstage = int(words[3])

        # Read values
        self.values = {}
        self.maxvalue = 0.0
        for line in fh:
            words = line.split()
            i, j, k, m = [int(item) for item in words[:4]]
            val = float(words[-1])
            self.values[i,j,k,m] = val
            if val > self.maxvalue:
                self.maxvalue = val

    def on_draw(self):
        """ Redraws the figure
        """

        # Get selected axial_level and stage
        axial_level = self.axial_level.currentIndex() + 1
        stage = self.stage.currentIndex() + 1

        # Set up matrix
        matrix = np.array([[self.values[i+1,j+1,axial_level,stage]
                            for i in range(self.nx)] for j in range(self.ny)])

        # clear the figure
        self.fig.clear()

        # Make figure
        self.axes = self.fig.add_subplot(111)
        self.axes.pcolor(matrix, vmin=0.0, vmax=self.maxvalue)
        self.axes.set_xlim(0,self.nx)
        self.axes.set_ylim(0,self.ny)
        self.axes.set_xticks([])
        self.axes.set_yticks([])
        self.axes.set_aspect('equal')

        # Set up colorbar
        cax = self.axes.imshow(matrix, vmin=0.0, vmax=self.maxvalue)
        self.fig.colorbar(cax)
        
        # Draw canvas
        self.canvas.draw()
    

def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
