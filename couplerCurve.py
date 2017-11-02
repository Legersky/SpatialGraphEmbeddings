from PyQt5 import QtCore
from PyQt5.QtWidgets import QAction,  QApplication, QInputDialog, QFileDialog, QErrorMessage, QAbstractSpinBox, QPlainTextEdit, QDialog, QVBoxLayout, QSizePolicy, QDockWidget
from PyQt5.QtGui import QIcon,  QKeySequence,  QTextCursor

from PyQt5.uic import loadUiType

import pickle
import ast
#import time 
import copy

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_agg import FigureCanvasAgg


from algRealEmbeddings import *
from graphCouplerCurve import *
from axel_vis import *

UI_MainWindow, MainWindow = loadUiType("couplerCurve.ui")
class MplWindow(UI_MainWindow, MainWindow):

    def __init__(self):
        super(MplWindow, self).__init__()
        self.verbose = 1
        self.setupUi(self)
        
        self.showMaximized()
        
        self.fig = Figure()
        
        self.graph_sequence = []
        self.graph_sequence_comments = []
        
        self._V6fromPHC = []
#        self._log = ''
#        self.setRequiresRecomputing()
        self.setActiveGraph(GraphCouplerCurve( window=self))
        self.init_widgets()
        self.init_mpl()
#        self.init_data()
        self.init_plot()
        #self.init_pick()
        self.dockLog.resize(self.dockLog.minimumWidth()*2, 300)
        self.dockLog.updateGeometry()

    
    def setActiveGraph(self, G):
        self.graph = G
        self.updateParameter()
        if self.graph.isComputed():
            self.setComputed()
        else:
            self.setRequiresRecomputing()

    def init_widgets(self):
        """Initialize the widgets, mainly connecting the relevant signals"""
        self.plainTextEdit.setReadOnly(True)
#        self.plainTextEdit.setCenterOnScroll(True)
#        self.plainTextEdit.setMaximumBlockCount(40)
        
        self.doubleSpinBoxLengths = {'12': self.doubleSpinBoxL12, 
                                    '13': self.doubleSpinBoxL13,
                                    '14': self.doubleSpinBoxL14,
                                    '15': self.doubleSpinBoxL15,
                                    '16': self.doubleSpinBoxL16, 
                                    '27': self.doubleSpinBoxL27,
                                    '37': self.doubleSpinBoxL37,
                                    '47': self.doubleSpinBoxL47,
                                    '57': self.doubleSpinBoxL57,
                                    '67': self.doubleSpinBoxL67,
                                    '23': self.doubleSpinBoxL23,
                                    '34': self.doubleSpinBoxL34,
                                    '45': self.doubleSpinBoxL45,
                                    '56': self.doubleSpinBoxL56,
                                    }

        self.update_graph2tabLengths()
        self.update_graph2phi()
        
        for e in self.doubleSpinBoxLengths:
            self.doubleSpinBoxLengths[e].valueChanged.connect(self.update_tabLengths2graph)
        
        self.update_graph2R26()
        self.doubleSpinBoxR26.valueChanged.connect(self.update_R26toGraph)
        self.doubleSpinBoxY2.valueChanged.connect(self.update_yV2toGraph)
        
        self.doubleSpinBoxPhi.valueChanged.connect(self.update_phi2graph)
        self.doubleSpinBoxTheta.valueChanged.connect(self.update_theta2graph)
        
        
        self.pushButtonPlot.clicked.connect(self.computeCouplerCurves)
        self.pushButtonPlot.setShortcut(QKeySequence("Ctrl+P"))
        
        self.checkBoxBranches = {
                                'orange': self.checkBoxOrange, 
                                'red': self.checkBoxRed, 
                                'green': self.checkBoxGreen, 
                                'blue': self.checkBoxBlue
                                }
        for c in self.checkBoxBranches:
            self.checkBoxBranches[c].stateChanged.connect(self.plotScene)
        
        self.checkBoxMirror.stateChanged.connect(self.plotScene)
        
#        self.spinBoxSamples.valueChanged.connect(self.updateSamples)
        
        self.comboBoxGraphFor.currentTextChanged.connect(self.updateDisplayedGraph)
        self.spinBoxParameter.valueChanged.connect(self.updateDisplayedGraph)
        self.spinBoxParameter.setWrapping(True)
        
        self.spinBoxParameterStep.valueChanged.connect(self.update_parameterStep)
        
        self.doc_tb = self.addToolBar("File")

        self.loadLengthsButton = QAction(QIcon.fromTheme("document-open"),"Load lengths",self)
        self.doc_tb.addAction(self.loadLengthsButton)
        self.loadLengthsButton.triggered.connect(self.loadLengths)
        action = self.menuFile.addAction('Load lengths')
        action.triggered.connect(self.loadLengths)
        
        self.saveLengthsButton = QAction(QIcon.fromTheme("document-save"),"Save lengths",self)
        self.doc_tb.addAction(self.saveLengthsButton)
        self.saveLengthsButton.triggered.connect(self.saveLengths)
        action = self.menuFile.addAction('Save lengths')
        action.triggered.connect(self.saveLengths)
        
        self.inOut_tb = self.addToolBar("Input/Output")
        
        self.insertLengthsButton = QAction(QIcon.fromTheme("insert-text"),"Insert lengths",self)
        self.inOut_tb.addAction(self.insertLengthsButton)
        self.insertLengthsButton.triggered.connect(self.insertLengths)
        action = self.menuInputOutput.addAction('Insert lengths')
        action.triggered.connect(self.insertLengths)
        
        self.insertLengthsByEmbeddingsButton = QAction("Insert lengths by embedding",self)
        self.inOut_tb.addAction(self.insertLengthsByEmbeddingsButton)
        self.insertLengthsByEmbeddingsButton.triggered.connect(self.insertLengthsByEmbeddings)
        action = self.menuInputOutput.addAction('Insert lengths by embedding')
        action.triggered.connect(self.insertLengthsByEmbeddings) 
        
        self.exportLengthsButton = QAction("Export lengths",self)
        self.inOut_tb.addAction(self.exportLengthsButton)
        self.exportLengthsButton.triggered.connect(self.exportLengths)
        action = self.menuInputOutput.addAction('Export lengths')
        action.triggered.connect(self.exportLengths)

        
        
        
        self.actionFullscreen.toggled.connect(self.fullScreen)
        
        self.updateParameter()
        
        self.buttonLoadSequence.clicked.connect(self.loadSequence)
        
        self.spinBoxImgInSeq.setMinimum(1)
        self.spinBoxImgInSeq.setMaximum(1)
        self.spinBoxImgInSeq.valueChanged.connect(self.plotGraphFromSequence)
        
        self.spinBoxSeqLength.setReadOnly(True)
        self.spinBoxSeqLength.setButtonSymbols(QAbstractSpinBox.NoButtons)
        
        self.noteToImgInSeq.setReadOnly(True)
        self.boxInfoImgInSeq.setVisible(False)
    
        self.buttonRunPHC.clicked.connect(self.runPHC)
        self.spinBoxNumberReal.setReadOnly(True)
        self.spinBoxNumberReal.setButtonSymbols(QAbstractSpinBox.NoButtons)
        
        self.buttonRotateVertices.clicked.connect(self.rotateVertices)
        
#        self.exportButton.clicked.connect(self.exportLengths)
        

        
        self.doubleSpinBox_x.valueChanged.connect(self.plotScene)
        self.doubleSpinBox_y.valueChanged.connect(self.plotScene)
        self.doubleSpinBox_z.valueChanged.connect(self.plotScene)

        dockWidgets= self.findChildren(QDockWidget)
        for dockWidget in dockWidgets:
            action = self.menuView.addAction(dockWidget.windowTitle())
            action.setCheckable(True)
            action.setChecked(True)
            action.triggered.connect(dockWidget.setVisible)
            dockWidget.visibilityChanged.connect(action.setChecked)
#            dockWidget.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
#            dockWidget.updateGeometry()
#        self.dockLog.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)
#        self.dockLog.updateGeometry()
        
        self.buttonSamplingPhiTheta.clicked.connect(self.runSamplingPhiTheta)
        
        self._possibleParametrizedVertices = {'[2, 1, 6, 3, 7]': [2, 1, 6, 3, 7],
                                                 '[2, 3, 1, 7, 6]': [2, 3, 1, 7, 6],
                                                 '[2, 6, 1, 7, 3]': [2, 6, 1, 7, 3],
                                                 '[2, 7, 6, 3, 1]': [2, 7, 6, 3, 1],
                                                 '[3, 1, 2, 4, 7]': [3, 1, 2, 4, 7],
                                                 '[3, 2, 1, 7, 4]': [3, 2, 1, 7, 4],
                                                 '[3, 4, 1, 7, 2]': [3, 4, 1, 7, 2],
                                                 '[3, 7, 2, 4, 1]': [3, 7, 2, 4, 1],
                                                 '[4, 1, 3, 5, 7]': [4, 1, 3, 5, 7],
                                                 '[4, 3, 1, 7, 5]': [4, 3, 1, 7, 5],
                                                 '[4, 5, 1, 7, 3]': [4, 5, 1, 7, 3],
                                                 '[4, 7, 3, 5, 1]': [4, 7, 3, 5, 1],
                                                 '[5, 1, 4, 6, 7]': [5, 1, 4, 6, 7],
                                                 '[5, 4, 1, 7, 6]': [5, 4, 1, 7, 6],
                                                 '[5, 6, 1, 7, 4]': [5, 6, 1, 7, 4],
                                                 '[5, 7, 4, 6, 1]': [5, 7, 4, 6, 1],
                                                 '[6, 1, 5, 2, 7]': [6, 1, 5, 2, 7],
                                                 '[6, 2, 1, 7, 5]': [6, 2, 1, 7, 5],
                                                 '[6, 5, 1, 7, 2]': [6, 5, 1, 7, 2],
                                                 '[6, 7, 5, 2, 1]': [6, 7, 5, 2, 1]
                                                 }
        self.comboBoxParamVert.setInsertPolicy(6)
        for comb in self._possibleParametrizedVertices:
            self.comboBoxParamVert.addItem(comb)
        
        
        self.buttonFindMore.clicked.connect(self.findMoreEmbeddings)
        
        self.tabifyDockWidget(self.dockBranches, self.dockSceneShift)
        self.tabifyDockWidget(self.dockSceneShift, self.dockSequenceNavigation)
        
        

    def fullScreen(self, bool):
        if bool:
            self.showFullScreen()
        else:
            self.showMaximized()

    def loadLengths(self):
        options = QFileDialog.Options()
#        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Load lengths", "./examples","Python Objects (*.p);;All Files (*)", options=options)
        if fileName:
            try:
                self.printLog('Lengths loaded from: '+fileName)
                lengths,  R26 = pickle.load(open(fileName,'rb'))
                lengths['26'] = R26
                self.graph.setLengthsAndUpdateFixedTriangle(lengths)
#                self.graph.setR26(R26)
                blocked = self.doubleSpinBoxR26.blockSignals(True)
                self.update_graph2R26()
                self.update_graph2tabLengths()
                self.update_graph2phi()
                self.doubleSpinBoxR26.blockSignals(blocked)
                
                self.computeCouplerCurves()
            except Exception as e:
                self.showError('Some problem with loading: \n'+str(e))

    def saveLengths(self):
        options = QFileDialog.Options()
#        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Save lengths","./examples","Python Objects (*.p);;All Files (*)", options=options)
        if fileName:
            self.printLog('Lengths saved to: '+fileName)
            pickle.dump([self.graph._lengths, self.graph.getR26()], open(fileName, 'wb'))

    def insertLengths(self):
        text, ok = QInputDialog.getMultiLineText(self, 'Insert lengths', 'Insert lengths as dictionary or list of squares of lengths:')
        if ok:
            try:
                evaluated_expr = ast.literal_eval(str(text))
                if not type(evaluated_expr) is dict:
                    g = evaluated_expr
                    lengths = {
                    '12': np.sqrt(g[0]), 
                    '13': np.sqrt(g[1]), 
                    '14': np.sqrt(g[2]), 
                    '15': np.sqrt(g[3]), 
                    '16': np.sqrt(g[4]), 
                    '27': np.sqrt(g[5]), 
                    '37': np.sqrt(g[6]), 
                    '47': np.sqrt(g[7]), 
                    '57': np.sqrt(g[8]), 
                    '67': np.sqrt(g[9]), 
                    '23': np.sqrt(g[10]), 
                    '34': np.sqrt(g[11]), 
                    '45': np.sqrt(g[12]), 
                    '56': np.sqrt(g[13]),
                    '26': np.sqrt(g[14])
                    }
#                    lengths = {
#                    '12': np.sqrt(g[0]), 
#                    '13': np.sqrt(g[1]), 
#                    '14': np.sqrt(g[2]), 
#                    '15': np.sqrt(g[3]), 
#                    '16': np.sqrt(g[4]), 
#                    '27': np.sqrt(g[7]), 
#                    '37': np.sqrt(g[9]), 
#                    '47': np.sqrt(g[11]), 
#                    '57': np.sqrt(g[13]), 
#                    '67': np.sqrt(g[14]), 
#                    '23': np.sqrt(g[5]), 
#                    '34': np.sqrt(g[8]), 
#                    '45': np.sqrt(g[10]), 
#                    '56': np.sqrt(g[12]),
#                    '26': np.sqrt(g[6])
#                    }
                else:
                    lengths = evaluated_expr
                
                if type(lengths)==dict:
                    self.printLog('Inserted lengths: ')
                    self.graph.setLengthsAndUpdateFixedTriangle(lengths)
#                    self.graph.setR26(R26)
                    self.update_graph2R26()
                    self.update_graph2tabLengths()
                    self.update_graph2phi()
                    self.computeCouplerCurves()
                else:
                    self.showError('Input must be list containing dictionary of lengths and float R26 or list of squares of lengths')
            except Exception as e:
                 self.showError('Problem with input: \n'+str(e))

    def insertLengthsByEmbeddings(self):
        text, ok = QInputDialog.getMultiLineText(self, 'Insert embedding', 'Insert list of coordinates of v1,..., v7:')
        if ok:
            try:
                evaluated_expr = ast.literal_eval(str(text))
                v1, v2, v3, v4, v5, v6, v7 = evaluated_expr
                lengths = {
                '12': self.graph.dist(v1, v2), 
                '13': self.graph.dist(v1, v3), 
                '14': self.graph.dist(v1, v4), 
                '15': self.graph.dist(v1, v5), 
                '16': self.graph.dist(v1, v6), 
                '27': self.graph.dist(v2, v7), 
                '37': self.graph.dist(v3, v7), 
                '47': self.graph.dist(v4, v7), 
                '57': self.graph.dist(v5, v7), 
                '67': self.graph.dist(v6, v7), 
                '23': self.graph.dist(v2, v3), 
                '34': self.graph.dist(v3, v4), 
                '45': self.graph.dist(v4, v5), 
                '56': self.graph.dist(v5, v6),
                '26': self.graph.dist(v2, v6)
                }
            
                self.printLog('Inserted lengths: ')
                self.graph.setLengthsAndUpdateFixedTriangle(lengths)
                self.update_graph2R26()
                self.update_graph2tabLengths()
                self.update_graph2phi()
                self.computeCouplerCurves()
            except Exception as e:
                 self.showError('Problem with input: \n'+str(e))

    def showError(self, s):
        msg = QErrorMessage(self)
#        msg.setWindowModality(QtCore.Qt.WindowModal)
        msg.showMessage(s)
        self.printLog(s)

    def init_plot(self):
        """Initialize the plot of the axes"""

        self.fig.subplots_adjust(left=0.0,top=1,right=1,bottom=0.0)

        self._branches_plot = self.fig.add_subplot(111,  projection='3d')
        self._branches_plot.axis(aspect='equal')
        self._branches_plot.auto_scale_xyz([-1, 1], [-1, 1], [-1, 1])
        self._firstPlot = True

    def init_mpl(self):
        """Initialize the mpl canvas and connect the relevant parts"""

        self.canvas = FigureCanvas(self.fig)
#        self.canvas.updateGeometry()
        self.canvas.draw()
        self.figLayout.addWidget(self.canvas)

        NavigationToolbar.toolitems = [t for t in NavigationToolbar.toolitems
            if t[0] in ("Save", )]#,"Pan","Zoom","Subplots")]
        self.img_tb = NavigationToolbar(self.canvas, self, coordinates=False)
        actions = self.img_tb.findChildren(QAction)
        for a in actions:
            if a.text() == 'Customize':
                self.img_tb.removeAction(a)
                break

        self.addToolBar(self.img_tb)
        self.buttonAxelVisualisation = QAction('Axel visualisation', self)
        self.img_tb.addAction(self.buttonAxelVisualisation)
        self.buttonAxelVisualisation.triggered.connect(self.axelVisualisation) 
        action = self.menuInputOutput.addAction('Axel visualisation')
        action.triggered.connect(self.axelVisualisation)
        
        
        
#        self.centralWidget().layout().insertWidget(0,self.img_tb)

#    def init_data(self):
#        """Initialize the lengths"""
#        self.update_tabLengths2graph()
#        self.updateParameter()
#        self.update_R26toGraph()
 
    def updateParameter(self):
        N = self.graph.getSamples()
        
        self.spinBoxParameter.setMaximum(N)
        
        blocked = self.spinBoxParameter.blockSignals(True)
        self.spinBoxParameter.setValue(N/2)
        self.spinBoxParameter.blockSignals(blocked)
        
        self.spinBoxParameterStep.setValue(N/20)

    
    def update_parameterStep(self):
        self.spinBoxParameter.setSingleStep(self.spinBoxParameterStep.value())
    
    def updateDisplayedGraph(self):
        if self.isComputed():
            self.plotScene()
        else:
            self.showError('Recomputation of coupler curve needed!')

    
    def update_graph2tabLengths(self):
        self.printLog('Tab Lengths updated from graph', verbose=1)
        for e in self.doubleSpinBoxLengths:
            blocked = self.doubleSpinBoxLengths[e].blockSignals(True)
            self.doubleSpinBoxLengths[e].setValue(self.graph.getEdgeLength(e))
            self.doubleSpinBoxLengths[e].blockSignals(blocked)
            self.printLog(e+': '+str(self.doubleSpinBoxLengths[e].value()), verbose=2)
    
    def update_tabLengths2graph(self):
        self.printLog('Graph updated from tab Lengths', verbose=1)
        lengths = {
                '12': self.doubleSpinBoxL12.value(), 
                '13': self.doubleSpinBoxL13.value(), 
                '14': self.doubleSpinBoxL14.value(), 
               '15': self.doubleSpinBoxL15.value(), 
               '16': self.doubleSpinBoxL16.value(), 
               '27': self.doubleSpinBoxL27.value(), 
               '37': self.doubleSpinBoxL37.value(), 
               '47': self.doubleSpinBoxL47.value(), 
               '57': self.doubleSpinBoxL57.value(), 
               '67': self.doubleSpinBoxL67.value(), 
               '23': self.doubleSpinBoxL23.value(), 
               '34': self.doubleSpinBoxL34.value(), 
               '45': self.doubleSpinBoxL45.value(), 
               '56': self.doubleSpinBoxL56.value(), 
               '26': self.doubleSpinBoxR26.value()
                 }
        self.graph.setLengthsAndUpdateFixedTriangle(lengths)
        self.update_graph2yV2()
        self.update_graph2phi()
        self.setRequiresRecomputing()
    
    def update_phi2graph(self):
        self.graph.setPhiDegree(self.doubleSpinBoxPhi.value())
        
        blocked = self.doubleSpinBoxY2.blockSignals(True)
        self.doubleSpinBoxY2.setValue(self.graph.getyV2())
        self.doubleSpinBoxY2.blockSignals(blocked)
        
        self.update_yV2toGraph()
    
    def update_theta2graph(self):
        self.graph.setThetaDegree(self.doubleSpinBoxTheta.value())
        blocked = self.doubleSpinBoxR26.blockSignals(True)
        self.doubleSpinBoxR26.setValue(self.graph.getR26())
        self.doubleSpinBoxR26.blockSignals(blocked)
        self.update_R26toGraph()
    
    def update_R26toGraph(self):
        self.printLog('Graph updated from R26')
        self.graph.setR26(self.doubleSpinBoxR26.value())
        self.labelRecomputePHC.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')
        self.update_graph2theta()
        if self.isComputed():
            self.graph.computeIntersections()
            self.plotScene()
        else:
            self.showError('Recomputation of coupler curve needed!')

    def update_yV2toGraph(self):
        self.printLog('v2 changed')
        self.labelRecomputePHC.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')
        self.graph.setyV2(self.doubleSpinBoxY2.value())
        self.update_graph2tabLengths()
        self.update_graph2phi()
        if self.isComputed():
            self.graph.computeIntersections()
            self.plotScene()
        else:
            self.showError('Recomputation of coupler curve needed!')



    def update_graph2R26(self):
        blocked = self.doubleSpinBoxR26.blockSignals(True)
        self.doubleSpinBoxR26.setValue(self.graph.getR26())
        self.doubleSpinBoxR26.blockSignals(blocked)
        self.update_graph2yV2()
        self.update_graph2theta()

    def update_graph2yV2(self):
        blocked = self.doubleSpinBoxY2.blockSignals(True)
        self.doubleSpinBoxY2.setValue(self.graph.getyV2())
        self.doubleSpinBoxY2.blockSignals(blocked)
        
        self.update_graph2phi()
        

    def update_graph2phi(self):
        self.update_graph2theta()
        blocked = self.doubleSpinBoxPhi.blockSignals(True)
        self.doubleSpinBoxPhi.setValue(self.graph.getPhiDegree())
        self.doubleSpinBoxPhi.blockSignals(blocked)
    
    def update_graph2theta(self):
        blocked = self.doubleSpinBoxTheta.blockSignals(True)
        self.doubleSpinBoxTheta.setValue(self.graph.getThetaDegree())
        self.doubleSpinBoxTheta.blockSignals(blocked)
    
    def setRequiresRecomputing(self):
        self.graph.setRequiresRecomputing(propagate=False)
        self.labelComputed.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')
        self._V6fromPHC = []
        self.labelRecomputePHC.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')

    def setComputed(self):
        self.graph.setComputed(propagate=False)
        self.labelComputed.setText('Ready')
        c_x, c_y, c_z = self.graph.getCenterOfGravity()
        self.doubleSpinBox_x.setValue(-c_x)
        self.doubleSpinBox_y.setValue(-c_y)
        self.doubleSpinBox_z.setValue(-c_z)
    
    def isComputed(self):
        return self.graph.isComputed()

    def computeCouplerCurves(self):
        N = self.spinBoxSamples.value()
        
        self.printLog('Computing coupler curve')
        self.pushButtonPlot.setEnabled(False)
        
        self.graph.computeCouplerCurve(N)
#            self.graph.computeIntersections()
        self.updateParameter()
        self.updateDisplayedGraph()
        self.pushButtonPlot.setEnabled(True)
    
#    def tabChanged(self):
#        if self.tabWidget.currentWidget()==self.tabLengths:
#            self.printLog('Tab changed to Lengths')
#        elif self.tabWidget.currentWidget()==self.tabSequence:
#            self.plotGraphFromSequence()

    def plotScene(self):
#        self.printLog('Updating branches plot')
#        def drawLine(points,  color='black'):
#            self._branches_plot.plot([x for x, y, z in points], [y for x, y, z in points], [z for x, y, z in points],  color=color)
#        
        c_x, c_y, c_z = self.doubleSpinBox_x.value(), self.doubleSpinBox_y.value(), self.doubleSpinBox_z.value()
        def draw(points, style='black'):
            self._branches_plot.plot([x+c_x for x, y, z in points], [y+c_y for x, y, z in points], [z+c_z for x, y, z in points],  style)
        
        if self.isComputed():
            if not self._firstPlot:
                xlim = self._branches_plot.get_xlim()
                ylim = self._branches_plot.get_ylim()
                zlim = self._branches_plot.get_zlim()
                minbound = min([xlim[0], ylim[0], zlim[0]])
                maxbound = max([xlim[1], ylim[1], zlim[1]])
            else:
                minbound = self.graph.minbound
                maxbound = self.graph.maxbound
                self._firstPlot = False
            
            self._branches_plot.cla()
            
            if self.comboBoxGraphFor.currentText()=='orange':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 0, 0)
            elif self.comboBoxGraphFor.currentText()=='red':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 1, 0)
            elif self.comboBoxGraphFor.currentText()=='green':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 0, 1)
            elif self.comboBoxGraphFor.currentText()=='blue':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 1, 1)
            else:
                pos = None
            if pos:
                v1, v2, v3, v4, v5, v6, v7 = pos
                draw([v2, v3, v7, v2, v1, v3, v4, v5, v1, v4, v7, v6, v5, v7])
                draw([v1, v6])
                draw([v1, v2, v3], 'ko')
                draw([v6, v6], 'y^')
                for i, v in enumerate(pos):
                    self._branches_plot.text(v[0]+0.1+c_x, v[1]+c_y, v[2]+c_z, 'v'+str(i+1))
            
            draw(self.graph.intersections, 'ro')
            for c in self.checkBoxBranches:
                if self.checkBoxBranches[c].checkState():
                    for part in self.graph.getBranch(c):
                        draw(part, c)
            if self.checkBoxMirror.checkState():
                draw(self.graph.intersections_mirror, 'ro')
                for c in self.checkBoxBranches:
                    if self.checkBoxBranches[c].checkState():
                        for part in self.graph.getMirrorBranch(c):
                            draw(part, 'dark'+c)
        
            self._branches_plot.set_xlim3d(minbound, maxbound)
            self._branches_plot.set_ylim3d(minbound, maxbound)
            self._branches_plot.set_zlim3d(minbound, maxbound)
            self._branches_plot.get_proj()
        else:
            self.showError('Recomputation of coupler curve needed!')
        if self._V6fromPHC:
            draw(self._V6fromPHC, 'go')
#        self._branches_plot.auto_scale_xyz([minbound, maxbound], [minbound, maxbound], [minbound, maxbound])
        
#        x, y, z = self.graph.getSphere(100)
#        self._branches_plot.plot_surface(x, y, z, alpha=0.2, color='grey',  shade=True, linewidth=0)

        self._branches_plot.set_axis_off()

        self.canvas.draw()

    def printLog(self,s, verbose=0, newLine=True):
#        self._log += s +'\n'
#        self.plainTextEdit.setPlainText(self._log)
        if verbose<=self.verbose:
            if newLine:
                print s
                self.plainTextEdit.appendPlainText(str(s))
            else:
                print s,
                self.plainTextEdit.insertPlainText(str(s))
        self.plainTextEdit.moveCursor(QTextCursor.End)
        self.plainTextEdit.ensureCursorVisible()
        QApplication.processEvents()
    
    def setGraphSequence(self, seq,  seq_comments):
        self.graph_sequence = seq
        self.graph_sequence_comments = seq_comments
        
        self.spinBoxSeqLength.setValue(len(seq))
        
        self.spinBoxImgInSeq.setMinimum(1)
        self.spinBoxImgInSeq.setMaximum(len(self.graph_sequence))
        self.noteToImgInSeq.setPlainText(self.graph_sequence_comments[self.spinBoxImgInSeq.value()-1])

        self.boxInfoImgInSeq.setVisible(True)

        self.pushButtonPlot.setEnabled(False)
        N = self.spinBoxSamples.value()
        for graph in self.graph_sequence:
            graph.computeCouplerCurve(N)
        self.plotGraphFromSequence()
        self.pushButtonPlot.setEnabled(True)
#   
#[48, [
#          [
#           [['', 0], [[5, 6, 1, 7, 4], 1], [[4, 5, 1, 7, 3], 1], [[5, 6, 1, 7, 4], 1], [[3, 2, 1, 7, 4], 1]],
#        [{(1, 2): 1.99993774567597, (2, 7): 10.53609172287933, (4, 7): 10.53572330314948, (2, 6): 1.001987710974071, (6, 7): 10.53647884635266, 
#          (5, 6): 4.504386360535721, (5, 7): 11.26064825798889, (1, 4): 2.003436460984393, (1, 5): 4.45587900107854, (1, 3): 1.99476987780024, 
#          (1, 6): 2.000134247468136, (4, 5): 4.1810358189222256, (3, 7): 10.53631716364608, (3, 4): 1.0036864448806, (2, 3): 0.999614322089483}, 
#          
#          {(1, 2): 1.99993774567597, (2, 7): 10.53609172287933, (4, 7): 12.831335036787282, (2, 6): 1.001987710974071, (6, 7): 10.53647884635266, 
#          (5, 6): 4.504386360535721, (5, 7): 11.26064825798889, (1, 4): 7.588587051351411, (1, 5): 4.45587900107854, (1, 3): 1.99476987780024, 
#          (1, 6): 2.000134247468136, (4, 5): 11.306894660621534, (3, 7): 10.53631716364608, (3, 4): 7.2894600859795649, (2, 3): 0.999614322089483}, 
#          
#          {(1, 2): 1.99993774567597, (2, 7): 10.53609172287933, (4, 7): 12.831335036787282, (2, 6): 1.001987710974071, (6, 7): 10.53647884635266, 
#          (5, 6): 3.674491243171125, (5, 7): 10.99293129764172, (1, 4): 7.588587051351411, (1, 5): 3.7261764241301862, (1, 3): 1.99476987780024, 
#          (1, 6): 2.000134247468136, (4, 5): 6.5316212079357383, (3, 7): 10.53631716364608, (3, 4): 7.289460085979565, (2, 3): 0.999614322089483}, 
#          
#          {(1, 2): 1.99993774567597, (2, 7): 10.53609172287933, (4, 7): 12.831335036787282, (2, 6): 1.001987710974071, (6, 7): 10.53647884635266, 
#          (5, 6): 3.674491243171125, (5, 7): 10.99293129764172, (1, 4): 7.588587051351411, (1, 5): 3.7261764241301862, (1, 3): 1.9347054477961303, 
#          (1, 6): 2.000134247468136, (4, 5): 6.531621207935738, (3, 7): 10.524592180929439, (3, 4): 7.3736846540081995, (2, 3): 0.5698832598284466}]]
#      ]]   
   
    def loadSequence(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"Load sequence", "./sequences","Python Objects (*.p);;All Files (*)", options=options)
        if fileName:
            try:
                self.printLog('Sequence loaded from: '+fileName)
                seq = pickle.load(open(fileName,'rb'))
                graph_sequence = []
                graph_sequence_comments = []
                if type(seq[0])==list:
                    for g in seq:
                        lengths = {
                                '12': np.sqrt(g[0]), 
                                '13': np.sqrt(g[1]), 
                                '14': np.sqrt(g[2]), 
                                '15': np.sqrt(g[3]), 
                                '16': np.sqrt(g[4]), 
                                '27': np.sqrt(g[7]), 
                                '37': np.sqrt(g[9]), 
                                '47': np.sqrt(g[11]), 
                                '57': np.sqrt(g[13]), 
                                '67': np.sqrt(g[14]), 
                                '23': np.sqrt(g[5]), 
                                '34': np.sqrt(g[8]), 
                                '45': np.sqrt(g[10]), 
                                '56': np.sqrt(g[12]), 
                                '26': np.sqrt(g[6])
                                }
                    graph_sequence.append(GraphCouplerCurve(lengths=lengths, window=self))
                    try:
                        graph_sequence_comments.append(str(g[15]))
                    except:
                        graph_sequence_comments.append(str(' '))
                else:
                    self.printLog(seq)
                    for len_seq in seq[1]:
                        for lengths in len_seq[1]:
                            graph_sequence.append(GraphCouplerCurve(lengths=lengths, window=self))
                            graph_sequence_comments.append(str(' '))

                self.setGraphSequence(graph_sequence, graph_sequence_comments)
            except Exception as e:
                self.showError('Some problem with loading: \n'+str(e))

    def plotGraphFromSequence(self):
        if self.graph_sequence:
            self.setActiveGraph(self.graph_sequence[self.spinBoxImgInSeq.value()-1])
            self._V6fromPHC = []
            self.noteToImgInSeq.setPlainText(self.graph_sequence_comments[self.spinBoxImgInSeq.value()-1])
            self.update_graph2tabLengths()
            self.update_graph2R26()           
            self.updateParameter()
            self.updateDisplayedGraph()
            self.labelRecomputePHC.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')
    
    def runPHC(self):
        self.buttonRunPHC.setEnabled(False)
        self._V6fromPHC = self.graph.getSolutionsForV6(usePrev=False)
        num_sol = len(self._V6fromPHC)
        self.printLog('Number of real solutions by PHC:')
        self.printLog(str(num_sol))
        self.spinBoxNumberReal.setValue(num_sol)
        self.plotScene()
        self.labelRecomputePHC.setText('Green points OK')
        self.buttonRunPHC.setEnabled(True)
        return num_sol

    
    def rotateVertices(self):
        self.buttonRotateVertices.setEnabled(False)
        self.printLog('Rotating labeling of vertices:')
        rotated_lengths = {
                        '12' : self.graph.getEdgeLength('13'),
                        '13' : self.graph.getEdgeLength('14'),
                        '14' : self.graph.getEdgeLength('15'),
                        '15' : self.graph.getEdgeLength('16'),
                        '16' : self.graph.getEdgeLength('12'),
                        '27' : self.graph.getEdgeLength('37'),
                        '37' : self.graph.getEdgeLength('47'),
                        '47' : self.graph.getEdgeLength('57'),
                        '57' : self.graph.getEdgeLength('67'),
                        '67' : self.graph.getEdgeLength('27'),
                        '23' : self.graph.getEdgeLength('34'),
                        '34' : self.graph.getEdgeLength('45'),
                        '45' : self.graph.getEdgeLength('56'),
                        '56' : self.graph.getEdgeLength('26'), 
                        '26' : self.graph.getEdgeLength('23')
                        }
        self.graph.setLengthsAndUpdateFixedTriangle(rotated_lengths)

        blocked = self.doubleSpinBoxR26.blockSignals(True)
        self.update_graph2R26()
        self.update_graph2tabLengths()
        self.doubleSpinBoxR26.blockSignals(blocked)
        
        self.computeCouplerCurves()
        self.buttonRotateVertices.setEnabled(True)

    def showDialog(self, texts):
        dialog = QDialog(self)

        dialog.textBrowser = QPlainTextEdit(dialog)
        dialog.verticalLayout = QVBoxLayout(dialog)
        dialog.verticalLayout.addWidget(dialog.textBrowser)
        dialog.setMinimumSize(600, 300)
        
        for s in texts:
            dialog.textBrowser.appendPlainText(str(s))
        
        dialog.show()
    
    def exportLengths(self):
        len_vect = []
        for e in ['12', '13', '14', '15', '16', '23']:
            len_vect.append(self.graph.getEdgeLength(e))
        len_vect.append(self.graph.getR26())
        for e in ['27', '34', '37', '45', '47', '56', '57', '67']:
            len_vect.append(self.graph.getEdgeLength(e))
        squared=True
        if squared:
            tmp = copy.copy(len_vect)
            len_vect = []
            for r in tmp:
                len_vect.append(r**2)
#        print len_vect
        
        dialog = QDialog(self)
#        dialog.buttonBox = QDialogButtonBox(self)
#        dialog.buttonBox.setOrientation(QtCore.Qt.Horizontal)
#        dialog.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)

        dialog.textBrowser = QPlainTextEdit(dialog)
        dialog.textBrowser.appendPlainText(str(len_vect))
        len_dict = copy.deepcopy(self.graph._lengths)
        len_dict[(2, 6)] = self.graph.getR26()
        dialog.textBrowser.appendPlainText(str(len_dict))
        

        dialog.verticalLayout = QVBoxLayout(dialog)
        dialog.verticalLayout.addWidget(dialog.textBrowser)
#        dialog.verticalLayout.addWidget(dialog.buttonBox)
        dialog.setMinimumSize(600, 300)
        
        dialog.exec_()
        
#        msg = QMessageBox(self)
##        msg.setIcon(QMessageBox.Information)
#        msg.setStandardButtons(QMessageBox.Close)
##        msg.setText("Lengths:")
#        msg.setWindowTitle("Lengths export")
#        msg.setDetailedText(str(len_vect))
#        msg.show()
    def axelVisualisation(self):       
#        #colors
        blue = [0,0,255]
        cyan = [0,255,255]
        green = [0,255,0]
        magenta=[255,0,255]
        dark_yellow=[100,100,0]
        yellow=[255,255,0]
        dark_green=[0,190,0]
        red = [255,0,0]
        black = [0,0,0]
        
        colors = {
                  'orange': yellow,
                  'red': red,
                  'green': green,
                  'blue': blue, 
                  'darkorange': dark_yellow,
                  'darkred': cyan,
                  'darkgreen': dark_green,
                  'darkblue': magenta
                  }    
        
        if self.isComputed():
            v = VisualizationAxel("/home/jan/Programs/miniconda3/bin/axel-s")
            if self.comboBoxGraphFor.currentText()=='orange':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 0, 0)
            elif self.comboBoxGraphFor.currentText()=='red':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 1, 0)
            elif self.comboBoxGraphFor.currentText()=='green':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 0, 1)
            elif self.comboBoxGraphFor.currentText()=='blue':
                pos = self.graph.getPositionsOfVertices(self.spinBoxParameter.value()-1, 0, 1, 1)
            else:
                pos = None
            if pos:
                v1, v2, v3, v4, v5, v6, v7 = pos
                v.add_polyline([v2, v3, v7, v2, v1, v3, v4, v5, v1, v4, v7, v6, v5, v7], color=black)
                v.add_polyline([v1, v6], color=black)
                v.add_points([v1, v2, v3], color=black,  size=1.5)
#                for i, v in enumerate(pos):
#                    self._branches_plot.text(v[0]+0.1, v[1], v[2], 'v'+str(i+1))
            
            v.add_points(self.graph.intersections, color=red, size=1.0)
            for c in self.checkBoxBranches:
                if self.checkBoxBranches[c].checkState():
                    for part in self.graph.getBranch(c):
                        v.add_polyline(part, color=colors[c])
            
            if self.checkBoxMirror.checkState():
                v.add_points(self.graph.intersections_mirror, color=red, size=1.0)
                for c in self.checkBoxBranches:
                    if self.checkBoxBranches[c].checkState():
                        for part in self.graph.getMirrorBranch(c):
                            v.add_polyline(part, color=colors['dark'+c])
            
            if self._V6fromPHC:
                v.add_points(self._V6fromPHC, color=dark_green, size=1.0)
            v.show()
        else:
            self.showError('Recomputation of coupler curve needed!')
    
    def runSamplingPhiTheta(self):
        self.computeCouplerCurves()
        n = 0
        first = True
        while n<48 and (not self.interrupt.checkState() or first):
            first = False
            self.printLog('Sampling phi and theta')
            alg = AlgRealEmbeddings('Max7vertices', window=self)
            alg.runSamplingPhiTheta(self.graph.getLengths(),
                                    self.spinBoxSamplesPhi.value(), self.spinBoxSamplesTheta.value(), 
                                    self._possibleParametrizedVertices[self.comboBoxParamVert.currentText()])
            self.printLog('Sampling finished, see sequence')
            
            if not self.interrupt.checkState():
                self.printLog('Rotating:')
                self.rotateVertices()
    
    def findMoreEmbeddings(self):
        self.computeCouplerCurves()
        self.printLog('Searching more embeddings:')
        alg = AlgRealEmbeddings('Max7vertices', window=self, num_phi=self.spinBoxSamplesPhi.value(), num_theta=self.spinBoxSamplesTheta.value(), name=self.lineEditName.text())
        alg.findMoreEmbeddings(self.graph.getLengths())
        
        if not self.interrupt.checkState():
            self.printLog('Rotating:')
            self.rotateVertices()
    
    def showClusters(self, clusters, centers):
        pass
        newDialog = QDialog(self)
        newDialog.figure = Figure()
        newDialog.canvas = FigureCanvas(newDialog.figure)
        
        layout = QVBoxLayout()
        layout.addWidget(newDialog.canvas)
        newDialog.setLayout(layout)
        
        ax = newDialog.figure.add_subplot(111)
        for cluster in clusters:
            ax.plot([x for x, y in cluster], [y for x, y in cluster], 'o')
#        print centers
#        print [x for x, y in centers]
        ax.plot([x for x, y in centers], [y for x, y in centers], 'ro')
        
        newDialog.canvas.draw()
        
        newDialog.show()

if __name__=="__main__":
    import sys

    app = QApplication(sys.argv)
    main = MplWindow()
    main.show()
    sys.exit(app.exec_())
