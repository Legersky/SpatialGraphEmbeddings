

from PyQt5 import QtCore
from PyQt5.QtWidgets import QAction,  QApplication, QInputDialog, QFileDialog, QErrorMessage, QAbstractSpinBox, QPlainTextEdit, QDialog, QVBoxLayout
from PyQt5.QtGui import QIcon,  QKeySequence,  QTextCursor

from PyQt5.uic import loadUiType

import pickle
import ast
#import time 

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_agg import FigureCanvasAgg


from graphEmbedding import *
from axel_vis import *

UI_MainWindow, MainWindow = loadUiType("couplerCurve.ui")
class MplWindow(UI_MainWindow, MainWindow):

    def __init__(self):
        super(MplWindow, self).__init__()
        self.setupUi(self)
        
        self.showMaximized()
        
        self.fig = Figure()
#        self.graph = GraphEmbedding(window=self)
        self.graph_sequence = []
        self.graph_sequence_num_intersections = []
        
        self._V6fromPHC = []
#        self._log = ''
#        self.setRequiresRecomputing()
        self.setActiveGraph(GraphEmbedding(window=self))
        self.init_widgets()
        self.init_mpl()
#        self.init_data()
        self.init_plot()
        #self.init_pick()
    
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
        
        for e in self.doubleSpinBoxLengths:
            self.doubleSpinBoxLengths[e].valueChanged.connect(self.update_tabLengths2graph)
        
        self.update_graph2R26()
        self.doubleSpinBoxR26.valueChanged.connect(self.update_R26toGraph)
        
        self.doubleSpinBoxY2.valueChanged.connect(self.update_yV2toGraph)
        
        
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

        self.saveLengthsButton = QAction(QIcon.fromTheme("document-save"),"Save lengths",self)
        self.doc_tb.addAction(self.saveLengthsButton)
        self.saveLengthsButton.triggered.connect(self.saveLengths)
     
        self.insertLengthsButton = QAction(QIcon.fromTheme("insert-text"),"Insert lengths",self)
        self.doc_tb.addAction(self.insertLengthsButton)
        self.insertLengthsButton.triggered.connect(self.insertLengths)
        
        self.updateParameter()
        
        self.buttonLoadSequence.clicked.connect(self.loadSequence)
        
        self.spinBoxImgInSeq.setMinimum(1)
        self.spinBoxImgInSeq.setMaximum(1)
        self.spinBoxImgInSeq.valueChanged.connect(self.plotGraphFromSequence)
        
        self.spinBoxSeqLength.setReadOnly(True)
        self.spinBoxSeqLength.setButtonSymbols(QAbstractSpinBox.NoButtons)
        
        self.noteToImgInSeq.setReadOnly(True)
        self.boxInfoImgInSeq.setVisible(False)
        
#        self.tabWidget.currentChanged.connect(self.tabChanged)
        
        self.buttonRunPHC.clicked.connect(self.runPHC)
        self.spinBoxNumberReal.setReadOnly(True)
        self.spinBoxNumberReal.setButtonSymbols(QAbstractSpinBox.NoButtons)
        
        self.buttonRotateVertices.clicked.connect(self.rotateVertices)
        
        self.exportButton.clicked.connect(self.exportForVangelis)
        
        self.buttonAxelVisualisation.clicked.connect(self.axelVisualisation)
        
        self.doubleSpinBox_x.valueChanged.connect(self.plotScene)
        self.doubleSpinBox_y.valueChanged.connect(self.plotScene)
        self.doubleSpinBox_z.valueChanged.connect(self.plotScene)
        
#        self.dockSequenceNavigation.setVisible(False)
#        self.splitDockWidget(self.dockPHC, self.dockLog, QtCore.Qt.Horizontal)

    def loadLengths(self):
        options = QFileDialog.Options()
#        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Load lengths", "./examples","Python Objects (*.p);;All Files (*)", options=options)
        if fileName:
            try:
                self.printLog('Lengths loaded from: '+fileName)
                lengths,  R26 = pickle.load(open(fileName,'rb'))
                self.graph.setLengths(lengths)
                self.graph.setR26(R26)
                blocked = self.doubleSpinBoxR26.blockSignals(True)
                self.update_graph2R26()
                self.update_graph2tabLengths()
                self.computeCouplerCurves()
                self.doubleSpinBoxR26.blockSignals(blocked)
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
        text, ok = QInputDialog.getMultiLineText(self, 'Insert lengths', 'Insert lengths as [dictionary, R26] or list of squares of lengths:')
        if ok:
            try:
                evaluated_expr = ast.literal_eval(str(text))
                if not type(evaluated_expr[0]) is dict:
                    g = evaluated_expr
                    lengths = {
                    'L12': np.sqrt(g[0]), 
                    'L13': np.sqrt(g[1]), 
                    'L14': np.sqrt(g[2]), 
                    'L15': np.sqrt(g[3]), 
                    'L16': np.sqrt(g[4]), 
                    'L27': np.sqrt(g[7]), 
                    'L37': np.sqrt(g[9]), 
                    'L47': np.sqrt(g[11]), 
                    'L57': np.sqrt(g[13]), 
                    'L67': np.sqrt(g[14]), 
                    'L23': np.sqrt(g[5]), 
                    'L34': np.sqrt(g[8]), 
                    'L45': np.sqrt(g[10]), 
                    'L56': np.sqrt(g[12])
                    }
                    R26 = float(np.sqrt(g[6]))
                else:
                    lengths, R26 = evaluated_expr
                
                if type(lengths)==dict and type(R26)==float:
                    self.printLog('Inserted lengths: ')
                    self.graph.setLengths(lengths)
                    self.graph.setR26(R26)
                    self.update_graph2R26()
                    self.update_graph2tabLengths()
                    self.computeCouplerCurves()
                else:
                    self.showError('Input must be list containing dictionary of lengths and float R26 or list of squares of lengths')
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
        self.toolbar = NavigationToolbar(self.canvas, self, coordinates=False)
        actions = self.toolbar.findChildren(QAction)
        for a in actions:
            if a.text() == 'Customize':
                self.toolbar.removeAction(a)
                break
        self.addToolBar(self.toolbar)
#        self.centralWidget().layout().insertWidget(0,self.toolbar)

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
        self.printLog('Tab Lengths updated from graph')
        for e in self.doubleSpinBoxLengths:
            blocked = self.doubleSpinBoxLengths[e].blockSignals(True)
            self.doubleSpinBoxLengths[e].setValue(self.graph.getEdgeLength(e))
            self.doubleSpinBoxLengths[e].blockSignals(blocked)
            self.printLog(e+': '+str(self.doubleSpinBoxLengths[e].value()))
    
    def update_tabLengths2graph(self):
        self.printLog('Graph updated from tab Lengths')
        lengths = {
                'L12': self.doubleSpinBoxL12.value(), 
                'L13': self.doubleSpinBoxL13.value(), 
                'L14': self.doubleSpinBoxL14.value(), 
               'L15': self.doubleSpinBoxL15.value(), 
               'L16': self.doubleSpinBoxL16.value(), 
               'L27': self.doubleSpinBoxL27.value(), 
               'L37': self.doubleSpinBoxL37.value(), 
               'L47': self.doubleSpinBoxL47.value(), 
               'L57': self.doubleSpinBoxL57.value(), 
               'L67': self.doubleSpinBoxL67.value(), 
               'L23': self.doubleSpinBoxL23.value(), 
               'L34': self.doubleSpinBoxL34.value(), 
               'L45': self.doubleSpinBoxL45.value(), 
               'L56': self.doubleSpinBoxL56.value()
                 }
        self.graph.setLengths(lengths)
        self.update_graph2yV2()
        self.setRequiresRecomputing()
        
    def update_R26toGraph(self):
        self.printLog('Graph updated from R26')
        self.graph.setR26(self.doubleSpinBoxR26.value())
        self.labelRecomputePHC.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')
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

    def update_graph2yV2(self):
        blocked = self.doubleSpinBoxY2.blockSignals(True)
        self.doubleSpinBoxY2.setValue(self.graph.getyV2())
        self.doubleSpinBoxY2.blockSignals(blocked)
    
    def update_yV2toGraph(self):
        self.printLog('v2 changed manually')
        self.labelRecomputePHC.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')
        self.graph.setyV2(self.doubleSpinBoxY2.value())
        self.update_graph2tabLengths()
        if self.isComputed():
            self.graph.computeIntersections()
            self.plotScene()
        else:
            self.showError('Recomputation of coupler curve needed!')
    
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
        self.canvas.draw()

    def printLog(self,s):
#        self._log += s +'\n'
#        self.plainTextEdit.setPlainText(self._log)
        self.plainTextEdit.appendPlainText(s)
        self.plainTextEdit.moveCursor(QTextCursor.End)
        self.plainTextEdit.ensureCursorVisible()
        QApplication.processEvents()
    
    def loadSequence(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"Load sequence", "./sequences","Python Objects (*.p);;All Files (*)", options=options)
        if fileName:
            try:
                self.printLog('Sequence loaded from: '+fileName)
                seq = pickle.load(open(fileName,'rb'))
                self.graph_sequence = []
                self.graph_sequence_num_intersections = []
                for g in seq:
                    lengths = {
                            'L12': np.sqrt(g[0]), 
                            'L13': np.sqrt(g[1]), 
                            'L14': np.sqrt(g[2]), 
                            'L15': np.sqrt(g[3]), 
                            'L16': np.sqrt(g[4]), 
                            'L27': np.sqrt(g[7]), 
                            'L37': np.sqrt(g[9]), 
                            'L47': np.sqrt(g[11]), 
                            'L57': np.sqrt(g[13]), 
                            'L67': np.sqrt(g[14]), 
                            'L23': np.sqrt(g[5]), 
                            'L34': np.sqrt(g[8]), 
                            'L45': np.sqrt(g[10]), 
                            'L56': np.sqrt(g[12])
                            }
                    self.graph_sequence.append(GraphEmbedding(lengths=lengths,  r26=np.sqrt(g[6]), window=self))
                    try:
                        self.graph_sequence_num_intersections.append(str(g[15]))
                    except:
                        self.graph_sequence_num_intersections.append(str(' '))

                self.spinBoxImgInSeq.setMinimum(1)
                self.spinBoxImgInSeq.setMaximum(len(self.graph_sequence))
                self.noteToImgInSeq.setPlainText(self.graph_sequence_num_intersections[self.spinBoxImgInSeq.value()-1])

                self.boxInfoImgInSeq.setVisible(True)

                self.pushButtonPlot.setEnabled(False)
                N = self.spinBoxSamples.value()
                for graph in self.graph_sequence:
                    graph.computeCouplerCurve(N)
                self.plotGraphFromSequence()
                self.pushButtonPlot.setEnabled(True)
            except Exception as e:
                self.showError('Some problem with loading: \n'+str(e))

    def plotGraphFromSequence(self):
        if self.graph_sequence:
            self.setActiveGraph(self.graph_sequence[self.spinBoxImgInSeq.value()-1])
            self._V6fromPHC = []
            self.noteToImgInSeq.setPlainText(self.graph_sequence_num_intersections[self.spinBoxImgInSeq.value()-1])
            self.update_graph2R26()
            self.update_graph2tabLengths()
            self.updateParameter()
            self.updateDisplayedGraph()
            self.labelRecomputePHC.setText('<html><head/><body><p><span style=" color:#ff0000;">Recomputation needed</span></p></body></html>')
    
    def runPHC(self):
        self.buttonRunPHC.setEnabled(False)
        self._V6fromPHC = self.graph.getSolutionsForV6()
        num_sol = len(self._V6fromPHC)
        self.printLog('Number of real solutions by PHC:')
        self.printLog(str(num_sol))
        self.spinBoxNumberReal.setValue(num_sol)
        self.plotScene()
        self.labelRecomputePHC.setText('Green points OK')
        self.buttonRunPHC.setEnabled(True)

    
    def rotateVertices(self):
        self.buttonRotateVertices.setEnabled(False)
        self.printLog('Rotating labeling of vertices:')
        rotated_lengths = {
                        'L12' : self.graph.getEdgeLength('13'),
                        'L13' : self.graph.getEdgeLength('14'),
                        'L14' : self.graph.getEdgeLength('15'),
                        'L15' : self.graph.getEdgeLength('16'),
                        'L16' : self.graph.getEdgeLength('12'),
                        'L27' : self.graph.getEdgeLength('37'),
                        'L37' : self.graph.getEdgeLength('47'),
                        'L47' : self.graph.getEdgeLength('57'),
                        'L57' : self.graph.getEdgeLength('67'),
                        'L67' : self.graph.getEdgeLength('27'),
                        'L23' : self.graph.getEdgeLength('34'),
                        'L34' : self.graph.getEdgeLength('45'),
                        'L45' : self.graph.getEdgeLength('56'),
                        'L56' : self.graph.getR26()
                        }
        R26 = self.graph.getEdgeLength('23')
        self.graph.setLengths(rotated_lengths)
        self.graph.setR26(R26)
        blocked = self.doubleSpinBoxR26.blockSignals(True)
        self.update_graph2R26()
        self.update_graph2tabLengths()
        self.doubleSpinBoxR26.blockSignals(blocked)
        
        self.computeCouplerCurves()
        self.buttonRotateVertices.setEnabled(True)
    
    def exportForVangelis(self):
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
        len_dict['L26'] = self.graph.getR26()
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

       

if __name__=="__main__":
    import sys

    app = QApplication(sys.argv)
    main = MplWindow()
    main.show()
    sys.exit(app.exec_())
