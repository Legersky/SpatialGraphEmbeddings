from graphEmbeddings3D.couplerCurve import MplWindow
from PyQt5.QtWidgets import QApplication


if __name__=="__main__":
    import sys

    app = QApplication(sys.argv)
    main = MplWindow()
    main.show()
    sys.exit(app.exec_())
