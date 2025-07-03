# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'simulator_GUI_new.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QDialog, QFileDialog
import numpy as np
import yaml

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(888, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.formLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.formLayoutWidget.setGeometry(QtCore.QRect(90, 80, 261, 121))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.formLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setObjectName("formLayout")
        self.label_num_emitter = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_num_emitter.setObjectName("label_num_emitter")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_num_emitter)
        self.number_of_emitter = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.number_of_emitter.setObjectName("number_of_emitter")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.number_of_emitter)
        self.label_image_size = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_image_size.setObjectName("label_image_size")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_image_size)
        self.Image_size = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.Image_size.setObjectName("Image_size")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.Image_size)
        self.label_num_tiles = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_num_tiles.setObjectName("label_num_tiles")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_num_tiles)
        self.num_tiles = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.num_tiles.setObjectName("num_tiles")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.num_tiles)
        self.label_14 = QtWidgets.QLabel(self.formLayoutWidget)
        self.label_14.setObjectName("label_14")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_14)
        self.photon_count = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.photon_count.setReadOnly(True)
        self.photon_count.setObjectName("photon_count")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.photon_count)
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(90, 340, 160, 141))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.wavelengh_473 = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.wavelengh_473.setObjectName("wavelengh_473")
        self.verticalLayout.addWidget(self.wavelengh_473)
        self.wavelength_561 = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.wavelength_561.setObjectName("wavelength_561")
        self.verticalLayout.addWidget(self.wavelength_561)
        self.wavelength_647 = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.wavelength_647.setObjectName("wavelength_647")
        self.verticalLayout.addWidget(self.wavelength_647)
        self.wavelength_750 = QtWidgets.QCheckBox(self.verticalLayoutWidget)
        self.wavelength_750.setObjectName("wavelength_750")
        self.verticalLayout.addWidget(self.wavelength_750)
        self.formLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.formLayoutWidget_2.setGeometry(QtCore.QRect(90, 210, 261, 111))
        self.formLayoutWidget_2.setObjectName("formLayoutWidget_2")
        self.formLayout_2 = QtWidgets.QFormLayout(self.formLayoutWidget_2)
        self.formLayout_2.setContentsMargins(0, 0, 0, 0)
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_is_cell = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.label_is_cell.setObjectName("label_is_cell")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_is_cell)
        self.label_is_nucleus = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.label_is_nucleus.setObjectName("label_is_nucleus")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_is_nucleus)
        self.is_nucleus = QtWidgets.QComboBox(self.formLayoutWidget_2)
        self.is_nucleus.setObjectName("is_nucleus")
        self.is_nucleus.addItem("")
        self.is_nucleus.addItem("")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.is_nucleus)
        self.is_cell = QtWidgets.QComboBox(self.formLayoutWidget_2)
        self.is_cell.setObjectName("is_cell")
        self.is_cell.addItem("")
        self.is_cell.addItem("")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.is_cell)
        self.label_15 = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.label_15.setObjectName("label_15")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_15)
        self.comboBox = QtWidgets.QComboBox(self.formLayoutWidget_2)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.comboBox)
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(310, 500, 113, 32))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.on_click) #!!!
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(390, 80, 311, 116))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label_3 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 2, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.cell_axis_c_max = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.cell_axis_c_max.setObjectName("cell_axis_c_max")
        self.gridLayout.addWidget(self.cell_axis_c_max, 3, 4, 1, 1)
        self.cell_axes_b_min = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.cell_axes_b_min.setObjectName("cell_axes_b_min")
        self.gridLayout.addWidget(self.cell_axes_b_min, 2, 2, 1, 1)
        self.cell_axis_a_max = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.cell_axis_a_max.setObjectName("cell_axis_a_max")
        self.gridLayout.addWidget(self.cell_axis_a_max, 1, 4, 1, 1)
        self.cell_axis_b_max = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.cell_axis_b_max.setObjectName("cell_axis_b_max")
        self.gridLayout.addWidget(self.cell_axis_b_max, 2, 4, 1, 1)
        self.cell_axes_c_min = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.cell_axes_c_min.setObjectName("cell_axes_c_min")
        self.gridLayout.addWidget(self.cell_axes_c_min, 3, 2, 1, 1)
        self.cell_axis_a_min = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.cell_axis_a_min.setObjectName("cell_axis_a_min")
        self.gridLayout.addWidget(self.cell_axis_a_min, 1, 2, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 1, 1, 1, 1)
        self.num_cells = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.num_cells.setObjectName("num_cells")
        self.gridLayout.addWidget(self.num_cells, 0, 2, 1, 1)
        self.label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 3, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 1, 3, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 2, 1, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_8.setObjectName("label_8")
        self.gridLayout.addWidget(self.label_8, 3, 1, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_9.setObjectName("label_9")
        self.gridLayout.addWidget(self.label_9, 2, 3, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_10.setObjectName("label_10")
        self.gridLayout.addWidget(self.label_10, 3, 3, 1, 1)
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(390, 210, 311, 102))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.browse_psf = QtWidgets.QPushButton(self.gridLayoutWidget_2)
        self.browse_psf.setObjectName("browse_psf")
        self.browse_psf.clicked.connect(self.browsefiles2)
        self.gridLayout_2.addWidget(self.browse_psf, 0, 2, 1, 1)
        self.psf_filename = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.psf_filename.setObjectName("psf_filename")
        self.gridLayout_2.addWidget(self.psf_filename, 0, 1, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_12.setObjectName("label_12")
        self.gridLayout_2.addWidget(self.label_12, 0, 0, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_11.setObjectName("label_11")
        self.gridLayout_2.addWidget(self.label_11, 1, 0, 1, 1)
        self.codebook_filename = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.codebook_filename.setObjectName("codebook_filename")
        self.gridLayout_2.addWidget(self.codebook_filename, 1, 1, 1, 1)
        self.browse_codebook = QtWidgets.QPushButton(self.gridLayoutWidget_2)
        self.browse_codebook.setObjectName("browse_codebook")
        self.browse_codebook.clicked.connect(self.browsefiles1)
        self.gridLayout_2.addWidget(self.browse_codebook, 1, 2, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(600, 310, 231, 211))
        self.label_13.setText("")
        self.label_13.setPixmap(QtGui.QPixmap("../resources/images/anglerfish copy.png"))
        self.label_13.setObjectName("label_13")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 888, 24))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.number_of_emitter.setPlaceholderText("300")
        self.num_tiles.setPlaceholderText("2")
        self.Image_size.setPlaceholderText("200")
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
    
    def browsefiles1(self):
        fname = QFileDialog.getOpenFileName(None, 'Open file', 'C:\Desktop')
        self.codebook_filename.setText(fname[0])

    def browsefiles2(self):
        fname = QFileDialog.getOpenFileName(None, 'Open file', 'C:\Desktop')
        self.psf_filename.setText(fname[0])

    def on_click(self):
        filepath = "config.yml"
        with open(filepath, "r") as f:
            config = yaml.safe_load(f)
        print(config)

        number_emitter = self.number_of_emitter.text() #str
        if (number_emitter != ""):
            config["simulation"]["emitter_count"] = int(number_emitter)

        number_cell = self.num_cells.text()
        config["simulation"]["cells"]["count"] = number_cell

        
        image_size = self.Image_size.text() #str
        if(image_size != ""):
            config["simulation"]["image_size"] = int(image_size)
    
        number_tiles = self.num_tiles.text() #str
        if(number_tiles != ""):
            config["simulation"]["tile_count"] = int(float(number_tiles))
        
        psf_file = self.psf_filename.text()
        if(psf_file != ""):
            config["data"]["psf_file"] = psf_file
        
        codebook_file = self.codebook_filename.text()
        if(codebook_file != ""):
            config["data"]["codebook_file"] = codebook_file
        
        is_cell = self.is_cell.currentText() #str
        if (is_cell == "True"):
            is_cell = True
        else:
            is_cell = False
        config["simulation"]["is_cell"] = is_cell
        
        is_nucleus = self.is_nucleus.currentText()
        if (is_nucleus == "True"):
            is_nucleus = True
        else:
            is_nucleus = False
        config["simulation"]["is_nucleus"] = is_nucleus

        wavelength_seq = []
    
        wavelength_473 = bool(self.wavelengh_473.checkState())
        if wavelength_473:
            wavelength_seq.append(473)
        wavelength_561 = bool(self.wavelength_561.checkState())
        if wavelength_561:
            wavelength_seq.append(561)
        wavelength_647 = bool(self.wavelength_647.checkState())
        if wavelength_647:
            wavelength_seq.append(647)
        wavelength_750 = bool(self.wavelength_750.checkState())
        if wavelength_750:
            wavelength_seq.append(750)

        dict_axes = config["simulation"]["cells"]["axes"]
        dict_axes["a"] = [int(self.cell_axis_a_min.text()), int(self.cell_axis_a_max.text())]
        dict_axes["b"] = [int(self.cell_axes_b_min.text()), int(self.cell_axis_b_max.text())]
        dict_axes["c"] = [int(self.cell_axes_c_min.text()), int(self.cell_axis_c_max.text())]

        if (len(wavelength_seq) != 0):
            config["simulation"]["signal_wavelengths"] = wavelength_seq

        print(config)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_num_emitter.setText(_translate("MainWindow", "Number of emitters:"))
        self.label_image_size.setText(_translate("MainWindow", "Image size:               "))
        self.label_num_tiles.setText(_translate("MainWindow", "Number of Tiles:      "))
        self.label_14.setText(_translate("MainWindow", "photon count:          "))
        self.photon_count.setText(_translate("MainWindow", "50000000"))
        self.wavelengh_473.setText(_translate("MainWindow", "473 nm"))
        self.wavelength_561.setText(_translate("MainWindow", "561 nm"))
        self.wavelength_647.setText(_translate("MainWindow", "647 nm"))
        self.wavelength_750.setText(_translate("MainWindow", "750 nm"))
        self.label_is_cell.setText(_translate("MainWindow", "Is_cell ?                                  "))
        self.label_is_nucleus.setText(_translate("MainWindow", "Is_nucleus?                            "))
        self.is_nucleus.setItemText(0, _translate("MainWindow", "True"))
        self.is_nucleus.setItemText(1, _translate("MainWindow", "False"))
        self.is_cell.setItemText(0, _translate("MainWindow", "True"))
        self.is_cell.setItemText(1, _translate("MainWindow", "False"))
        self.label_15.setText(_translate("MainWindow", "Objective magnification:       "))
        self.comboBox.setItemText(0, _translate("MainWindow", "60X"))
        self.pushButton.setText(_translate("MainWindow", "Done"))
        self.label_3.setText(_translate("MainWindow", "cell axis b:"))
        self.label_2.setText(_translate("MainWindow", "cell axis a:"))
        self.cell_axis_c_max.setText(_translate("MainWindow", "60"))
        self.cell_axes_b_min.setText(_translate("MainWindow", "60"))
        self.cell_axis_a_max.setText(_translate("MainWindow", "40"))
        self.cell_axis_b_max.setText(_translate("MainWindow", "70"))
        self.cell_axes_c_min.setText(_translate("MainWindow", "50"))
        self.cell_axis_a_min.setText(_translate("MainWindow", "30"))
        self.label_5.setText(_translate("MainWindow", "min"))
        self.num_cells.setText(_translate("MainWindow", "3"))
        self.label.setText(_translate("MainWindow", "Number of cells"))
        self.label_4.setText(_translate("MainWindow", "cell axis c:"))
        self.label_6.setText(_translate("MainWindow", "max"))
        self.label_7.setText(_translate("MainWindow", "min"))
        self.label_8.setText(_translate("MainWindow", "min"))
        self.label_9.setText(_translate("MainWindow", "max"))
        self.label_10.setText(_translate("MainWindow", "max"))
        self.browse_psf.setText(_translate("MainWindow", "Browse"))
        self.label_12.setText(_translate("MainWindow", "psf file"))
        self.label_11.setText(_translate("MainWindow", "codebook file"))
        self.browse_codebook.setText(_translate("MainWindow", "Browse"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

