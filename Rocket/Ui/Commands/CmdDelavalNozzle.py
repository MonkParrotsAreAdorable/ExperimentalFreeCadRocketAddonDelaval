import FreeCAD
import FreeCADGui
from PySide import QtGui, QtCore
import numpy as np
from scipy.optimize import fsolve
import Part
# we need to review the formulas form the wiki
def mach_from_pressure_ratio(Pe_Po, gamma):
    def func(M):
        return (1 + (gamma - 1)/2 * M**2) ** (-gamma/(gamma - 1)) - Pe_Po
    return fsolve(func, 2.0)[0]
def area_ratio_from_mach(M, gamma):
    term = (2/(gamma + 1)) * (1 + (gamma - 1)/2 * M**2)
    exponent = (gamma + 1) / (2*(gamma - 1))
    return (1/M) * term**exponent
# Freecad uses millimeters by default . dont forget that again .
def throat_diameter(De, Pe, Po, gamma):
    if Pe <= 0 or Po <= 0:
        raise ValueError("Pressures must be positive")
    if Pe >= Po:
        raise ValueError("Exit pressure must be less than stagnation pressure")
    Pe_Po = Pe / Po
    Me = mach_from_pressure_ratio(Pe_Po, gamma)
    area_ratio = area_ratio_from_mach(Me, gamma)   
    Ae = np.pi * (De/2)**2
    At = Ae / area_ratio
    Dt = 2 * np.sqrt(At / np.pi)
    return Dt, Me, area_ratio
# using a bslpine is better than the previous method . 
def generate_nozzle_profile(x_inlet, x_throat, x_exit,r_inlet, r_throat, r_exit,n_points=200):
    x_conv = np.linspace(x_inlet, x_throat, n_points//2 + 1)
    t_conv = (x_conv - x_inlet) / (x_throat - x_inlet)
    r_conv = r_throat + (r_inlet - r_throat) * 0.5 * (1 + np.cos(np.pi * t_conv))
    x_div = np.linspace(x_throat, x_exit, n_points//2 + 1)
    xi = (x_div - x_throat) / (x_exit - x_throat)
    r_div = r_throat + (r_exit - r_throat) * (3 - 2*xi) * xi**2
    x = np.concatenate([x_conv[:-1], x_div])
    r = np.concatenate([r_conv[:-1], r_div])
    return x, r
class FeatureDeLavalNozzle:
    def __init__(self, obj):
        obj.Proxy = self
        self.addProperties(obj)
    def addProperties(self, obj):
        # add some verification logic later , we cant have the nozzle chocke due to pressure diffrences , ask prof. for further mathematical proods or help.
        obj.addProperty("App::PropertyFloat", "Gamma","Isentropic", "Specific heat ratio (γ)").Gamma = 1.4 #need to see docs.
        obj.addProperty("App::PropertyFloat", "StagnationPressure","Isentropic", "Stagnation pressure Po (Pa)").StagnationPressure = 101325.0
        obj.addProperty("App::PropertyFloat", "ExitPressure","Isentropic", "Exit static pressure Pe (Pa)").ExitPressure = 5000.0
        obj.addProperty("App::PropertyLength", "ExitDiameter_mm","Geometry", "Exit diameter De (mm)").ExitDiameter_mm = 150.0
        obj.addProperty("App::PropertyFloat", "InletRatio","Geometry", "Inlet radius / throat radius").InletRatio = 2.5
        obj.addProperty("App::PropertyFloat", "ConvergentLengthRatio","Geometry", "Convergent length / throat diameter").ConvergentLengthRatio = 1.5
        obj.addProperty("App::PropertyFloat", "DivergentLengthRatio","Geometry", "Divergent length / (exit radius - throat radius)").DivergentLengthRatio = 8.0
        obj.addProperty("App::PropertyInteger", "ProfilePoints","Geometry", "Number of profile points").ProfilePoints = 200
        obj.addProperty("App::PropertyLength", "ThroatDiameter_mm","Results", "Computed throat diameter (mm)").ThroatDiameter_mm = 0.0
        obj.addProperty("App::PropertyFloat", "ExitMach","Results", "Exit Mach number").ExitMach = 0.0
    def execute(self, obj):
        #also, it is better if we like , generate it as a solid body , not a hollow one . Users should have the liberty to create their own housing instead of some thin one no ? 
        try:
            gamma = obj.Gamma
            Po = obj.StagnationPressure
            Pe = obj.ExitPressure
            De_mm = obj.ExitDiameter_mm.Value   
            De_m = De_mm / 1000.0    # HIGHLY LIKERY THIS IS NOT ACCURATE , NEED TO DOUBLE CHECK           
            inlet_ratio = obj.InletRatio
            Lc_ratio = obj.ConvergentLengthRatio
            Ld_ratio = obj.DivergentLengthRatio
            n_points = obj.ProfilePoints
            Dt_m, Me, _ = throat_diameter(De_m, Pe, Po, gamma)
            obj.ExitMach = Me
            obj.ThroatDiameter_mm = Dt_m * 1000.0
            r_exit_m = De_m / 2.0
            r_throat_m = Dt_m / 2.0
            r_inlet_m = inlet_ratio * r_throat_m
            x_inlet_m = 0.0
            x_throat_m = Lc_ratio * Dt_m
            if r_exit_m > r_throat_m:
                L_div_m = Ld_ratio * (r_exit_m - r_throat_m)
            else:
                L_div_m = Dt_m
            x_exit_m = x_throat_m + L_div_m
            x_m, r_m = generate_nozzle_profile(x_inlet_m, x_throat_m, x_exit_m,r_inlet_m, r_throat_m, r_exit_m,n_points)
            scale = 1000.0
            x = x_m * scale
            r = r_m * scale
            r_exit = r_exit_m * scale
            r_throat = r_throat_m * scale
            r_inlet = r_inlet_m * scale
            x_inlet = x_inlet_m * scale
            x_throat = x_throat_m * scale
            x_exit = x_exit_m * scale
            outer_pts = [FreeCAD.Vector(xi, ri, 0) for xi, ri in zip(x, r)]
            bspline = Part.BSplineCurve()
            bspline.interpolate(outer_pts)
            outer_edge = bspline.toShape()
            exit_pt = FreeCAD.Vector(x_exit, r_exit, 0)
            axis_exit = FreeCAD.Vector(x_exit, 0, 0)
            bottom_edge1 = Part.makeLine(exit_pt, axis_exit)
            axis_inlet = FreeCAD.Vector(x_inlet, 0, 0)
            bottom_edge2 = Part.makeLine(axis_exit, axis_inlet)
            inlet_pt = FreeCAD.Vector(x_inlet, r_inlet, 0)
            bottom_edge3 = Part.makeLine(axis_inlet, inlet_pt)
            wire = Part.Wire([outer_edge, bottom_edge1, bottom_edge2, bottom_edge3])
            if not wire.isValid():
                raise Exception("Wire is not valid")
            face = Part.Face(wire)
            if not face.isValid():
                raise Exception("Face from wire is invalid")
            solid = face.revolve(FreeCAD.Vector(0,0,0), FreeCAD.Vector(1,0,0), 360)
            if not solid.isValid():
                raise Exception("Revolved solid is invalid")
            obj.Shape = solid
        except Exception as e:
            FreeCAD.Console.PrintError(f"Nozzle recompute failed: {e}\n")
            obj.Shape = None
class ViewProviderDeLavalNozzle:
    def __init__(self, vobj):
        vobj.Proxy = self
    def attach(self, vobj):
        self.Object = vobj.Object
    def updateData(self, obj, prop):
        return
    def onChanged(self, vobj, prop):
        return
    def doubleClicked(self, vobj):
        doc = FreeCADGui.getDocument(vobj.Object.Document)
        if not doc.getInEdit():
            doc.setEdit(vobj.Object.Name)
        return True
    def getIcon(self):
        return ""
    def __getstate__(self):
        return None
    def __setstate__(self, state):
        return None
class DeLavalNozzleTaskPanel:
    def __init__(self):
        self.form = QtGui.QWidget()
        self.form.setWindowTitle("De Laval Nozzle Parameters")
        layout = QtGui.QVBoxLayout(self.form)
        self.gamma = self._add_field(layout, "Specific heat ratio γ:", 1.4, decimals=3)
        self.Po = self._add_field(layout, "Stagnation pressure Po (Pa):", 101325.0, decimals=0)
        self.Pe = self._add_field(layout, "Exit static pressure Pe (Pa):", 5000.0, decimals=0)
        self.De_mm = self._add_field(layout, "Exit diameter De (mm):", 150.0, decimals=1)
        self.inlet_ratio = self._add_field(layout, "Inlet radius / throat radius:", 2.5, decimals=2)
        self.Lc_ratio = self._add_field(layout, "Convergent length / throat diameter:", 1.5, decimals=2)
        self.Ld_ratio = self._add_field(layout, "Divergent length / (exit radius - throat radius):", 8.0, decimals=2)
        self.n_points = self._add_int_field(layout, "Number of profile points:", 200)
        btn_box = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel)
        btn_box.accepted.connect(self.accept)
        btn_box.rejected.connect(self.reject)
        layout.addWidget(btn_box)
    def _add_field(self, layout, label, default, decimals=1):
        hbox = QtGui.QHBoxLayout()
        lbl = QtGui.QLabel(label)
        spin = QtGui.QDoubleSpinBox()
        spin.setRange(0.001, 1e9)
        spin.setValue(default)
        spin.setDecimals(decimals)
        hbox.addWidget(lbl)
        hbox.addWidget(spin)
        layout.addLayout(hbox)
        return spin
#    def _add_int_field(self, layout, label, default):
#        hbox = QtGui.QHBoxLayout()
#        lbl = QtGui.QLabel(label)
#        spin = QtGui.QSpinBox()
#        spin.setRange(10, 1000)
#        spin.setValue(default)
#        hbox.addWidget(lbl)
#        hbox.addWidget(spin)
#        layout.addLayout(hbox)
#        return spin
    def _add_int_field(self, layout, label, default):
        hbox = QtGui.QHBoxLayout()
        lbl = QtGui.QLabel(label)
        spin = QtGui.QSpinBox()
        spin.setRange(10, 1000)
        spin.setValue(default)
        hbox.addWidget(lbl)
        hbox.addWidget(spin)
        layout.addLayout(hbox)
        return spin
    def accept(self):
        gamma = self.gamma.value()
        Po = self.Po.value()
        Pe = self.Pe.value()
        De_mm = self.De_mm.value()          # mm
        inlet_ratio = self.inlet_ratio.value()
        Lc_ratio = self.Lc_ratio.value()
        Ld_ratio = self.Ld_ratio.value()
        n_points = self.n_points.value()
        doc = FreeCAD.ActiveDocument
        if not doc:
            doc = FreeCAD.newDocument("NozzleDesign")
        obj = doc.addObject("Part::FeaturePython", "DeLavalNozzle")
        FeatureDeLavalNozzle(obj)
        if FreeCAD.GuiUp:
            ViewProviderDeLavalNozzle(obj.ViewObject)
        obj.Gamma = gamma
        obj.StagnationPressure = Po
        obj.ExitPressure = Pe
        obj.ExitDiameter_mm = De_mm
        obj.InletRatio = inlet_ratio
        obj.ConvergentLengthRatio = Lc_ratio
        obj.DivergentLengthRatio = Ld_ratio
        obj.ProfilePoints = n_points
        doc.recompute()
        FreeCADGui.Control.closeDialog()
        FreeCADGui.Selection.clearSelection()
        FreeCADGui.Selection.addSelection(obj)
    def reject(self):
        FreeCADGui.Control.closeDialog()
class _DeLavalNozzle:
    def GetResources(self):
        return {
            "Pixmap": FreeCAD.getUserAppDataDir() + "Mod/Rocket/Resources/icons/DeLavalNozzle.svg",  
            "MenuText": "De Laval Nozzle",
            "ToolTip": "Create a de Laval nozzle using isentropic flow relations"
        }
    def Activated(self):
        QtGui.QMessageBox.information(None, "De Laval Nozzle", "Pressure choke systems have not been implemented.\n\nThe generated nozzle is a solid body , so users can boolean substract it from another body on their own.\n\nThis is only a proof of concept to demonstrate the learning capabilities of my group and I\n\nDocumentation can be provided upon request")
        panel = DeLavalNozzleTaskPanel()
        FreeCADGui.Control.showDialog(panel)
#    def IsActive(self):
#       return FreeCAD.ActiveDocument is not Nones
FreeCADGui.addCommand('De_Laval_Nozzle', _DeLavalNozzle())
