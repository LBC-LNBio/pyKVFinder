# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This is the source code for the custom box operations by PyMOL pyKVFinder Tools.
Changes in this file are not advised, as it controls interactions with
pyKVFinder and PyQt5.
"""

import toml
import pymol
from math import pi, sin, cos
from pymol import cmd

from PyQt5.QtCore import QCoreApplication
from PyQt5.QtWidgets import QPushButton, QDoubleSpinBox, QMessageBox, QMainWindow


class Box(object):

    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0

    def _set_box(
        self,
        button_draw_box: QPushButton,
        button_redraw_box: QPushButton,
        min_x: QDoubleSpinBox,
        max_x: QDoubleSpinBox,
        min_y: QDoubleSpinBox,
        max_y: QDoubleSpinBox,
        min_z: QDoubleSpinBox,
        max_z: QDoubleSpinBox,
        angle1: QDoubleSpinBox,
        angle2: QDoubleSpinBox,
        padding: QDoubleSpinBox,
    ) -> None:
        """
        Create box coordinates, enable 'Delete Box' and 'Redraw Box' buttons
        and call _draw_box function.

        Parameters
        ----------
        button_draw_box : QPushButton
            The button to draw the box.
        button_redraw_box : QPushButton
            The button to redraw the box.
        min_x : QDoubleSpinBox
            The minimum x value of the box.
        max_x : QDoubleSpinBox
            The maximum x value of the box.
        min_y : QDoubleSpinBox
            The minimum y value of the box.
        max_y : QDoubleSpinBox
        min_z : QDoubleSpinBox
            The minimum z value of the box.
        max_z : QDoubleSpinBox
            The maximum z value of the box.
        angle1 : QDoubleSpinBox
            The first angle value of the box.
        angle2 : QDoubleSpinBox
            The second angle value of the box.
        padding : QDoubleSpinBox
            The padding value of the box.
        """

        # Delete Box object in PyMOL
        if "box" in cmd.get_names("all"):
            cmd.delete("box")
        # Get dimensions of selected residues
        selection = "sele"
        if selection in cmd.get_names("selections"):
            ([minx, miny, minz], [maxx, maxy, maxz]) = cmd.get_extent(selection)
        else:
            ([minx, miny, minz], [maxx, maxy, maxz]) = cmd.get_extent("")

        # Get center of each dimension (x, y, z)
        self.x = (minx + maxx) / 2
        self.y = (miny + maxy) / 2
        self.z = (minz + maxz) / 2

        # Set Box variables in interface
        min_x.setValue(round(self.x - (minx - padding.value()), 1))
        max_x.setValue(round((maxx + padding.value()) - self.x, 1))
        min_y.setValue(round(self.y - (miny - padding.value()), 1))
        max_y.setValue(round((maxy + padding.value()) - self.y, 1))
        min_z.setValue(round(self.z - (minz - padding.value()), 1))
        max_z.setValue(round((maxz + padding.value()) - self.z, 1))
        angle1.setValue(0)
        angle2.setValue(0)

        # Setting background box values
        self.min_x_set = min_x.value()
        self.max_x_set = max_x.value()
        self.min_y_set = min_y.value()
        self.max_y_set = max_y.value()
        self.min_z_set = min_z.value()
        self.max_z_set = max_z.value()
        self.angle1_set = angle1.value()
        self.angle2_set = angle2.value()
        self.padding_set = padding.value()

        # Draw box
        self._draw_box(min_x, max_x, min_y, max_y, min_z, max_z, angle1, angle2)

        # Enable/Disable buttons
        button_draw_box.setEnabled(False)
        button_redraw_box.setEnabled(True)
        min_x.setEnabled(True)
        min_y.setEnabled(True)
        min_z.setEnabled(True)
        max_x.setEnabled(True)
        max_y.setEnabled(True)
        max_z.setEnabled(True)
        angle1.setEnabled(True)
        angle2.setEnabled(True)

    def _draw_box(
        self,
        min_x: QDoubleSpinBox,
        max_x: QDoubleSpinBox,
        min_y: QDoubleSpinBox,
        max_y: QDoubleSpinBox,
        min_z: QDoubleSpinBox,
        max_z: QDoubleSpinBox,
        angle1: QDoubleSpinBox,
        angle2: QDoubleSpinBox,
    ) -> None:
        """
        Draw box in PyMOL interface.

        Parameters
        ----------
        min_x : QDoubleSpinBox
            The minimum x value of the box.
        max_x : QDoubleSpinBox
            The maximum x value of the box.
        min_y : QDoubleSpinBox
            The minimum y value of the box.
        max_y : QDoubleSpinBox
        min_z : QDoubleSpinBox
            The minimum z value of the box.
        max_z : QDoubleSpinBox
            The maximum z value of the box.
        angle1 : QDoubleSpinBox
            The first angle value of the box.
        angle2 : QDoubleSpinBox
            The second angle value of the box.
        """
        # Convert angle
        a1 = (angle1.value() / 180.0) * pi
        a2 = (angle2.value() / 180.0) * pi

        # Get positions of box vertices
        # P1
        x1 = (
            -min_x.value() * cos(a2)
            - (-min_y.value()) * sin(a1) * sin(a2)
            + (-min_z.value()) * cos(a1) * sin(a2)
            + self.x
        )

        y1 = -min_y.value() * cos(a1) + (-min_z.value()) * sin(a1) + self.y

        z1 = (
            min_x.value() * sin(a2)
            + min_y.value() * sin(a1) * cos(a2)
            - min_z.value() * cos(a1) * cos(a2)
            + self.z
        )

        # P2
        x2 = (
            max_x.value() * cos(a2)
            - (-min_y.value()) * sin(a1) * sin(a2)
            + (-min_z.value()) * cos(a1) * sin(a2)
            + self.x
        )

        y2 = (-min_y.value()) * cos(a1) + (-min_z.value()) * sin(a1) + self.y

        z2 = (
            (-max_x.value()) * sin(a2)
            - (-min_y.value()) * sin(a1) * cos(a2)
            + (-min_z.value()) * cos(a1) * cos(a2)
            + self.z
        )

        # P3
        x3 = (
            (-min_x.value()) * cos(a2)
            - max_y.value() * sin(a1) * sin(a2)
            + (-min_z.value()) * cos(a1) * sin(a2)
            + self.x
        )

        y3 = max_y.value() * cos(a1) + (-min_z.value()) * sin(a1) + self.y

        z3 = (
            -(-min_x.value()) * sin(a2)
            - max_y.value() * sin(a1) * cos(a2)
            + (-min_z.value()) * cos(a1) * cos(a2)
            + self.z
        )

        # P4
        x4 = (
            (-min_x.value()) * cos(a2)
            - (-min_y.value()) * sin(a1) * sin(a2)
            + max_z.value() * cos(a1) * sin(a2)
            + self.x
        )

        y4 = (-min_y.value()) * cos(a1) + max_z.value() * sin(a1) + self.y

        z4 = (
            -(-min_x.value()) * sin(a2)
            - (-min_y.value()) * sin(a1) * cos(a2)
            + max_z.value() * cos(a1) * cos(a2)
            + self.z
        )

        # P5
        x5 = (
            max_x.value() * cos(a2)
            - max_y.value() * sin(a1) * sin(a2)
            + (-min_z.value()) * cos(a1) * sin(a2)
            + self.x
        )

        y5 = max_y.value() * cos(a1) + (-min_z.value()) * sin(a1) + self.y

        z5 = (
            (-max_x.value()) * sin(a2)
            - max_y.value() * sin(a1) * cos(a2)
            + (-min_z.value()) * cos(a1) * cos(a2)
            + self.z
        )

        # P6
        x6 = (
            max_x.value() * cos(a2)
            - (-min_y.value()) * sin(a1) * sin(a2)
            + max_z.value() * cos(a1) * sin(a2)
            + self.x
        )

        y6 = (-min_y.value()) * cos(a1) + max_z.value() * sin(a1) + self.y

        z6 = (
            (-max_x.value()) * sin(a2)
            - (-min_y.value()) * sin(a1) * cos(a2)
            + max_z.value() * cos(a1) * cos(a2)
            + self.z
        )

        # P7
        x7 = (
            (-min_x.value()) * cos(a2)
            - max_y.value() * sin(a1) * sin(a2)
            + max_z.value() * cos(a1) * sin(a2)
            + self.x
        )

        y7 = max_y.value() * cos(a1) + max_z.value() * sin(a1) + self.y

        z7 = (
            -(-min_x.value()) * sin(a2)
            - max_y.value() * sin(a1) * cos(a2)
            + max_z.value() * cos(a1) * cos(a2)
            + self.z
        )

        # P8
        x8 = (
            max_x.value() * cos(a2)
            - max_y.value() * sin(a1) * sin(a2)
            + max_z.value() * cos(a1) * sin(a2)
            + self.x
        )

        y8 = max_y.value() * cos(a1) + max_z.value() * sin(a1) + self.y

        z8 = (
            (-max_x.value()) * sin(a2)
            - max_y.value() * sin(a1) * cos(a2)
            + max_z.value() * cos(a1) * cos(a2)
            + self.z
        )

        # Create box object
        pymol.stored.list = []
        if "box" in cmd.get_names("all"):
            cmd.iterate("box", "stored.list.append((name, color))", quiet=1)
        list_color = pymol.stored.list
        cmd.delete("box")
        if len(list_color) > 0:
            for item in list_color:
                at_name = item[0]
                at_c = item[1]
                cmd.set_color(at_name + "color", cmd.get_color_tuple(at_c))
        else:
            for at_name in [
                "v2",
                "v3",
                "v4",
                "v5",
                "v6",
                "v7",
                "v8",
                "v1x",
                "v1y",
                "v1z",
                "v2x",
                "v3y",
                "v4z",
            ]:
                cmd.set_color(at_name + "color", [0.86, 0.86, 0.86])

        # Create vertices
        cmd.pseudoatom("box", name="v2", pos=[x2, y2, z2], color="v2color")
        cmd.pseudoatom("box", name="v3", pos=[x3, y3, z3], color="v3color")
        cmd.pseudoatom("box", name="v4", pos=[x4, y4, z4], color="v4color")
        cmd.pseudoatom("box", name="v5", pos=[x5, y5, z5], color="v5color")
        cmd.pseudoatom("box", name="v6", pos=[x6, y6, z6], color="v6color")
        cmd.pseudoatom("box", name="v7", pos=[x7, y7, z7], color="v7color")
        cmd.pseudoatom("box", name="v8", pos=[x8, y8, z8], color="v8color")

        # Connect vertices
        cmd.select("vertices", "(name v3,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v5,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v3,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v6,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v7,v8)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("box", name="v1x", pos=[x1, y1, z1], color="red")
        cmd.pseudoatom("box", name="v2x", pos=[x2, y2, z2], color="red")
        cmd.select("vertices", "(name v1x,v2x)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("box", name="v1y", pos=[x1, y1, z1], color="forest")
        cmd.pseudoatom("box", name="v3y", pos=[x3, y3, z3], color="forest")
        cmd.select("vertices", "(name v1y,v3y)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("box", name="v4z", pos=[x4, y4, z4], color="blue")
        cmd.pseudoatom("box", name="v1z", pos=[x1, y1, z1], color="blue")
        cmd.select("vertices", "(name v1z,v4z)")
        cmd.bond("vertices", "vertices")
        cmd.delete("vertices")

    def _delete_box(
        self,
        button_draw_box: QPushButton,
        button_redraw_box: QPushButton,
        min_x: QDoubleSpinBox,
        max_x: QDoubleSpinBox,
        min_y: QDoubleSpinBox,
        max_y: QDoubleSpinBox,
        min_z: QDoubleSpinBox,
        max_z: QDoubleSpinBox,
        angle1: QDoubleSpinBox,
        angle2: QDoubleSpinBox,
    ) -> None:
        """
        Delete box object, disable 'Delete Box' and 'Redraw Box' buttons and
        enable 'Draw Box' button.

        Parameters
        ----------
        button_draw_box : QPushButton
            The button to draw the box.
        button_redraw_box : QPushButton
            The button to redraw the box.
        min_x : QDoubleSpinBox
            The minimum x value of the box.
        max_x : QDoubleSpinBox
            The maximum x value of the box.
        min_y : QDoubleSpinBox
            The minimum y value of the box.
        max_y : QDoubleSpinBox
        min_z : QDoubleSpinBox
            The minimum z value of the box.
        max_z : QDoubleSpinBox
            The maximum z value of the box.
        angle1 : QDoubleSpinBox
            The first angle value of the box.
        angle2 : QDoubleSpinBox
            The second angle value of the box.
        """
        from pymol import cmd

        # Reset all box variables
        self.x = 0
        self.y = 0
        self.z = 0

        # Delete Box and Vertices objects in PyMOL
        cmd.delete("vertices")
        cmd.delete("box")

        # Set Box variables in the interface
        min_x.setValue(0.0)
        max_x.setValue(0.0)
        min_y.setValue(0.0)
        max_y.setValue(0.0)
        min_z.setValue(0.0)
        max_z.setValue(0.0)
        angle1.setValue(0)
        angle2.setValue(0)

        # Change state of buttons in the interface
        button_draw_box.setEnabled(True)
        button_redraw_box.setEnabled(False)
        min_x.setEnabled(False)
        min_y.setEnabled(False)
        min_z.setEnabled(False)
        max_x.setEnabled(False)
        max_y.setEnabled(False)
        max_z.setEnabled(False)
        angle1.setEnabled(False)
        angle2.setEnabled(False)

    def _redraw_box(
        self,
        min_x: QDoubleSpinBox,
        max_x: QDoubleSpinBox,
        min_y: QDoubleSpinBox,
        max_y: QDoubleSpinBox,
        min_z: QDoubleSpinBox,
        max_z: QDoubleSpinBox,
        angle1: QDoubleSpinBox,
        angle2: QDoubleSpinBox,
        padding: QDoubleSpinBox,
    ) -> None:
        """
        Redraw box in PyMOL interface. If selection is provided, the box will
        be drawn around the selected residues. If no selection is provided,
        the box will be drawn around the entire structure.

        Parameters
        ----------
        button_draw_box : QPushButton
            The button to draw the box.
        button_redraw_box : QPushButton
            The button to redraw the box.
        min_x : QDoubleSpinBox
            The minimum x value of the box.
        max_x : QDoubleSpinBox
            The maximum x value of the box.
        min_y : QDoubleSpinBox
            The minimum y value of the box.
        max_y : QDoubleSpinBox
        min_z : QDoubleSpinBox
            The minimum z value of the box.
        max_z : QDoubleSpinBox
            The maximum z value of the box.
        angle1 : QDoubleSpinBox
            The first angle value of the box.
        angle2 : QDoubleSpinBox
            The second angle value of the box.
        padding : QDoubleSpinBox
            The padding value of the box.
        """

        # Provided a selection
        if "sele" in cmd.get_names("selections"):
            # Get dimensions of selected residues
            ([minx, miny, minz], [maxx, maxy, maxz]) = cmd.get_extent("sele")

            if (
                min_x.value() != self.min_x_set
                or max_x.value() != self.max_x_set
                or min_y.value() != self.min_y_set
                or max_y.value() != self.max_y_set
                or min_z.value() != self.min_z_set
                or max_z.value() != self.max_z_set
                or angle1.value() != self.angle1_set
                or angle2.value() != self.angle2_set
            ):
                self.min_x_set = min_x.value()
                self.max_x_set = max_x.value()
                self.min_y_set = min_y.value()
                self.max_y_set = max_y.value()
                self.min_z_set = min_z.value()
                self.max_z_set = max_z.value()
                self.angle1_set = angle1.value()
                self.angle2_set = angle2.value()
            # Padding or selection altered
            else:
                # Get center of each dimension (x, y, z)
                self.x = (minx + maxx) / 2
                self.y = (miny + maxy) / 2
                self.z = (minz + maxz) / 2

                # Set background box values
                self.min_x_set = (
                    round(self.x - (minx - padding.value()), 1)
                    + min_x.value()
                    - self.min_x_set
                )
                self.max_x_set = (
                    round((maxx + padding.value()) - self.x, 1)
                    + max_x.value()
                    - self.max_x_set
                )
                self.min_y_set = (
                    round(self.y - (miny - padding.value()), 1)
                    + min_y.value()
                    - self.min_y_set
                )
                self.max_y_set = (
                    round((maxy + padding.value()) - self.y, 1)
                    + max_y.value()
                    - self.max_y_set
                )
                self.min_z_set = (
                    round(self.z - (minz - padding.value()), 1)
                    + min_z.value()
                    - self.min_z_set
                )
                self.max_z_set = (
                    round((maxz + padding.value()) - self.z, 1)
                    + max_z.value()
                    - self.max_z_set
                )
                self.angle1_set = 0 + angle1.value()
                self.angle2_set = 0 + angle2.value()
                self.padding_set = padding.value()
        # Not provided a selection
        else:
            if (
                min_x.value() != self.min_x_set
                or max_x.value() != self.max_x_set
                or min_y.value() != self.min_y_set
                or max_y.value() != self.max_y_set
                or min_z.value() != self.min_z_set
                or max_z.value() != self.max_z_set
                or angle1.value() != self.angle1_set
                or angle2.value() != self.angle2_set
            ):
                self.min_x_set = min_x.value()
                self.max_x_set = max_x.value()
                self.min_y_set = min_y.value()
                self.max_y_set = max_y.value()
                self.min_z_set = min_z.value()
                self.max_z_set = max_z.value()
                self.angle1_set = angle1.value()
                self.angle2_set = angle2.value()

            if self.padding_set != padding.value():
                # Prepare dimensions without old padding
                minx = self.padding_set - self.min_x_set
                maxx = self.max_x_set - self.padding_set
                miny = self.padding_set - self.min_y_set
                maxy = self.max_y_set - self.padding_set
                minz = self.padding_set - self.min_z_set
                maxz = self.max_z_set - self.padding_set

                # Get center of each dimension (x, y, z)
                self.x = (minx + maxx) / 2
                self.y = (miny + maxy) / 2
                self.z = (minz + maxz) / 2

                # Set background box values
                self.min_x_set = round(self.x - (minx - padding.value()), 1)
                self.max_x_set = round((maxx + padding.value()) - self.x, 1)
                self.min_y_set = round(self.y - (miny - padding.value()), 1)
                self.max_y_set = round((maxy + self.padding.value()) - self.y, 1)
                self.min_z_set = round(self.z - (min_z - padding.value()), 1)
                self.max_z_set = round((maxz + padding.value()) - self.z, 1)
                self.angle1_set = angle1.value()
                self.angle2_set = angle2.value()
                self.padding_set = padding.value()

        # Set Box variables in the interface
        min_x.setValue(self.min_x_set)
        max_x.setValue(self.max_x_set)
        min_y.setValue(self.min_y_set)
        max_y.setValue(self.max_y_set)
        min_z.setValue(self.min_z_set)
        max_z.setValue(self.max_z_set)
        angle1.setValue(self.angle1_set)
        angle2.setValue(self.angle2_set)

        # Redraw box
        self._draw_box(min_x, max_x, min_y, max_y, min_z, max_z, angle1, angle2)

    def _create_box_file(self, filename: str) -> dict:
        """
        Create a box file (box.<basename>.pdb) with the parameters from the GUI.

        Returns
        -------
        dict
            The box parameters.
        """
        # Get box parameters
        min_x = self.min_x_set
        max_x = self.max_x_set
        min_y = self.min_y_set
        max_y = self.max_y_set
        min_z = self.min_z_set
        max_z = self.max_z_set
        angle1 = self.angle1_set
        angle2 = self.angle2_set

        # Convert angle
        angle1 = (angle1 / 180.0) * pi
        angle2 = (angle2 / 180.0) * pi

        # Get positions of box vertices
        # P1
        x1 = (
            -min_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + self.x
        )

        y1 = -min_y * cos(angle1) + (-min_z) * sin(angle1) + self.y

        z1 = (
            min_x * sin(angle2)
            + min_y * sin(angle1) * cos(angle2)
            - min_z * cos(angle1) * cos(angle2)
            + self.z
        )

        # P2
        x2 = (
            max_x * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + self.x
        )

        y2 = (-min_y) * cos(angle1) + (-min_z) * sin(angle1) + self.y

        z2 = (
            (-max_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P3
        x3 = (
            (-min_x) * cos(angle2)
            - max_y * sin(angle1) * sin(angle2)
            + (-min_z) * cos(angle1) * sin(angle2)
            + self.x
        )

        y3 = max_y * cos(angle1) + (-min_z) * sin(angle1) + self.y

        z3 = (
            -(-min_x) * sin(angle2)
            - max_y * sin(angle1) * cos(angle2)
            + (-min_z) * cos(angle1) * cos(angle2)
            + self.z
        )

        # P4
        x4 = (
            (-min_x) * cos(angle2)
            - (-min_y) * sin(angle1) * sin(angle2)
            + max_z * cos(angle1) * sin(angle2)
            + self.x
        )

        y4 = (-min_y) * cos(angle1) + max_z * sin(angle1) + self.y

        z4 = (
            -(-min_x) * sin(angle2)
            - (-min_y) * sin(angle1) * cos(angle2)
            + max_z * cos(angle1) * cos(angle2)
            + self.z
        )

        # Create points
        p1 = [x1,y1,z1]
        p2 = [x2,y2,z2]
        p3 = [x3,y3,z3]
        p4 = [x4,y4,z4]
        box = {"p1": p1, "p2": p2, "p3": p3, "p4": p4}

        with open(filename, "w") as f:
            toml.dump(box, f)

        return filename

    def _help(self, mainwindow: QMainWindow) -> None:
        """
        Display help information about the box operations.
        """
        text = QCoreApplication.translate(
            "pyKVFinder",
            '<html><head/><body><p align="justify"><span style=" font-weight:600; text-decoration: underline;">Box Adjustment mode:</span></p><p align="justify">- Create a selection (optional);</p><p align="justify">- Define a <span style=" font-weight:600;">Padding</span> (optional);</p><p align="justify">- Click on <span style=" font-weight:600;">Draw Box</span> button.</p><p align="justify"><br/><span style="text-decoration: underline;">Customize your <span style=" font-weight:600;">box</span></span>:</p><p align="justify">- Change one item at a time (e.g. <span style=" font-style:italic;">Padding</span>, <span style=" font-style:italic;">Minimum X</span>, <span style=" font-style:italic;">Maximum X</span>, ...);</p><p align="justify">- Click on <span style=" font-weight:600;">Redraw Box</span> button.<br/></p><p><span style=" font-weight:400; text-decoration: underline;">Delete </span><span style=" text-decoration: underline;">box</span><span style=" font-weight:400; text-decoration: underline;">:</span></p><p align="justify">- Click on <span style=" font-weight:600;">Delete Box</span> button.<br/></p><p align="justify"><span style="text-decoration: underline;">Colors of the <span style=" font-weight:600;">box</span> object:</span></p><p align="justify">- <span style=" font-weight:600;">Red</span> corresponds to <span style=" font-weight:600;">X</span> axis;</p><p align="justify">- <span style=" font-weight:600;">Green</span> corresponds to <span style=" font-weight:600;">Y</span> axis;</p><p align="justify">- <span style=" font-weight:600;">Blue</span> corresponds to <span style=" font-weight:600;">Z</span> axis.</p></body></html>',
            None,
        )
        help_information = QMessageBox(mainwindow)
        help_information.setText(text)
        help_information.setWindowTitle("Help")
        help_information.setStyleSheet("QLabel{min-width:500 px;}")
        help_information.exec_()
