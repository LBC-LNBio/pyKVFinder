import numpy
from _gridprocessing import igrid, fgrid, dgrid, cgrid

# FIXME: What is the best way to do this?
# Idea: Create + Fill + Surf

class Grid(object):

    def __init__(self, dx: int, dy: int, dz: int, kind: str):
        self._matrix = self.create(dx, dy, dz, kind)

    def create(self, dx, dy, dz, kind):
        if kind == "i":
            return igrid(dx * dy * dz).reshape(dx, dy, dz)
        elif kind == "f":
            return fgrid(dx * dy * dz).reshape(dx, dy, dz)
        elif kind == "d":
            return dgrid(dx * dy * dz).reshape(dx, dy, dz)
        elif kind == "c":
            return cgrid(dx * dy * dz).reshape(dx, dy, dz)
        else:
            raise TypeError("kind should be: i, f, d or c.")

    def fill(self):
        pass

        

