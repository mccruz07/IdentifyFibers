"""Index path bounding rectangles without changing geometric decisions."""
from collections import defaultdict


class PathBoxes:
    def __init__(self, cell_size=64):
        self.cell_size = cell_size
        self.cells = defaultdict(set)
        self.bounds = {}

    def _cells(self, box):
        y0, x0, y1, x1 = box
        for row in range(int(y0 // self.cell_size), int(y1 // self.cell_size) + 1):
            for col in range(int(x0 // self.cell_size), int(x1 // self.cell_size) + 1):
                yield row, col

    def insert(self, number, path):
        y0, x0 = path.min(axis=0)
        y1, x1 = path.max(axis=0)
        box = self.bounds[number] = (int(y0), int(x0), int(y1), int(x1))
        for cell in self._cells(box):
            self.cells[cell].add(number)

    def remove(self, number):
        for cell in self._cells(self.bounds.pop(number)):
            members = self.cells[cell]
            members.remove(number)
            if not members:
                del self.cells[cell]

    def neighbors(self, number, margin):
        y0, x0, y1, x1 = self.bounds[number]
        box = (y0-margin, x0-margin, y1+margin, x1+margin)
        candidates = set()
        for cell in self._cells(box):
            candidates.update(self.cells.get(cell, ()))
        candidates.discard(number)
        for other in sorted(candidates):
            a, b, c, d = self.bounds[other]
            if a <= box[2] and b <= box[3] and c >= box[0] and d >= box[1]:
                yield other
