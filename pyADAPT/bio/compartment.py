class Compartment(object):
    def __init__(self, sbcmp):
        self.id = sbcmp.id
        self.name = sbcmp.name
        self.size = sbcmp.size
        self.parent = sbcmp.outside

class CmpTopo(object):
    def __init__(self, *args, **kwargs):
        self.dc = dict()
        for c in args:
            self.dc[c.id] = c

    def getParent(self, c):
        return self.dc[c.parent]

    def getChildren(self, c):
        r = []
        for _c in self.dc:
            if _c.parent == c.id:
                r.append(_c)
        return r


