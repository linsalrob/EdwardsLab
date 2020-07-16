from networkx import Graph

class mydepth(Graph):

    def __init__(self, node, graph):
        self.node = node
        self.graph = graph

from networkx import DiGraph

class anotherDepth(DiGraph):

    def __init__(self, node, graph):
        self.node = node
        self.graph = graph

