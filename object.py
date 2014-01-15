from point import Point3f, Point2f

class Object:
    def __init__(self, v=[], t=[], f=[]):
        self.v = v # vertex
        self.t = t # texture points
        self.f = f # face
    def add_neighbour(self, u, v):
        if v not in self.adj[u]:
            self.adj[u].append(v)
        if u not in self.adj[v]:
            self.adj[v].append(u)
    def update(self):
        self.adj = [[]] * len(v) # neibour indexes
        for f in self.f:
            add_neighbour(f[0][0], f[0][1])
            add_neighbour(f[0][1], f[0][2])
            add_neighbour(f[0][2], f[0][0])

def parse(file):
    o = Object()
    for line in file:
        e = line.split()
        if e[0] == 'v':
            o.v.append([float(x) for x in e[1:]])
        elif e[0] == 'vt':
            o.t.append([float(x) for x in e[1:]])
        elif e[0] == 'f':
            o.f.append([[int(x) - 1 for x in y.split('/')] for y in e[1:]])
            # remember to minus 1

import argparse
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='.obj file', type=str, default='')
    args = parser.parse_args()
    if not args.file:
        args.file = raw_input('Please input the filename:')
    f = open(args.file, 'r')
    parse(f)
    f.close()

if __name__ == '__main__':
    main()