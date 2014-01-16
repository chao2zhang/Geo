import numpy as np
import algorithm as algo
import time

class Object:
    def __init__(self, 
        v=np.zeros((0, 3), dtype=float), 
        t=np.zeros((0, 2), dtype=float), 
        f=[],
        text=[]):
        self.v = v # vertex
        self.t = t # texture points
        self.f = f # face
        self.text = text # other text lines
        self.update()
    def _clear(self):
        self.v = np.zeros((0, 3), dtype=float)
        self.t = np.zeros((0, 2), dtype=float)
        self.f = []
        self.text = []
    def _add_neighbour(self, u, v):
        self.adj[u].add(v)
        self.adj[v].add(u)
    def update(self):
        # neighbours
        # performance: replace 'list' with 'tuple'
        self.adj = [set() for x in range(self.v.shape[0])]
        for f in self.f:
            self._add_neighbour(f[0][0], f[1][0])
            self._add_neighbour(f[1][0], f[2][0])
            self._add_neighbour(f[2][0], f[0][0])
        self.adj = map(lambda a: np.array(tuple(a)), self.adj)
    def load(self, file):
        v = []
        t = []
        for line in file:
            e = line.split()
            if e[0] == 'v':
                v.append(map(float, e[1:]))
            elif e[0] == 'vt':
                t.append(map(float, e[1:]))
            elif e[0] == 'f':
                self.f.append([[int(x) - 1 for x in y.split('/')] for y in e[1:]])
                # remember to minus 1
            else:
                self.text.append(line)
        self.v = np.array(v)
        self.t = np.array(t)
        self.update()
    def save(self, file):
        file.writelines(self.text)
        for v in self.v:
            file.write('v %f %f %f\n' % tuple(v))
        for t in self.t:
            file.write('vt %f %f\n' % tuple(t))
        for f in self.f:
            file.write('f %i/%i %i/%i %i/%i\n' % (f[0][0] + 1, f[0][1] + 1, f[1][0] + 1, f[1][1] + 1, f[2][0] + 1, f[2][1] + 1))

import argparse
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='.obj file', type=str, default='')
    parser.add_argument('-o', '--output', help='.obj file', type=str, default='')
    args = parser.parse_args()
    if not args.input:
        args.input = raw_input('Please input the filename:')
    o = Object()
    t = time.time()
    print 'Loading...'
    with open(args.input, 'r') as f:
        o.load(f)
    print 'Loaded', time.time() - t

    t = time.time()
    print 'Executing laplacian HC smooth'
    algo.laplacian_hc_smooth(o,times=3)
    print 'Executed laplacian HC smooth', time.time() - t
    
    t = time.time()
    print 'Writing...'
    with open(args.output, 'w') as f:
        o.save(f)
    print 'Writed', time.time() - t

if __name__ == '__main__':
    main()