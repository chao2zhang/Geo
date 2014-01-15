class Object:
    def __init__(self, v=[], t=[], f=[]):
        self.v = v
        self.t = t
        self.f = f

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