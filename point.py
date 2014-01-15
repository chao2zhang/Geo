class Point3f:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
    def __add__(self, other):
        if isinstance(other, Point3f):
            return Point3f(self.x + other.x, self.y + other.y, self.z + other.z)
        elif isinstance(other, (int, float)):
            return Point3f(self.x + other, self.y + other, self.z + other)
        return self

    def __sub__(self, other):
        if isinstance(other, Point3f):
            return Point3f(self.x - other.x, self.y - other.y, self.z - other.z)
        elif isinstance(other, (int, float)):
            return Point3f(self.x - other, self.y - other, self.z - other)
        return self

    def __mul__(self, other):
        if isinstance(other, Point3f):
            return Point3f(self.x * other.x, self.y * other.y, self.z * other.z)
        elif isinstance(other, (int, float)):
            return Point3f(self.x * other, self.y * other, self.z * other)
        return self

    def __div__(self, other):
        if isinstance(other, Point3f):
            return Point3f(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, (int, float)):
            return Point3f(self.x / other, self.y / other, self.z / other)
        return self

    def __radd__(self, other):
        return self + other

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return Point3f(-self.x, -self.y, -self.z)

    def __abs__(self):
        return pow(self.x * self.x + self.y * self.y + self.z * self.z, 0.5)

class Point2f:    
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
    def __add__(self, other):
        if isinstance(other, Point2f):
            return Point2f(self.x + other.x, self.y + other.y)
        elif isinstance(other, (int, float)):
            return Point2f(self.x + other, self.y + other)
        return self

    def __sub__(self, other):
        if isinstance(other, Point2f):
            return Point2f(self.x - other.x, self.y - other.y)
        elif isinstance(other, (int, float)):
            return Point2f(self.x - other, self.y - other)
        return self

    def __mul__(self, other):
        if isinstance(other, Point2f):
            return Point2f(self.x * other.x, self.y * other.y)
        elif isinstance(other, (int, float)):
            return Point2f(self.x * other, self.y * other)
        return self

    def __div__(self, other):
        if isinstance(other, Point2f):
            return Point2f(self.x / other.x, self.y / other.y)
        elif isinstance(other, (int, float)):
            return Point2f(self.x / other, self.y / other)
        return self

    def __radd__(self, other):
        return self + other

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return Point2f(-self.x, -self.y)

    def __abs__(self):
        return pow(self.x * self.x + self.y * self.y, 0.5)