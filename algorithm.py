from object import Object

# http://www.informatikbuero.com/downloads/Improved_Laplacian_Smoothing_of_Noisy_Surface_Meshes.pdf
def laplacian_hc_smooth(obj, times=1):
    o = obj.v
    p = o
    for i in range(times):
        q = p
        for i in range(len(o)):
            if len(obj.adj[i]):
                p[i] = reduce(lambda x, y: x + y, )
