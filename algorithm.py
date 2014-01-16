import numpy as np
from object import Object

# http://www.informatikbuero.com/downloads/Improved_Laplacian_Smoothing_of_Noisy_Surface_Meshes.pdf
def laplacian_hc_smooth(obj, times=1, alpha=0, beta=0.5):
    o = obj.v.copy()
    p = o.copy()
    for i in range(times):
        b = np.empty(o.shape, dtype=float)
        q = p.copy()
        for i in range(o.shape[0]):
            if obj.adj[i].shape[0]:
                p[i,:] = np.sum(q[obj.adj[i]], axis=0) / obj.adj[i].shape[0]
            b[i,:] = p[i,:] - (alpha * o[i,:] + (1 - alpha) * q[i,:])
        for i in range(o.shape[0]):
            if obj.adj[i].shape[0]:
                p[i,:] = p[i,:] - (beta * b[i,:] + (1 - beta) * np.sum(b[obj.adj[i]], axis=0) / obj.adj[i].shape[0])
    obj.v = p
    obj.update()
