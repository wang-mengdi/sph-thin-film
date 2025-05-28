import colorpy.thinfilm as thinfilm
import colorpy.illuminants as illuminants
import numpy as np
import math

illuminant = illuminants.get_illuminant_D65()
illuminants.scale_illuminant (illuminant, 9.50)
n1 = 1.000293
n2 = 1.333
n3 = 1.000293
tf = thinfilm.thin_film(n1, n2, n3, 0)
filename = "table.npy"
table = np.load(filename)

def get_color(h):
    h = np.array(h)
    h = h.reshape((-1,8))
    start = h[:,:3]
    end = h[:,4:7]
    #print ("start is", start.shape)
    #print ("end is", end.shape)
    h = 2e8 * 1./2000. * np.linalg.norm((end-start), axis = 1)
    #print ("h is", h)
    h -= np.min(h)
    h /= np.max(h)
    h *= 1000
    #print (h)

    colors = np.ones((h.shape[0], 4))
    for i in range(colors.shape[0]):
        tf = thinfilm.thin_film(1.000293, 1.333, 1.000293, h[i])
        colors[i, :3] = tf.illuminated_color(illuminant)

    #print (colors)
    #print (colors)
    return colors

def get_color_single(h):
    tf.thickness_nm = h
    tf.phase_factor = -2.0 * tf.thickness_nm * 2.0 * math.pi * n2
    return tf.illuminated_color(illuminant)

def precompute_table():
    num = int (tf.max_thickness_nm)-1
    table = np.zeros((num, 3))
    for i in range(num):
        print ("i: ", i)
        table[i, :] = get_color_single(i)
    print (table.shape)
    np.save(filename, table)

def get_color_lookup(h):
    #print("what is h: ", h)
    if h < int (tf.max_thickness_nm)-1:
        return table[int(h)]
    else:
        return np.ones((1,3)).flatten()

def precompute_table_bin():
    #num = int (tf.max_thickness_nm)-1
    num = 1000
    table = np.zeros((num, 3))
    for i in range(num):
        print ("i: ", i)
        table[i, :] = get_color_single(i)
    print (table.shape)
    print (table.dtype)
    filename = "color_table.dat"
    fileobj = open(filename, mode='wb')
    table.tofile(fileobj)
    fileobj.close



if __name__ == "__main__":
    precompute_table_bin()
    # for i in range(1000):
    #     color = get_color_lookup(i)
    #     print (color)
