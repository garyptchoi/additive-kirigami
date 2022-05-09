import numpy as np


def norm(x):
    return np.linalg.norm(x)


def normalize(x):
    return x / norm(x)


def append_row(mat, row):
    return np.append(mat, [row], axis=0)


def get_num_rows(mat):
    return mat.shape[0]


def is_even(x):
    return bool((x + 1) % 2)


def is_odd(x):
    return bool(x % 2)


def rotate_90(v, ccw):
    rotated_v = np.array((0.0, 0.0))
    rotated_v[0], rotated_v[1] = ccw * -v[1], ccw * v[0]
    return rotated_v


def cyclic(x, a):
    return np.roll(x, a, axis=0)


def parse_vertex(v):
    return 'v ' + ' '.join([str(coord) for coord in v]) + ' 0.0\n'


def parse_face(f):
    return 'f ' + '// '.join([str(ind + 1) for ind in f]) + '//\n'


def empty_list_of_lists(n):
    return [[] for _ in range(n)]


def rotation_matrix(angle):
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


def identity_matrix(n):
    return np.eye(n)


def rotation_matrix_3d(angle):
    return np.array([[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])


def rotation_matrix_homog(angle):
    zero_row = np.array([0.0, 0.0])
    homog_col = np.array([[0.0], [0.0], [1.0]])
    return np.hstack([np.vstack([rotation_matrix(angle), zero_row]), homog_col])


def translation_matrix_homog(tx, ty):
    return np.array([[1.0, 0.0, tx],
                     [0.0, 1.0, ty],
                     [0.0, 0.0, 1.0]])


def multiply_matrices(mats):
    if len(mats) == 2:
        return np.dot(mats[0], mats[1])
    else:
        mats[-2] = np.dot(mats[-2], mats[-1])
        mats.pop()
        return multiply_matrices(mats)


def rotate_points(points, origin, angle):

    if len(points.shape) == 1:
        onedim = True
        points = np.array([points])
    else:
        onedim = False

    num_points = len(points)

    if points.shape[1] == 3:
        rot_mat = rotation_matrix_3d(angle)
    else:
        rot_mat = rotation_matrix(angle)

    rotated_points = points - np.tile(origin, (num_points, 1))

    for i, point in enumerate(rotated_points):
        rotated_points[i] = np.array((np.matrix(rot_mat) * np.matrix(point).T).T)
    rotated_points += np.tile(origin, (num_points, 1))

    if onedim:
        rotated_points = rotated_points[0]

    return rotated_points


def planar_cross(a, b):

    return a[0] * b[1] - a[1] * b[0]


def calculate_angle(a, b, c):

    a = np.array(a, copy=True)
    b = np.array(b, copy=True)
    c = np.array(c, copy=True)

    ab_hat = normalize(b - a)
    ac_hat = normalize(c - a)

    x = np.dot(ab_hat, ac_hat)
    y = planar_cross(ab_hat, ac_hat)

    atan2_angle = np.arctan2(y, x)

    return atan2_angle % (2.0 * np.pi)


def shift_points(points, shift):

    return points + np.tile(shift, (len(points), 1))


def plot_structure(points, quads, linkages, ax):
       
    for i, quad in enumerate(quads):
        x = points[quad, 0]
        y = points[quad, 1]
        ax.fill(x, y, color=(1, 229/255, 204/255),edgecolor='k', linewidth=2, alpha = 0.8)

    ax.axis('off');
    ax.set_aspect('equal')
    

def deployment_linkage2matrix(i, j, num_linkage_rows, num_linkage_cols):

    if 0 <= i < num_linkage_rows and 0 <= j < num_linkage_cols:  # bulk

        matrix_row_ind = i * num_linkage_cols + j
        label = 'bulk'

    else:

        num_bulk_linkages = num_linkage_rows * num_linkage_cols

        num_boundary_linkages = [num_linkage_rows,
                                 num_linkage_cols,
                                 num_linkage_rows,
                                 num_linkage_cols]

        if j == -1 and 0 <= i < num_linkage_rows:  # left
            side_ind = 0
            bound_ind = i
            label = 'left'

        elif i == num_linkage_rows and 0 <= j < num_linkage_cols:  # bottom
            side_ind = 1
            bound_ind = j
            label = 'bottom'

        elif j == num_linkage_cols and 0 <= i < num_linkage_rows:  # right
            side_ind = 2
            bound_ind = num_linkage_rows - 1 - i
            label = 'right'

        elif i == -1 and 0 <= j < num_linkage_cols:  # top
            side_ind = 3
            bound_ind = num_linkage_cols - 1 - j
            label = 'top'

        else:  # linkage DNE
            return None

        num_other_boundary_linkages = sum(num_boundary_linkages[:side_ind])

        matrix_row_ind = num_bulk_linkages + num_other_boundary_linkages + bound_ind

    return matrix_row_ind, label


def write_obj(filename, points, quads):

    obj = open(filename, 'w')
    obj.write('# {} vertices, {} faces\n'.format(get_num_rows(points), get_num_rows(quads)))
    str_points = [parse_vertex(point) for point in points]
    obj.writelines(str_points)
    str_quads = [parse_face(quad) for quad in quads]
    obj.writelines(str_quads)
    obj.close()


def read_obj(filename):

    obj = open(filename, 'r')

    points = []
    faces = []

    for line in obj:

        first_char = line[0]

        if first_char == 'v':
            point = [float(_) for _ in line.split(' ')[1:-1]]
            points.append(point)

        elif first_char == 'f':
            face = [int(_) for _ in line.replace('//', '').split(' ')[1:]]
            faces.append(face)

        else:
            continue

    obj.close()

    return np.array(points), np.array(faces)


def main():
    print('reloading Utils')
    return
