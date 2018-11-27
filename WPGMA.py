import numpy as np

species = ['Homo_sapiens', 'Dog', 'Zebrafish', 'Chicken', 'Rat', 'Mouse', 'Cattle', 'Chimpanzee']
mapping = {key: name for key, name in enumerate(species)}


def find_nearest_neighbours(score_matrix, labels):
    min_val = np.float('inf')
    nearest_neighbours = []
    rows = len(score_matrix)
    cols = len(score_matrix[0])
    for i in range(rows):
        for j in range(i + 1, cols):
            if score_matrix[i][j] < min_val:
                min_val = score_matrix[i][j]
                nearest_neighbours = [labels[i], labels[j]]
    return nearest_neighbours, min_val


def wpgma(score_matrix, labels):
    rows = len(score_matrix)
    cols = len(score_matrix[0])
    labels = labels
    cluster_seq = []
    while rows * cols > 1:
        # Find nearest neighbours
        nearest_neighbours, distance = find_nearest_neighbours(score_matrix, labels)
        # print nearest_neighbours

        cluster_seq.append([[nearest_neighbours[0], nearest_neighbours[1]], distance])

        neighbour_index_1 = labels.index(nearest_neighbours[0])
        neighbour_index_2 = labels.index(nearest_neighbours[1])

        old_labels = labels[::]

        for i in nearest_neighbours:
            labels.remove(i)
        labels.append(nearest_neighbours)

        if len(labels) == 1:
            return cluster_seq

        new_row = []
        for j in labels[:-1]:
            aggr = 0
            aggr += score_matrix[neighbour_index_1][old_labels.index(j)]
            aggr += score_matrix[neighbour_index_2][old_labels.index(j)]
            new_row.append(aggr / float(2))

        score_matrix = np.delete(score_matrix,
                                 (neighbour_index_1, neighbour_index_2), axis=0)
        score_matrix = np.delete(score_matrix,
                                 (neighbour_index_1, neighbour_index_2), axis=1)
        score_matrix = np.vstack((score_matrix, np.array(new_row)))
        new_row = np.array([np.append(new_row, 0)])
        score_matrix = np.hstack((score_matrix, new_row.T))
        rows = len(score_matrix)
        cols = len(score_matrix[0])


sample_score_matrix_1 = [[0, 180, 80, 200, 60],
                         [180, 0, 100, 40, 160],
                         [80, 100, 0, 120, 100],
                         [200, 40, 120, 0, 180],
                         [60, 160, 100, 180, 0]]

sample_score_matrix_2 = [[0, 17, 21, 31, 23],
                         [17, 0, 30, 34, 21],
                         [21, 30, 0, 28, 39],
                         [31, 34, 28, 0, 43],
                         [23, 21, 39, 43, 0]]

HTT_score_matrix = [[0.0000, 0.0917, 0.3500, 0.1893, 0.0996, 0.0993, 0.1196, 0.0094],
                    [0.0917, 0.0000, 0.3572, 0.1966, 0.1192, 0.1192, 0.1139, 0.0886],
                    [0.3500, 0.3572, 0.0000, 0.3072, 0.3527, 0.3572, 0.3705, 0.3455],
                    [0.1893, 0.1966, 0.3072, 0.0000, 0.1966, 0.1997, 0.2156, 0.1833],
                    [0.0996, 0.1192, 0.3527, 0.1966, 0.0000, 0.0250, 0.1393, 0.0951],
                    [0.0993, 0.1192, 0.3572, 0.1997, 0.0250, 0.0000, 0.1414, 0.0958],
                    [0.1196, 0.1139, 0.3705, 0.2156, 0.1393, 0.1414, 0.0000, 0.1203],
                    [0.0094, 0.0886, 0.3455, 0.1833, 0.0951, 0.0958, 0.1203, 0.0000]]

print mapping
a = np.array(HTT_score_matrix)
steps = wpgma(a, range(len(a)))
for i in range(len(steps)):
    print "Step %d -" % i, "nearest neighbors", steps[i][0], ",distance:", steps[i][1]