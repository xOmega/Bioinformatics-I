import numpy as np


def find_nearest_neighbours(score_matrix, labels):
    min_val = np.float('inf')
    nearest_neighbours = []
    final = []
    rows = len(score_matrix)
    cols = len(score_matrix[0])
    for i in range(rows):
        for j in range(i + 1, cols):
            if score_matrix[i][j] < min_val:
                min_val = score_matrix[i][j]
                nearest_neighbours = [labels[i], labels[j]]
    return nearest_neighbours, min_val


def upgma(score_matrix, labels):
    rows = len(score_matrix)
    cols = len(score_matrix[0])
    labels = labels
    base_score_table = np.copy(score_matrix)
    cluster_seq = []
    while rows * cols > 1:
        # Find nearest neighbours
        nearest_neighbours, distance = find_nearest_neighbours(score_matrix, labels)
        print nearest_neighbours
        final = []
        # Cluster them together
        score_matrix = np.delete(score_matrix, (labels.index(nearest_neighbours[0]), labels.index(nearest_neighbours[1])), axis=0)
        score_matrix = np.delete(score_matrix, (labels.index(nearest_neighbours[0]), labels.index(nearest_neighbours[1])), axis=1)
        for i in nearest_neighbours:
            if type(i) != list:
                final.extend([i])
            else:
                final.extend(i)
        cluster_seq.append([[nearest_neighbours[0], nearest_neighbours[1]], distance])

        for i in nearest_neighbours:
            labels.remove(i)
        labels.append(final)

        if len(labels) == 1:
            return cluster_seq

        new_row = []
        flag = False
        for i in labels[:-1]:
            aggr = 0
            for j in labels[-1]:
                if type(i) == list:
                    flag = True
                    for k in i:
                        aggr += base_score_table[k][j]
                else:
                    aggr += base_score_table[i][j]
            if not flag:
                new_row.append(aggr/len(labels[-1]))
            else:
                new_row.append(aggr/(len(labels[-1]) + len(i)))

        if flag:
            flag = False

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


a = np.array(sample_score_matrix_2)
print upgma(a, range(len(a)))
