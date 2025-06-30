import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(
    'Visualize the distribution of learned edges between residues.')
parser.add_argument('--num-residues', type=int, default=146,
                    help='Number of residues of the PDB.')
parser.add_argument('--windowsize', type=int, default=95,
                    help='window size')
parser.add_argument('--threshold', type=float, default=0.1,
                    help='threshold for plotting')
parser.add_argument('--dist-threshold', type=int, default=12,
                    help='threshold for shortest distance')
parser.add_argument('--filename', type=str, default='logs/out_probs_train.npy',
                    help='File name of the probs file.')
args = parser.parse_args()


def getEdgeResults(threshold=False):
    a = np.load(args.filename)
    b = a[:, :, 1]
    c = a[:, :, 2]
    d = a[:, :, 3]

    # There are four types of edges, eliminate the first type as the non-edge
    probs = b+c+d
    # For default residue number 77, residueR2 = 77*(77-1)=5852
    residueR2 = args.num_residues*(args.num_residues-1)
    probs = np.reshape(probs, (args.windowsize, residueR2))

    # Calculate the occurence of edges
    edges_train = probs/args.windowsize

    results = np.zeros((residueR2))
    for i in range(args.windowsize):
        results = results+edges_train[i, :]

    if threshold:
        # threshold, default 0.6
        index = results < (args.threshold)
        results[index] = 0

    # Calculate prob for figures
    edges_results = np.zeros((args.num_residues, args.num_residues))
    count = 0
    for i in range(args.num_residues):
        for j in range(args.num_residues):
            if not i == j:
                edges_results[i, j] = results[count]
                count += 1
            else:
                edges_results[i, j] = 0

    return edges_results


def getDomainEdges(edges_results, domainName):

    if domainName == 'loop':
        startLoc = 16
        endLoc = 30
    elif domainName == 'S1':
        startLoc = 39
        endLoc = 48
    elif domainName == 'S2':
        startLoc = 69
        endLoc = 79
    elif domainName == 'S3':
        startLoc = 79
        endLoc = 94
    if domainName == 'S4':
        startLoc = 105
        endLoc = 109
    elif domainName == 'S5':
        startLoc = 125
        endLoc = 132

    edges_results_loop = edges_results[16:30, startLoc:endLoc]   #用切片取构域，这里要根据实际情况进行修改
    edges_results_S1 = edges_results[39:48, startLoc:endLoc]
    edges_results_S2 = edges_results[69:79, startLoc:endLoc]
    edges_results_S3 = edges_results[79:94, startLoc:endLoc]
    edges_results_S4= edges_results[105:109, startLoc:endLoc]
    edges_results_S5= edges_results[125:132, startLoc:endLoc]

    edge_num_loop = edges_results_loop.sum(axis=0)
    edge_num_S1 = edges_results_Sα4.sum(axis=0)
    edge_num_S2 = edges_results_S1.sum(axis=0)
    edge_num_S3 = edges_results_Sβ4.sum(axis=0)
    edge_num_S4 = edges_results_Sβ6.sum(axis=0)
    edge_num_S5 = edges_results_S2.sum(axis=0)

    if domainName == 'loop':
        edge_average_loop = 0
    else:
        edge_average_loop = edge_num_loop.sum(axis=0)/(14*(endLoc-startLoc))
    if domainName == 'S1':
        edge_average_S1 = 0
    else:
        edge_average_S1 = edge_num_Sα4.sum(axis=0)/(9*(endLoc-startLoc))
    if domainName == 'S2':
        edge_average_S2 = 0
    else:
        edge_average_S2 = edge_num_S1.sum(axis=0)/(10*(endLoc-startLoc))
    if domainName == 'S3':
        edge_average_S3 = 0
    else:
        edge_average_S3 = edge_num_S3.sum(axis=0)/(15*(endLoc-startLoc)) #这里的数字代表这个结构域的残基数量，也要根据实际情况进行修改
    if domainName == 'S4':
        edge_average_S4 = 0
    else:
        edge_average_S4 = edge_num_Sβ6.sum(axis=0)/(4*(endLoc-startLoc))
    if domainName == 'S5':
        edge_average_S5= 0
    else: 
        edge_average_S5 = edge_num_S5.sum(axis=0)/(7*(endLoc-startLoc))
   
    edges_to_all = np.hstack((edge_average_loop, edge_average_S1,
                             edge_average_S2, edge_average_S3, edge_average_S4, edge_average_S5))
    return edges_to_all


# Load distribution of learned edges
edges_results_visual = getEdgeResults(threshold=True)
# Step 1: Visualize results
ax = sns.heatmap(edges_results_visual, linewidth=0.5,
                 cmap="Blues", vmax=1.0, vmin=0.0)
plt.savefig('logs/probs.png', dpi=600)
# plt.show()
plt.close()

# Step 2: Get domain specific results
# According to the distribution of learned edges between residues, we integrated adjacent residues as blocks for a more straightforward observation of the interactions.
# For example, the residues in SOD1 structure are divided into seven domains (β1, diml, disl, zl, β2, el, β3).

edges_results = getEdgeResults(threshold=False)
# SOD1 specific:
loop = getDomainEdges(edges_results, 'loop')
Sα4 = getDomainEdges(edges_results, 'S1')
S1 = getDomainEdges(edges_results, 'S2')
Sβ4 = getDomainEdges(edges_results, 'S3')
Sβ6 = getDomainEdges(edges_results, 'S4')
S2 = getDomainEdges(edges_results, 'S5')
edges_results = np.vstack((loop, S1, S2, S3, S4, S5))

# print(edges_results)
edges_results_T = edges_results.T
index = edges_results_T < (args.threshold)
edges_results_T[index] = 0

print('hello')
print(edges_results)

# 保存矩阵
np.savetxt('logs/matrix.txt', edges_results)
np.save('logs/matrix.npy', edges_results)

# Visualize
ax = sns.heatmap(edges_results_T, linewidth=1,
                 cmap="Blues", vmax=1.0, vmin=0.0)
ax.set_ylim([7, 0])
plt.savefig('logs/edges_domain.png', dpi=600)
# plt.show()
plt.close()
