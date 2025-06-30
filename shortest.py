import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# 精细同心圆布局，基于大小进行定位并确保所有节点都在内部，调整层密度
def refined_concentric_layout_final_adjustment(G, sizes, scale=1):
    pos = {}
    sorted_nodes = sorted(G.nodes, key=lambda n: sizes[n], reverse=True)
    num_layers = int(np.ceil(np.sqrt(len(sorted_nodes))))

    layer = 1
    count = 0
    print(len(sorted_nodes))
    while count < len(sorted_nodes) - 11:  # 确保外层有足够的节点
        nodes_in_layer = 6 * layer  # 减少内层节点以使其更稀疏
        nodes = sorted_nodes[count:count + nodes_in_layer]
        angle_step = 2 * np.pi / len(nodes)
        for i, node in enumerate(nodes):
            angle = i * angle_step
            r = scale * layer
            pos[node] = (r * np.cos(angle), r * np.sin(angle))
        count += nodes_in_layer
        layer += 1

   # 将剩余的节点放在最外层，使其更紧凑
    remaining_nodes = sorted_nodes[count:]
    print(remaining_nodes)
    if len(remaining_nodes) > 0:  # 检查剩余节点数量是否大于 0
        angle_step = 2 * np.pi / len(remaining_nodes)
        for i, node in enumerate(remaining_nodes):
            angle = i * angle_step
            r = scale * (layer - 0.8)  # 调整半径使最外层更紧凑
            pos[node] = (r * np.cos(angle), r * np.sin(angle))
    else:
        # 当剩余节点数量为 0 时，可将其集中在一个默认位置，例如中心
        center_x = 0
        center_y = 0
        for node in remaining_nodes:
            pos[node] = (center_x, center_y)

    return pos

# 自定义函数在边缘上绘制箭头
def draw_arrows(G, pos, path_edges, ax, arrow_color):
    for edge in path_edges:
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        # 调整坐标使箭头出现在节点外部
        arrowprops = dict(arrowstyle='->',
                          color=arrow_color,
                          lw=5,  # 加粗线条和箭头
                          alpha=0.8,  # 使箭头和目标节点颜色一致
                          shrinkA=20, shrinkB=20,  # 增加缩减量以避免箭头嵌入节点
                          connectionstyle='arc3,rad=0.1',
                          mutation_scale = 20)
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=arrowprops)

# 加载邻接矩阵
file_path = './logs/NAD_interactions.npy'
matrix = np.load(file_path)

# 从邻接矩阵创建图对象
G = nx.from_numpy_array(matrix)

# 计算度中心性并转换为节点大小
degree_centrality = nx.degree_centrality(G)
node_sizes = [v * 1500 for v in degree_centrality.values()]  # 调整节点大小

# 生成基于大小定位的精细同心圆布局，确保所有节点都在内部并调整层密度
pos = refined_concentric_layout_final_adjustment(G, degree_centrality, scale=3)

# 统一所有节点的大背景颜色
node_colors = ['#eaeae6' for _ in range(len(G.nodes))]  # 使用浅绿色作为背景颜色

# 固定起点与多个终点
loop_range = [22]
S1_range = [43]
S2_range = [71]
S3_range = [90]
S4_range = [107]
S5_range = [127]

# 定义颜色
loop_color = '#A24C4C'     # 柔和暗红
S1_color   = '#B58AB2'     # 淡紫
S2_color   = '#95C7A6'     # 柔绿
S3_color   = '#90BDE6'     # 浅蓝
S4_color   = '#E7B97A'     # 浅橙
S5_color   = '#90D3D3'     # 淡青蓝

# 同时建议将非重点节点颜色调得更柔和些：
node_colors = ['#dddddd' for _ in range(len(G.nodes))]  # 替代 '#eaeae6'

def compute_node_weights(G):
    return {node: sum(weight for _, _, weight in G.edges(node, data='weight', default=1)) for node in G.nodes()}

node_weights = compute_node_weights(G)

# 节点信息打印
def print_node_info(name, node_list):
    print(f"\n{name}区域节点信息：")
    for node in node_list:
        print(f"节点 {node + 1} - 度中心性: {degree_centrality[node]}, 权重: {node_weights[node]}")

print_node_info("loop", loop_range)
print_node_info("S1", S1_range)
print_node_info("S2", S2_range)
print_node_info("S3", S3_range)
print_node_info("S4", S4_range)
print_node_info("S5", S5_range)

def find_best_path(start_nodes, end_nodes):
    all_paths = []
    path_probabilities = []
    for start_node in start_nodes:
        for end_node in end_nodes:
            path = nx.shortest_path(G, source=start_node, target=end_node, weight='weight')
            path_probability = np.prod([G[u][v]['weight'] for u, v in zip(path[:-1], path[1:])])
            all_paths.append(path)
            path_probabilities.append(path_probability)
    best_path_index = np.argmax(path_probabilities)
    return all_paths[best_path_index], all_paths, path_probabilities

# 路径计算
start_nodes = loop_range
targets = [
    (S1_range, S1_color),
    (S2_range, S2_color),
    (S3_range, S3_color),
    (S4_range, S4_color),
    (S5_range, S5_color),
]

best_paths = []
all_path_infos = []

for s_range, color in targets:
    best_path, all_paths, path_probs = find_best_path(start_nodes, s_range)
    total_prob = sum(path_probs)
    all_path_infos.append((all_paths, path_probs, total_prob, s_range))
    best_paths.append((best_path, color))

# 打印路径信息
for idx, (all_paths, path_probs, total_prob, target_range) in enumerate(all_path_infos):
    print(f"\n路径到 S{idx+1}：")
    for path, prob in sorted(zip(all_paths, path_probs), key=lambda x: x[1], reverse=True):
        print(f"路径: {[node + 1 for node in path]}, 概率: {prob / total_prob:.4f}")

# 绘图
fig, ax = plt.subplots(figsize=(16, 16))
enlarged_nodes = set([path[-1] for path, _ in best_paths])
node_sizes = [base_node_sizes[n] * 2 if n in enlarged_nodes else base_node_sizes[n] for n in G.nodes]

nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                       edgecolors='white', alpha=0.8, ax=ax, linewidths=2)
nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray', ax=ax)

# 路径箭头绘制
for path, color in best_paths:
    path_edges = list(zip(path, path[1:]))
    draw_arrows(G, pos, path_edges, ax, arrow_color=color)

# 特殊点绘制
special_nodes = [
    (loop_range, loop_color),
    (S1_range, S1_color),
    (S2_range, S2_color),
    (S3_range, S3_color),
    (S4_range, S4_color),
    (S5_range, S5_color),
]

for nodelist, color in special_nodes:
    nx.draw_networkx_nodes(G, pos, nodelist=nodelist,
                           node_size=[node_sizes[n] for n in nodelist],
                           node_color=color,
                           edgecolors='black', linewidths=2, ax=ax)

# 节点标签
labels = {node: str(node + 1) for node in G.nodes()}
nx.draw_networkx_labels(G, pos, labels=labels,
                        font_size=12, font_color='black', font_weight='bold',
                        ax=ax, bbox=dict(facecolor='none', edgecolor='none', alpha=0.5))

ax.set_facecolor('#f7f7f7')
plt.title('', fontsize=17)
plt.savefig('NAD_shortest_paths_all_targets.png', dpi=300, bbox_inches='tight')
plt.show()
