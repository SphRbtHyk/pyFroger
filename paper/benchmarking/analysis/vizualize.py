import networkx as nx
import matplotlib.pyplot as plt

def parse_nested_structure(s):
    stack = []
    edges = []
    current_node = None
    i = 0
    n = len(s)

    while i < n:
        char = s[i]
        if char == '(':
            stack.append(current_node)
            current_node = None
            i += 1
        elif char == ')':
            if stack:
                parent = stack.pop()
                if current_node is not None and parent is not None:
                    edges.append((parent, current_node))
                current_node = parent
            i += 1
        elif char == ',':
            if current_node is not None and stack:
                parent = stack[-1]
                if parent is not None:
                    edges.append((parent, current_node))
                current_node = None
            i += 1
        else:
            if char.strip():
                current_node = char.strip()
                i += 1
            else:
                i += 1
    return edges

def plot_nested_graph(nested_structure):
    edges = parse_nested_structure(nested_structure)
    G = nx.DiGraph()
    G.add_edges_from(edges)

    pos = nx.spring_layout(G, seed=42)
    plt.figure(figsize=(12, 8))
    nx.draw(G, pos, with_labels=True, node_size=1000, node_color='skyblue',
            font_size=10, font_weight='bold', arrows=True)
    plt.title("Nested Structure as NetworkX Graph")
    plt.show()

nested_structure = "((((Ae,S),W),T,Ba,((Bb,Bd),Be),O,P,V,Ca),(F,(((E,Cd,C),(X,H)),N)),(Cb,((G,Cc),(R,Ab,Ce),(((L,K,A),M,Cf,B),Ac)),Z,Ad)))"
plot_nested_graph(nested_structure)