"""
Layered lattice plot with per-node readings — pure matplotlib.

Tiers by set-size; within-tier order minimises edge crossings and edges
routed through nodes. Each node shows its set label plus a compact summary
of its readings (value x count). Reading-sets with no matching lattice node
are drawn detached with a dashed outline.
"""
import random
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


data = {
    '1': [[{'child': {'MS_C'}, 'parent': {'MS_D', 'MS_C'}}],
          [{'child': {'MS_D'}, 'parent': {'MS_D', 'MS_C'}}]],
    '2': [[{'child': {'MS_A', 'MS_B'}, 'parent': {'MS_A', 'MS_B', 'MS_D'}}],
          [{'child': {'MS_D', 'MS_C'}, 'parent': {'MS_A', 'MS_B', 'MS_D', 'MS_C'}}]],
}

readings = {
    frozenset({'MS_A', 'MS_B', 'MS_D'}): ['bon'],
    frozenset({'MS_C'}): ['vaillant', 'vit'],
    frozenset({'MS_A', 'MS_B'}): ['chevauche', 'trouva'],
    frozenset({'MS_D', 'MS_C'}): ['dirige'],
    frozenset({'MS_D'}): ['apperceut'],
}


GAP_X, GAP_Y = 3.2, 1.6


def label(s):
    return "{" + ", ".join(sorted(s)) + "}"


def summarize(vals):
    cnt = Counter(vals)
    return ", ".join(f"{k}\u00d7{cnt[k]}" if cnt[k] > 1 else k
                     for k in sorted(cnt))


# --- Parse relationships ---------------------------------------------------
edges, parents_of, node_size = set(), {}, {}
for groups in data.values():
    for group in groups:
        for rel in group:
            c, p = label(rel['child']), label(rel['parent'])
            edges.add((p, c))                         # arrow: parent -> child
            parents_of.setdefault(c, set()).add(p)
            node_size[c] = len(rel['child'])
            node_size[p] = len(rel['parent'])

red_nodes = {c for c, ps in parents_of.items() if len(ps) > 1}
edges = sorted(edges)

# --- Attach readings -------------------------------------------------------
reading_label, orphan_set = {}, set()
for fs, vals in readings.items():
    lbl = label(set(fs))
    reading_label[lbl] = summarize(vals)
    if lbl not in node_size:                          # reading w/o a node
        node_size[lbl] = len(fs)
        orphan_set.add(lbl)
if orphan_set:
    print("reading-only nodes (no relations):", orphan_set)

tiers = {}
for n, sz in node_size.items():
    tiers.setdefault(sz, []).append(n)
sizes = sorted(tiers)


# --- Crossing / through-node minimising ordering ---------------------------
def coords(order):
    pos = {}
    for sz in sizes:
        row = order[sz]
        for i, n in enumerate(row):
            pos[n] = ((i - (len(row) - 1) / 2) * GAP_X, sz * GAP_Y)
    return pos


def cost(order):
    pos = coords(order)

    def ccw(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

    cross = 0
    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            a, b = edges[i]
            c, d = edges[j]
            if len({a, b, c, d}) < 4:
                continue
            p1, p2, p3, p4 = pos[a], pos[b], pos[c], pos[d]
            if (ccw(p3, p4, p1) > 0) != (ccw(p3, p4, p2) > 0) and \
               (ccw(p1, p2, p3) > 0) != (ccw(p1, p2, p4) > 0):
                cross += 1

    through = 0
    for a, b in edges:
        (xa, ya), (xb, yb) = pos[a], pos[b]
        lo, hi = sorted((ya, yb))
        for sz in sizes:
            ylev = sz * GAP_Y
            if lo < ylev < hi:
                ex = xa + (xb - xa) * ((ylev - ya) / (yb - ya))
                for m in order[sz]:
                    if m not in (a, b) and abs(ex - pos[m][0]) < GAP_X * 0.45:
                        through += 1
    return cross * 10 + through


best_order, best_cost = None, float('inf')
for seed in range(30):
    rng = random.Random(seed)
    order = {sz: list(tiers[sz]) for sz in sizes}
    for sz in sizes:
        rng.shuffle(order[sz])
    improved = True
    while improved:
        improved = False
        for sz in sizes:
            row = order[sz]
            for i in range(len(row) - 1):
                base = cost(order)
                row[i], row[i + 1] = row[i + 1], row[i]
                if cost(order) < base:
                    improved = True
                else:
                    row[i], row[i + 1] = row[i + 1], row[i]
    c = cost(order)
    if c < best_cost:
        best_cost, best_order = c, {sz: list(order[sz]) for sz in sizes}

pos = coords(best_order)
print("red nodes:", red_nodes, "| cost:", best_cost)

# --- Draw ------------------------------------------------------------------
xs = [p[0] for p in pos.values()]
ys = [p[1] for p in pos.values()]
fig, ax = plt.subplots(figsize=(12, 7.5))
ax.set_xlim(min(xs) - 2.2, max(xs) + 2.2)
ax.set_ylim(min(ys) - 1.0, max(ys) + 0.9)
ax.axis('off')
fig.canvas.draw()
inv = ax.transData.inverted()
rend = fig.canvas.get_renderer()


def data_size(text):
    bb = text.get_window_extent(renderer=rend)
    (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
    return abs(x1 - x0), abs(y1 - y0)


boxes = {}
for n, (x, y) in pos.items():
    is_red = n in red_nodes
    is_orphan = n in orphan_set
    main = ax.text(x, y, n, ha='center', va='center', fontsize=10,
                   fontweight='bold', zorder=4,
                   color='white' if is_red else 'black')
    wm, hm = data_size(main)
    if n in reading_label:                            # stack readings below
        sub = ax.text(x, y, reading_label[n], ha='center', va='center',
                      fontsize=8, style='italic', zorder=4,
                      color='#f4f4f4' if is_red else '#16607a')
        ws, hs = data_size(sub)
        gap = 0.06
        total = hm + gap + hs
        main.set_y(y + total / 2 - hm / 2)
        sub.set_y(y - total / 2 + hs / 2)
        w, h = max(wm, ws), total
    else:
        w, h = wm, hm
    w, h = w + 0.34, h + 0.30
    boxes[n] = (w, h)

    if is_orphan:
        fc, ec, ls = '#dfe6e9', '#7f8c8d', (0, (4, 3))
    elif is_red:
        fc, ec, ls = '#e74c3c', '#1b2631', 'solid'
    else:
        fc, ec, ls = '#5dade2', '#1b2631', 'solid'
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle='round,pad=0.02,rounding_size=0.10', linewidth=1.4,
                 edgecolor=ec, facecolor=fc, linestyle=ls, zorder=2))

for p, c in edges:
    (x1, y1), (x2, y2) = pos[p], pos[c]
    ax.add_patch(FancyArrowPatch(
        (x1, y1 - boxes[p][1] / 2), (x2, y2 + boxes[c][1] / 2),
        arrowstyle='-|>', mutation_scale=16, color='#34495e',
        lw=1.5, shrinkA=2, shrinkB=2, zorder=1))

plt.title("Microservice subset lattice (parent \u2192 child)\n"
          "red = child with multiple parents   \u00b7   "
          "italic = readings (value\u00d7count)   \u00b7   dashed = no relations",
          fontsize=10.5)
plt.show()