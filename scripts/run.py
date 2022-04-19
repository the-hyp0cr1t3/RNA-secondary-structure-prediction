#!/usr/bin/env python3

import sys
import os
import subprocess
import matplotlib.pyplot as plt
import networkx as nx

# OUT_PATH = "../build/run_output.txt"
# EXEC_CMD = [ "../bin/app", "-nc", "-o" , OUT_PATH, "-i" ]
EXEC_CMD = ["../bin/app"]

def getSequence(path):
    try:
        with open(path, 'r') as f:
            sequence = f.read()
    except:
        with open(os.path.join("../data", path), 'r') as f:
            sequence = f.read()

    return sequence

def getMatches(output):
    # return [(int(x.lstrip().split()[0]), int(x.lstrip().split()[1])) for x in output.stdout.decode().split("\n")[2:(2+int(output.stdout.decode().split("\n")[1].lstrip()))]]\

    decoded = output.stdout.decode()
    n = int(decoded.split("\n")[1].lstrip())
    return [(int(x.lstrip().split()[0]), int(x.lstrip().split()[1])) for x in decoded.split("\n")[2:2+n]]


def plotData(sequence, matches):
    plt.figure(figsize=(max(4, len(sequence)/8), max(4, len(sequence)/8)))

    node_colors = {
        "A": "#FF6B6B",
        "C": "#FFD93D",
        "G": "#6BCB77",
        "U": "#4D96FF"
    }

    G = nx.Graph()
    for i, c in enumerate(sequence):
        G.add_node(i+1, base=c, color=node_colors[c])
    for i in range(1, len(sequence)):
        G.add_edge(i, i+1, color="#00092C", style="-")
    for e in matches:
        G.add_edge(e[0], e[1], color="#FF5F00", style="--")

    nx.draw_kamada_kawai(
        G,
        labels=nx.get_node_attributes(G, "base"),
        edge_color=nx.get_edge_attributes(G, "color").values(),
        style=nx.get_edge_attributes(G, "style").values(),
        node_color=nx.get_node_attributes(G, "color").values(),
        font_color="#151D3B",
        font_weight="bold",
        linewidths=3
    )
    plt.gca().collections[0].set_edgecolor("#00092C")
    plt.gca().collections[0].set_linewidth(1)

    plt.show()

def main():
    output = subprocess.run(EXEC_CMD + [sys.argv[1]], capture_output=True)
    print(output.stdout.decode())

    print("Plotting graph...")
    plotData(getSequence(sys.argv[1]), getMatches(output))

if __name__ == "__main__":
    main()