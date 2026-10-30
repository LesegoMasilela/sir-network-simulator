"""
SIR Disease Spread Simulator on Contact Networks
=================================================
Author : Lesego Masilela
GitHub : github.com/LesegoMasilela

Simulates the spread of infectious disease using the SIR (Susceptible-
Infected-Recovered) model on three types of contact networks:
  - Random (Erdos-Renyi)
  - Small-world (Watts-Strogatz)
  - Scale-free (Barabasi-Albert)

Each network type produces different epidemic dynamics, illustrating how
network topology influences disease propagation thresholds and outbreak size.
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import defaultdict
import random


# ── Simulation parameters ────────────────────────────────────────────────────

N          = 500    # Number of nodes (individuals)
BETA       = 0.3    # Transmission probability per contact per step
GAMMA      = 0.05   # Recovery probability per step
STEPS      = 100    # Number of time steps
SEED       = 42     # Random seed for reproducibility
INITIAL_INFECTED = 5  # Number of initially infected nodes


# ── Network constructors ──────────────────────────────────────────────────────

def make_random_network(n, p=0.02, seed=SEED):
    """Erdos-Renyi random graph. Each pair of nodes connected with probability p."""
    return nx.erdos_renyi_graph(n, p, seed=seed)


def make_smallworld_network(n, k=6, p=0.1, seed=SEED):
    """Watts-Strogatz small-world graph. High clustering, short path lengths."""
    return nx.watts_strogatz_graph(n, k, p, seed=seed)


def make_scalefree_network(n, m=3, seed=SEED):
    """Barabasi-Albert scale-free graph. Few highly connected hubs."""
    return nx.barabasi_albert_graph(n, m, seed=seed)


# ── SIR model ─────────────────────────────────────────────────────────────────

def run_sir(G, beta, gamma, initial_infected, steps, seed=SEED):
    """
    Run a discrete-time SIR epidemic on graph G.

    States:
      S = 0  (Susceptible)
      I = 1  (Infected)
      R = 2  (Recovered)

    At each time step:
      - Each infected node transmits to each susceptible neighbour
        with probability beta.
      - Each infected node recovers with probability gamma.

    Returns
    -------
    history : dict with keys 'S', 'I', 'R' — lists of counts over time.
    """
    random.seed(seed)
    np.random.seed(seed)

    # Initialise all nodes as susceptible
    state = {node: 0 for node in G.nodes()}

    # Seed initial infections
    infected_start = random.sample(list(G.nodes()), initial_infected)
    for node in infected_start:
        state[node] = 1

    history = defaultdict(list)

    for _ in range(steps):
        counts = {0: 0, 1: 0, 2: 0}
        for s in state.values():
            counts[s] += 1
        history['S'].append(counts[0])
        history['I'].append(counts[1])
        history['R'].append(counts[2])

        new_state = state.copy()

        for node in G.nodes():
            if state[node] == 1:
                # Try to infect susceptible neighbours
                for neighbour in G.neighbors(node):
                    if state[neighbour] == 0:
                        if random.random() < beta:
                            new_state[neighbour] = 1
                # Try to recover
                if random.random() < gamma:
                    new_state[node] = 2

        state = new_state

    # Record final counts
    counts = {0: 0, 1: 0, 2: 0}
    for s in state.values():
        counts[s] += 1
    history['S'].append(counts[0])
    history['I'].append(counts[1])
    history['R'].append(counts[2])

    return dict(history)


# ── Network statistics helper ─────────────────────────────────────────────────

def network_stats(G, name):
    """Print key structural properties of a network."""
    degrees = [d for _, d in G.degree()]
    print(f"\n{name}")
    print(f"  Nodes          : {G.number_of_nodes()}")
    print(f"  Edges          : {G.number_of_edges()}")
    print(f"  Avg degree     : {np.mean(degrees):.2f}")
    print(f"  Max degree     : {max(degrees)}")
    print(f"  Avg clustering : {nx.average_clustering(G):.4f}")
    try:
        print(f"  Avg path length: {nx.average_shortest_path_length(G):.4f}")
    except nx.NetworkXError:
        print(f"  Avg path length: N/A (disconnected graph)")


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_results(results, network_names):
    """
    Plot SIR dynamics for all three network types side by side,
    plus a summary comparison of peak infection and total outbreak size.
    """
    colours = {'S': '#378ADD', 'I': '#D85A30', 'R': '#1D9E75'}
    fig = plt.figure(figsize=(16, 10))
    fig.patch.set_facecolor('white')
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    time = range(STEPS + 1)

    # ── Row 1: SIR curves per network ────────────────────────────────────────
    for i, (name, history) in enumerate(zip(network_names, results)):
        ax = fig.add_subplot(gs[0, i])
        ax.plot(time, history['S'], color=colours['S'], linewidth=2, label='Susceptible')
        ax.plot(time, history['I'], color=colours['I'], linewidth=2, label='Infected')
        ax.plot(time, history['R'], color=colours['R'], linewidth=2, label='Recovered')
        ax.set_title(f"{name} Network", fontsize=13, fontweight='bold', pad=10)
        ax.set_xlabel("Time steps", fontsize=10)
        ax.set_ylabel("Number of individuals", fontsize=10)
        ax.set_xlim(0, STEPS)
        ax.set_ylim(0, N * 1.05)
        ax.legend(fontsize=9, loc='center right')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3)

        # Annotate peak infection
        peak_i = max(history['I'])
        peak_t = history['I'].index(peak_i)
        ax.annotate(
            f"Peak: {peak_i}",
            xy=(peak_t, peak_i),
            xytext=(peak_t + 5, peak_i + 15),
            fontsize=8,
            color=colours['I'],
            arrowprops=dict(arrowstyle='->', color=colours['I'], lw=1)
        )

    # ── Row 2: Comparison bar charts ─────────────────────────────────────────
    peak_infections = [max(h['I']) for h in results]
    total_recovered = [h['R'][-1] for h in results]
    bar_colours     = ['#7F77DD', '#D85A30', '#1D9E75']

    # Peak infection comparison
    ax2 = fig.add_subplot(gs[1, 0:2])
    bars = ax2.bar(network_names, peak_infections, color=bar_colours, width=0.5, edgecolor='white')
    ax2.set_title("Peak simultaneous infections by network type", fontsize=12, fontweight='bold')
    ax2.set_ylabel("Peak infected individuals", fontsize=10)
    ax2.set_ylim(0, N * 1.1)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, peak_infections):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 5, str(val),
                 ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Total outbreak size
    ax3 = fig.add_subplot(gs[1, 2])
    bars2 = ax3.bar(network_names, total_recovered, color=bar_colours, width=0.5, edgecolor='white')
    ax3.set_title("Total outbreak size\n(recovered by end)", fontsize=12, fontweight='bold')
    ax3.set_ylabel("Total recovered", fontsize=10)
    ax3.set_ylim(0, N * 1.1)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.grid(axis='y', alpha=0.3)
    for bar, val in zip(bars2, total_recovered):
        pct = round(val / N * 100)
        ax3.text(bar.get_x() + bar.get_width() / 2, val + 5, f"{pct}%",
                 ha='center', va='bottom', fontsize=11, fontweight='bold')

    fig.suptitle(
        f"SIR Disease Spread on Contact Networks  |  N={N}, β={BETA}, γ={GAMMA}",
        fontsize=14, fontweight='bold', y=1.01
    )

    plt.savefig("sir_results.png", dpi=150, bbox_inches='tight', facecolor='white')
    print("\nPlot saved to sir_results.png")
    plt.show()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("SIR Disease Spread Simulator on Contact Networks")
    print(f"Parameters: N={N}, β={BETA}, γ={GAMMA}, Steps={STEPS}")
    print("=" * 60)

    # Build networks
    print("\nBuilding networks...")
    networks = [
        make_random_network(N),
        make_smallworld_network(N),
        make_scalefree_network(N),
    ]
    names = ["Random (Erdos-Renyi)", "Small-World (Watts-Strogatz)", "Scale-Free (Barabasi-Albert)"]

    # Print structural stats
    print("\nNetwork structural properties:")
    print("-" * 40)
    for G, name in zip(networks, names):
        network_stats(G, name)

    # Run simulations
    print("\nRunning SIR simulations...")
    results = []
    for G, name in zip(networks, names):
        history = run_sir(G, BETA, GAMMA, INITIAL_INFECTED, STEPS)
        results.append(history)
        final_r  = history['R'][-1]
        peak_i   = max(history['I'])
        print(f"  {name}")
        print(f"    Peak infected : {peak_i} ({round(peak_i/N*100)}% of population)")
        print(f"    Total outbreak: {final_r} ({round(final_r/N*100)}% of population)")

    # Plot
    print("\nGenerating plots...")
    plot_results(results, ["Random", "Small-World", "Scale-Free"])

    print("\nDone.")


if __name__ == "__main__":
    main()
