"""Plot the reviewed leads against the inherited coarse future-water scenario."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import json

from .siting import HERE, load_bound, run


def main():
    spec = json.loads((HERE / 'siting_inputs.json').read_text())
    atlas = load_bound(spec, 'atlas')
    result = run()
    wet = atlas['water_label'] > 0
    world = np.roll(wet, wet.shape[1] // 2, axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={'width_ratios': [1.8, 1]})
    cmap = ListedColormap(['#ede7da', '#7ea9c2'])
    for ax in axes:
        ax.imshow(world, extent=(-180, 180, -90, 90), origin='upper', cmap=cmap,
                  vmin=0, vmax=1, interpolation='nearest', aspect='auto')
        ax.grid(alpha=.2, linewidth=.6)
        ax.set_xlabel('Moon-fixed longitude (degrees east)')
        ax.set_ylabel('Moon-fixed latitude (degrees north)')
        for n, row in enumerate(result['candidates'], 1):
            x = (row['longitude_east_deg'] + 180) % 360 - 180
            y = row['latitude_deg']
            lead = n <= 4
            ax.scatter(x, y, s=75 if lead else 24, c='#a43b23' if lead else '#333333',
                       edgecolors='white', linewidths=.7, zorder=3)
            if lead and (ax is axes[0] or (-70 < x < -35 and 0 < y < 50)):
                offsets = {1:(7, -11), 2:(-15, -12), 3:(5, 7), 4:(-15, 6)}
                ax.annotate(str(n), (x, y), xytext=offsets[n], textcoords='offset points',
                            fontsize=11, weight='bold', color='#672410')
    axes[0].set(xlim=(-180, 180), ylim=(-90, 90), title='Global context: four priority regions and five other leads')
    axes[0].set_xticks(np.arange(-180, 181, 60))
    axes[0].set_yticks(np.arange(-90, 91, 30))
    axes[1].set(xlim=(-70, -35), ylim=(0, 50), title='Western nearside detail')
    axes[1].annotate('South Marius', (-56.9322, 13.0188), xytext=(8, -18),
                     textcoords='offset points', fontsize=9, arrowprops={'arrowstyle':'-', 'lw':.7})
    fig.suptitle('Reconnaissance leads for four lava-tube cities', fontsize=15, weight='bold')
    fig.text(.05, .115, '1  Tranquillitatis     2  Marius Hills     3  West Rima Mairan     4  Schröter extension', fontsize=10)
    fig.text(.05, .067, 'Blue: inherited 28%-water atlas, 0.25° cells (~7.6 km north–south). Beige: ground above common water level.', fontsize=9)
    fig.text(.05, .028, 'Points mark pits or approximate geophysical leads; they establish neither cavity footprints nor approved city sites.', fontsize=9)
    fig.subplots_adjust(left=.06, right=.98, top=.84, bottom=.24, wspace=.24)
    output = HERE / 'results/siting_map.png'
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == '__main__':
    main()
