import seaborn as sb
import numpy as np

from matplotlib import pyplot as plt, lines

import argparse

from statistics import mean, median


parser = argparse.ArgumentParser(description="Generate plot from given depth data.")
parser.add_argument("input", help="Input file (text)")
parser.add_argument("title", help="The title of the figure")
parser.add_argument("output", help="Output file (image)")
parser.add_argument("-s", "--show", help="Show the result", action="store_true", default=False)
# parser.add_argument("-a", "--avg", dest="method", help="Compute the average per x bases (set with -n)")
# parser.add_argument("-m", "--med", dest="method", help="Compute the median per x bases (set with -n)")
# parser.add_argument("-p", "--per", type=int, dest="method", help="Compute the n-th percentile per x bases (set with -n)")
# parser.add_argument("-n", "--num", type=int, help="Set the amount of bases to plot on one point (default is 1)", default=1)
parser.add_argument("-n", nargs="*", metavar=('n', 'method'), default=[1], help="Set the amout of bases to plot on one point (default is 1) and set the computing method (currently supported are: mean, median, quantile n-th).")

args = parser.parse_args()


def parse(filename, n, method):
    """
    Depth-file parser.
    """
    depths = []
    refs = set()

    with open(filename, "r") as file:
        if n == 1:   # minor optimization
            for row in file:
                genome_id, pos, depth = row.split()
                
                refs.add(genome_id)
                if len(refs) > 1:
                    raise Exception("Only one genome contig is allowed!")

                depths.append(int(depth))

        else:
            depths_tmp = []
            for i, row in enumerate(file):
                genome_id, pos, depth = row.split()
                
                refs.add(genome_id)
                if len(refs) > 1:
                    raise Exception("Only one genome contig is allowed!")

                depths_tmp.append(int(depth))

                if i % n == 0:
                    if method == "mean":
                        depths.append(mean(depths_tmp))
                    elif method == "median":
                        depths.append(median(depths_tmp))
                    elif method == "quantile":
                        if len(args.n) < 3:
                            raise Exception("USAGE: -n basecount quantile n-th")
                        depths.append(np.quantile(np.array(depths_tmp), float(args.n[2])))
                    else:
                        raise Exception(f"Unknown computing function: {method}")

                    depths_tmp = []

    return depths


def plot_depth(filename, title, outputfile, n, method):
    
    depth = parse(filename, n, method)
    y_label = "Depth"
    x_label = "Genome position (bp)"

    sb.set(color_codes=True)
    fig = plt.figure()
    plt.title(title)
    ax = plt.subplot(111)

    sb_plot = sb.lineplot(x=range(0, len(depth)*n, n), y=depth)
    sb_plot.set(xlabel=x_label, ylabel=y_label)


    plt.savefig(outputfile)
    if args.show:
        plt.show()
    # plt.savefig(outputfile, dpi=400)
    plt.close()

    print("Done!")

def annotation_line( ax, xmin, xmax, y, text, ytext=0, linecolor='black', linewidth=1, fontsize=12 ):

    ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords='data', textcoords='data',
            arrowprops={'arrowstyle': '|-|', 'color':linecolor, 'linewidth':linewidth})
    ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords='data', textcoords='data',
            arrowprops={'arrowstyle': '<->', 'color':linecolor, 'linewidth':linewidth})

    xcenter = xmin + (xmax-xmin)/2
    if ytext==0:
        ytext = y + ( ax.get_ylim()[1] - ax.get_ylim()[0] ) / 20

    ax.annotate( text, xy=(xcenter,ytext), ha='center', va='center', fontsize=fontsize)


if __name__ == "__main__":
    #plot_depth("output-view.txt", "Simulated Medulloblastoma Sample", "plot.png")
    n = int(args.n[0])
    method = args.n[1] if len(args.n) > 1 else "mean"
    print(f"Generating read-depth plot \"{args.title}\". Graph per {n} base{'s using '+method if n>1 else ''} from {args.input} to {args.output}.")
    plot_depth(args.input, args.title, args.output, n, method)