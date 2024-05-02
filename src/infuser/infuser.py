from .single_tree import single_tree

# @app.command()
# def singletree(tree_path: str, sample_file: str, output_dir: str, chrom_sizes: str, chromlist: List[str],\
#     res: int = 10000, subset: str = None, dist: int = 0, n_values: int = 9, min: float = -4, max: float = 4,\
#     column: int = 4, transform: List[str] = ["Z-score"], balance: bool = True, n_jobs: int = 4):
#     '''
#     Given a differentiation tree and data at its leaves, infer the data at the internal nodes of the tree.
#     '''
#     single_tree(tree_path, sample_file, output_dir, chrom_sizes, chromlist,\
#     res, subset, dist, n_values, min, max,\
#     column, transform, balance, n_jobs)

# if __name__ == "__main__":
#     app()

import click
@click.group()
def cli():
    pass

@cli.command()
def info():
    print("InfUSER verison 1.0.0 (https://github.com/AudreyBaguette/InfUSER)")

@cli.command()
@click.argument('treepath', nargs=1)
@click.argument('samplefile', nargs=1)
@click.argument('outdir', nargs=1)
@click.argument('chromsizes', nargs=1)
@click.argument('chromlist', nargs=1)
@click.option('-r', '--resolution', type=int, help='the resolution to consider', default=10000)
@click.option('-s', '--subset', type=str, help='the path to the file containing the regions to subset, if any', default='')
@click.option('-d', '--dist', type=int, help='the distance to consider. All interactions beyond that distance will be ignored. If set to 0, all interactions are kept', default=0)
@click.option('-nv', '--nvalues', type=int, help='the number of values that can be stored in the nodes', default=9)
@click.option('--min', type=float, help='the minimal Z-score value to consider', default=-4)
@click.option('--max', type=float, help='the maximal Z-score value to consider', default=4)
@click.option('-c', '--column', type=int, help='the column conting the score to consider. The first column is column 1. Ignored if the input files are .mcool files', default=4)
@click.option('-t', '--transform', type=list, help='the transformation(s) to apply to the matrix. "OE". "log1p" and "Z-score" are supported.', default=['Z-score'])
@click.option('-b', '--balance', type=bool, help='should the balanced weights be used', default=True)
@click.option('-nj', '--njobs', type=int, help='paralleliation of pixel computation, how many jobs should be run in parallel', default=4)
def singletree(treepath, samplefile, outdir, chromsizes, chromlist, resolution, subset, dist, nvalues, \
min, max, column, transform, balance, njobs):
    '''
    Run InfUSER with a single data type

    TREEPATH     the path to the file with the tree topology

    SAMPLEFILE   the path to the file with the paths to the samples

    OUTDIR       the path to the output directory

    CHROMSIZES   the path to the file containing the size (in bp) of each chromosome

    CHROMLIST    the names of the chromosomes to consider. This list is ignored if subset is not null
 
    '''
    chromlist = chromlist.split(',')
    single_tree(treepath, samplefile, outdir, chromsizes, chromlist, resolution, subset, dist, nvalues, \
        min, max, column, transform, balance, njobs)
