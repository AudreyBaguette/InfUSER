import typer
from typing import List
from single_tree import single_tree

app = typer.Typer()

@app.command()
def info():
    print("InfUSER verison 1.0.0 (https://github.com/AudreyBaguette/InfUSER)")

@app.command()
def singletree(tree_path: str, sample_file: str, output_dir: str, chrom_sizes: str, chromlist: List[str],\
    res: int = 10000, subset: str = None, dist: int = 0, n_values: int = 9, min: float = -4, max: float = 4,\
    column: int = 4, transform: List[str] = ["Z-score"], balance: bool = True, n_jobs: int = 4):
    '''
    Given a differentiation tree and data at its leaves, infer the data at the internal nodes of the tree.
    '''
    single_tree(tree_path, sample_file, output_dir, chrom_sizes, chromlist,\
    res, subset, dist, n_values, min, max,\
    column, transform, balance, n_jobs)

if __name__ == "__main__":
    app()