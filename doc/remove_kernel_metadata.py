import glob
import nbformat


def remove_kernel_metadata(nb):
    for cell in nb.cells:
        if cell.metadata.get('kernel'):
            cell.metadata.kernel.pop('id', None)
            cell.metadata.kernel.pop('name', None)
    return nb


# Get a list of all the Jupyter notebook files in the directory
notebook_files = glob.glob('*.ipynb')


def removeK():
    # Register the preprocessor for each notebook file
    nbsphinx_preprocessors = [remove_kernel_metadata for _ in range(len(notebook_files))]
