For developers:

In order to generate new(er) documentation, follow the following steps:

    1. Make sure you have sphinx installed: `conda install sphinx`.
    2. In case of new .py files (modules):
        - Add module.rst file (copy another one and modify accordingly)
        - Add module to modules.rst
    3. Change version name in conf.py.
    4. Run `make html` for html docs or `make latexpdf` (make sure to 
       have LaTeX installed) 
    5. In case it gives error that 'sphinx-rtd-theme' is not installed, 
       install this with `pip install sphinx-rtd-theme`.
    6. Copy generated documentation to documentation folder. The (html) 
       home page is called index.html and the modules page modules.html.


For more information on Sphinx look [here](https://www.sphinx-doc.org/en/master/usage/quickstart.html).
