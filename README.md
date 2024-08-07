# MotifSpotter

A regular expressions supported motif searching script


## Getting Started

### Dependencies

This script required the [`regex`](https://github.com/mrabarnett/mrab-regex) and [`PyExcelerate`](https://github.com/kz26/PyExcelerate?tab=readme-ov-file) packages to run. Please run the following code to install the package before running the script.
```bash
$ pip install regex
$ pip install pyexcelerate
```
### Usage

```bash
motifspotter [-h] -f FILE -m MOTIF -e INTEGER -t {dna,rna,amino} -o {csv,excel} [-v]
```
> [!NOTE]
> As for Windows machine, the backslash appears in the file path (`\`) should replace with double backslash (`\\`) to avoid error.

> [!NOTE]
> In some cases, you might keep getting a warning about lack of dependence even after properly installing `regex`. This might because your python environment is mixed and packages was installed in a wrong path. This should be solve by running `pip3 install regex`. However this could possibly trigger another error, `error: externally-managed-environment`, because of python3.X's packages management policy. To avoid this error, you can just simply create a virtual environment and run the script from it.

```bash
$ cd $SCRIPTDIR
$ mkdir .venv
$ python3 -m venv .venv
$ source .venv/bin/activate
(venv.) $ python3 -m pip install regex
...
```

### Example

The following command is an example of searching for p53 response elements, allowing a maximum of three mismatches.
```bash
python3 motifspotter -f './p21.fasta' -m RRRCWWGYYYN{0,15}RRRCWWGYYY -e 3 -t dna -o excel
```


## Todo

- [x] Search specific motif of a given FASTA file 
- [x] Support regular expression
- [x] Support fuzzy searching
- [x] Mark unmatched base
- [x] Testing on multi running environment (Tested on Windows 11 and macOS 12.4)
- [ ] Reset span counting at every set of sequence
- [ ] Support static webpage output
- [ ] Get rid of double backslash


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code
of conduct, and the process for submitting pull requests to us.


## Authors

  - **Gaiki Ryo** - *The code* -
    [im-gaikiryo](https://github.com/im-gaikiryo)

## License

This project is licensed under the [MIT](LICENSE.md)
License - see the [LICENSE.md](LICENSE.md) file for
details

## Acknowledgments

  - [StackOverflow Thread](https://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string)

