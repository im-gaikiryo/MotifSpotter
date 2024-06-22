# MotifSpotter

A regular expressions supported motif searching script


## Getting Started

```bash
motifspotter [-h] -f FILE -m MOTIF -t {dna,rna,amino} [-o OUTPUT] [-v]
```

### Prerequisites

This script required `Regex` package to run. Please run the following code to install the package before running the script.
```bash
$ pip install regex
```

## Todo

- [x] Search specific motif of a given FASTA file 
- [x] Support regular expression
- [] > Support fuzzy searching
- [] Mark unmatched base
- [] Reset span counting at every set of sequence
- [] Testing on multi running environment

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

