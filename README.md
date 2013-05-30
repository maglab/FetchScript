# FetchScript

This is a basic python script to fetch details from multiple databases and return them as a simple python dict.

It can also translate in a limited fashion ID's from one biological database to another.

## Basic usage:

```
import fetchscript.fetch

f = fetch.FetchReference()
print f.fetchPubmed(23039964, withTerms=True)

f = fetch.FetchDetails()
f.fetchDetailsFromNucleotide('NM_000163.4')
```
