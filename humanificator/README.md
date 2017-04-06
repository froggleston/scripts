# Prerequisites:

You'll need Python 3, and to pip install nltk, feedparser and scholarly.

# Running:

``` python humanificator.py -q "Moore's Law in Single Cell Transcriptomics."```

# Fun:

If the tool finds references to the human genome, and the title doesn't contain 'human', it'll add ```, in humans.``` to the end of the title.

# Stuff to add:

As suggested by @jennifergardy on Twitter:

```We need one that looks to see if a paper title/abstract containing "genome" actually contains sequenced genome(s)!``` - Working on it! :)
