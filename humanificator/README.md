# Prerequisites:

You'll need Python 3, and to pip install nltk, bs4 (BeautifulSoup), feedparser and scholarly.

# Usage:

``` python humanificator.py -g -w genome-wordlist.txt -q "The complete genome sequence of a Neanderthal from the Altai Mountains"```

This will search for the article by its title, extract abstract and full content if possible, and search for words in your wordlist. With the ```-g``` option, this example will look for genome sequence accessions.

# Flags:

```-v``` - verbose

```-q``` - query title

```-w``` - wordlist to use for tokenising and phase construction

```-d``` - search for reference to humans (this is hard - it's not a great solution)

```-g``` - detect phrases like 'complete genome sequence' and see if there are any accessions listed

# Fun:

Using ```-d```, if the tool finds references to the human genome, and the title doesn't contain 'human', it'll add ```, in humans.``` to the end of the title.

# Stuff to add:

Better phase collocation functionality - the extraction of phases attributed to words in the wordlist isn't perfect.

# Suggestions completed:

As suggested by @jennifergardy on Twitter:

```We need one that looks to see if a paper title/abstract containing "genome" actually contains sequenced genome(s)!``` - Pretty much done! (2017-04-12)
