import re, argparse, sys, scholarly, urllib, feedparser, nltk
nltk.download('punkt')
nltk.download('wordnet')

from nltk import word_tokenize, sent_tokenize, corpus
from nltk.text import Text
from nltk.collocations import *
from nltk.corpus import wordnet
from Bio import Entrez

email = 'your.email.here@you.com'

wnl = nltk.WordNetLemmatizer()

def search(query):
    Entrez.email = email
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='20',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = email
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

## create a filter that consults a scored dictionary and a list of query words, from http://stackoverflow.com/questions/23479179/combining-filters-in-nltk-collocations
def create_filter_minfreq_inwords(scored, words, minfreq):
    def bigram_filter(w1, w2):
        return (w1 not in words and w2 not in words) and (
                (w1, w2) in scored and scored[w1, w2] <= minfreq)
    return bigram_filter

## collect word phases, from https://simply-python.com/2014/03/14/saving-output-of-nltk-text-concordance/
def word_phases(target_word, query_text, left_margin = 10, right_margin = 10):
    """
        Function to get all the phases that contain the target word in a text/passage tar_passage.
         
        str target_word, str tar_passage int left_margin int right_margin --> list of str
        left_margin and right_margin allocate the number of words/pununciation before and after target word
        Left margin will take note of the beginning of the text
    """
    ## Create list of tokens using nltk function
    #tokens = nltk.word_tokenize(abstract)
     
    ## Create the text of tokens
    #text = nltk.Text(tokens)
 
    ## Collect all the index or offset position of the target word
    c = nltk.ConcordanceIndex(query_text.tokens, key = lambda s: s.lower())
 
    ## Collect the range of the words that is within the target word by using text.tokens[start;end].
    ## The map function is use so that when the offset position - the target range < 0, it will be default to zero
    concordance_txt = ([query_text.tokens[list(map(lambda x: x-5 if (x-left_margin)>0 else 0,[offset]))[0]:offset+right_margin] for offset in c.offsets(target_word)])
                         
    ## join the sentences for each of the target phrase and return it
    return [''.join([x+' ' for x in con_sub]) for con_sub in concordance_txt]
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query')
    parser.add_argument('-w', '--wordfile')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    if args.wordfile is None:
        if args.verbose:
            print('No wordfile. Defaulting to human stuff')
        wordlist = ['GCRh37', 'GCRh38', 'human', 'genome', 'variant']
    else:
        wordlist = []
        if args.verbose:
            print('Wordfile specified. Loading...')
        with open(args.wordfile) as f:
            for t in set(f.read().splitlines()):
#                wordlist.append(wordnet.morphy(str.lower(t)))
                wordlist.append(wnl.lemmatize(t, 'n'))

    results = search(args.query)
    id_list = results['IdList']
    abstract = None
    meshterms = []

    if not id_list:
        print('Nothing in pubmed, looking in arXiv and Google Scholar...')
        # try arxiv
        url = 'http://export.arxiv.org/api/query'
        q = 'ti:"{0}"'.format(args.query)
        params = {'search_query':q, 'start':0, 'max_results':1}
        data = urllib.parse.urlencode(params).encode('ascii')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as response:
            res = response.read()
            feed = feedparser.parse(res)
            if not feed.entries:
                print('Not found in arXiv\n')
            else:
                abstract = feed.entries[0].summary
                if args.verbose:
                    print('Found in arXiv:\n' + abstract + '\n')

        # try scholar
        search_query = scholarly.search_pubs_query('%22'+args.query+'%22')
        for scholar_result in search_query:
            sr = scholar_result.fill()
            if not sr.bib:
                print('Not found in Google Scholar\n')
            else:
                if not sr.bib["abstract"]:
                    # ignore
                    continue
                abstract = scholar_result.bib["abstract"]
                if args.verbose:
                    print('Found Google Scholar records:\n' + abstract + '\n')
                break
    else:
        papers = fetch_details(id_list)
        for i, paper in enumerate(papers["PubmedArticle"]):
            if paper["MedlineCitation"]['MeshHeadingList']:
                for mesh in paper["MedlineCitation"]['MeshHeadingList']:
                    for term in nltk.word_tokenize(mesh["DescriptorName"]):
                        meshterms.append(wnl.lemmatize(str.lower(term), 'n'))
                if args.verbose:
                    print(meshterms)
            abstract = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]

    if not abstract:
        print('No abstracts available')
    else:
        tokens = nltk.word_tokenize(abstract)
        tokens.extend(meshterms)
        print(tokens)
        text = nltk.Text(tokens)
        
        all_phases = []
        for word in wordlist:
            phases = word_phases(word, text)
            if phases:
                all_phases.append(phases)

        if args.verbose:
            print('Located potential phases:\n{0}'.format(all_phases))

        # initialize finder object with the tokens
        finder = nltk.collocations.BigramCollocationFinder.from_words(tokens)

        # build a dictionary with bigrams and their frequencies
        bigram_measures = nltk.collocations.BigramAssocMeasures()
        scored = dict(finder.score_ngrams(bigram_measures.raw_freq))

        # create the filter...
        myfilter = create_filter_minfreq_inwords(scored, wordlist, 0.1)
        before_filter = list(finder.score_ngrams(bigram_measures.raw_freq))
        print('Before filter:\n', before_filter)

        # apply filter
        finder.apply_ngram_filter(myfilter)
        after_filter = list(finder.score_ngrams(bigram_measures.raw_freq))
        print('\nAfter filter:\n', after_filter)
        
        if not after_filter:
            print('No humans detected')
        else:
            print(args.query.rstrip('.') + ', in humans.')

