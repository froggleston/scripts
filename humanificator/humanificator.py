import re, argparse, sys, scholarly, urllib, feedparser, nltk, bs4, string
nltk.download('punkt')
nltk.download('wordnet')
nltk.download('averaged_perceptron_tagger')
nltk.download('maxent_ne_chunker')
nltk.download('words')
nltk.download('stopwords')

from nltk import word_tokenize, sent_tokenize, corpus
from nltk.text import Text
from nltk import chunk
from nltk.chunk import *
from nltk.chunk.util import *
from nltk.chunk.regexp import *
from nltk.tree import Tree
from nltk.collocations import *
from nltk.corpus import wordnet, stopwords
from Bio import Entrez
from bs4 import BeautifulSoup

email = 'your.email.here@you.com'

wnl = nltk.WordNetLemmatizer()

translator = str.maketrans('', '', string.punctuation)

refseq_pattern = re.compile(r'\b([A-Z]{2}_[\d]+)')
insdc_pattern = re.compile(r'\b([SED]R[APRSXZ]\d{6,7})')
#genbank_pattern = re.compile(r'\b([A-Z]{1}[0-9]{5}|[A-Z]{2}[0-9]{6}|[A-Z]{4}[0-9]{8,9}|[A-Z]{5}[0-9]{7})(\.[0-9]{1,3})*')
genbank_pattern = re.compile(r'\b([HNTRWUBIGS]{1}[0-9]{5}|(?:A[AIWFY]{1}|B[EFGIMQU]{1}|C[ABDFKNOVX]{1}|D[NRTVWYQ]{1}|E[BCEGHLSVWXYFU]{1}|F[CDEFGKLJ]{1}|G[DEHORTWQU]{1}|H[OSMQ]{1}|J[GKZFNQX]{1}|K[CFJMPRTUXY]{1})[0-9]{6}|[A-Z]{4}[0-9]{8,9}|[A-Z]{5}[0-9]{7}\.[0-9]{1,3})')
fallback_json = re.compile(r'\"text\":\"([A-Z0-9]+)\"')

def grab_pmid(pmid=4304705):
    print("Grabbing {0}".format(pmid))
#    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={0}'.format(pmid)
    url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/{0}/'.format(pmid)
    headers = { 'User-Agent' : 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:12.0) Gecko/20100101 Firefox/12.0' }
    req = urllib.request.Request(url=url, data=b'None', headers=headers)
    with urllib.request.urlopen(req) as response:
        html = response.read()
        soup = BeautifulSoup(html, "lxml")
        ## warning - all this is NCBI specific HTML munging
        # grab article div
        article = soup.find('div', attrs={'class':'article'})
        
        # scrub refs
        article.find('div', attrs={'id':'reference-list'}).decompose()

        # scrub article citation header
        article.find('div', attrs={'class':'fm-sec'}).decompose()

        # scrub figures and tables
        article.find('div', attrs={'class':'fig'}).decompose()
        article.find('div', attrs={'class':'table-wrap'}).decompose()
        
        # remove scripts, inline styling, headings and inline anchors/hrefs
        for script in article(["script", "style", "h1", "h2", "h3", "h4", "a"]):
            script.decompose()
        text = article.get_text()

        # break into lines and remove leading and trailing space on each
        chunks = (line.strip() for line in text.splitlines())

        ## this seemed to get rid of some trailing full-stops and whitespace helpful for tokenising
        # break multi-headlines into a line each
#        chunks = (phrase.strip() for line in lines for phrase in line.split("  "))

        # drop blank lines
        return '\n'.join(chunk for chunk in chunks if chunk)

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
def create_bigram_filter_minfreq_inwords(scored, words, minfreq):
    def bigram_filter(w1, w2):
        return (w1 not in words and w2 not in words) and (
                (w1, w2) in scored and scored[w1, w2] <= minfreq)
    return bigram_filter

def create_trigram_filter_minfreq_inwords(scored, words, minfreq):
    def trigram_filter(w1, w2, w3):
        return (w1 not in words and w2 not in words and w3 not in words) and (
                (w1, w3, w3) in scored and scored[w1, w2, w3] <= minfreq)
    return trigram_filter

## collect word phases, from https://simply-python.com/2014/03/14/saving-output-of-nltk-text-concordance/
def word_phases(target_word, query_text, left_margin = 10, right_margin = 10):
    """
        Function to get all the phases that contain the target word in a text/passage tar_passage.
         
        str target_word, str tar_passage int left_margin int right_margin --> list of str
        left_margin and right_margin allocate the number of words/pununciation before and after target word
        Left margin will take note of the beginning of the text
    """
    ## Collect all the index or offset position of the target word
    c = nltk.ConcordanceIndex(query_text.tokens, key = lambda s: s.lower())
 
    ## Collect the range of the words that is within the target word by using text.tokens[start;end].
    ## The map function is use so that when the offset position - the target range < 0, it will be default to zero
    concordance_txt = ([query_text.tokens[list(map(lambda x: x-5 if (x-left_margin)>0 else 0,[offset]))[0]:offset+right_margin] for offset in c.offsets(target_word)])
                         
    ## join the sentences for each of the target phrase and return it
    return [''.join([x+' ' for x in con_sub]) for con_sub in concordance_txt]

def extract_phases(tokens, wordlist):
    all_phases = []
    text = nltk.Text(tokens)
    for word in wordlist:
        phases = word_phases(word, text)
        if phases:
            all_phases.append(phases)
    return all_phases

def chunk_entities(tokens):
    ## produce tags
    tagged = nltk.pos_tag(tokens)

    ## identify named entities, and don't expand NE types
    return nltk.chunk.ne_chunk(tagged, binary=True)

def find_and_filter_bigrams(tokens):
    # initialize finder object with the tokens
    finder = nltk.collocations.BigramCollocationFinder.from_words(tokens)

    # build a dictionary with bigrams and their frequencies
    bigram_measures = nltk.collocations.BigramAssocMeasures()
    scored = dict(finder.score_ngrams(bigram_measures.raw_freq))

    # create the filter...
    myfilter = create_bigram_filter_minfreq_inwords(scored, wordlist, 0.1)
#    before_filter = list(finder.score_ngrams(bigram_measures.raw_freq))
#    if args.verbose:
#        print('Before filter:\n', before_filter)

    # apply filter
    finder.apply_ngram_filter(myfilter)
    
    # remove empty
    finder.apply_word_filter(lambda w: w in (''))

    # remove low freq
    finder.apply_freq_filter(2)

    #after_filter = list(finder.score_ngrams(bigram_measures.raw_freq))
    best_filter = sorted(finder.nbest(bigram_measures.raw_freq, 5))
    if args.verbose:
        #print('\nAfter filter:\n', after_filter)
        print('\nTop bigrams:\n', best_filter)
    return best_filter

def find_and_filter_trigrams(tokens):
    # initialize finder object with the tokens
    finder = nltk.collocations.TrigramCollocationFinder.from_words(tokens)

    # build a dictionary with bigrams and their frequencies
    trigram_measures = nltk.collocations.TrigramAssocMeasures()
    scored = dict(finder.score_ngrams(trigram_measures.raw_freq))

    # create the filter...
    myfilter = create_trigram_filter_minfreq_inwords(scored, wordlist, 0.1)

    # apply filter
    finder.apply_ngram_filter(myfilter)

    # remove empty
    finder.apply_word_filter(lambda w: w in (''))
    
    # remove low freq
    finder.apply_freq_filter(2)

    #after_filter = list(finder.score_ngrams(trigram_measures.raw_freq))
    best_filter = sorted(finder.nbest(trigram_measures.raw_freq, 5))
    if args.verbose:
        #print('\nAfter filter:\n', after_filter)
        print('\nTop trigrams:\n', best_filter)
    return best_filter

def score_chunks(chunkparser, tree, trace=0):
    chunkscore = chunk.ChunkScore()
    
    tokens = tree.leaves()
    gold = chunkparser.parse(Tree('S', tokens), trace)
    chunkscore.score(gold, tree)

    a = chunkscore.accuracy()*100
    p = chunkscore.precision()*100
    r = chunkscore.recall()*100
    f = chunkscore.f_measure()*100
    i = chunkscore.incorrect()
    m = chunkscore.missed()

    if args.verbose:
        print('Testing: {0}'.format(tree))
        print('Guessed: {0}'.format(chunkscore.guessed()))
        print('\n/'+('='*75)+'\\')
        print('Scoring', chunkparser)
        print(('-'*77))
        print('Accuracy: %5.1f%%' % (a), ' '*4, end=' ')
        print('Precision: %5.1f%%' % (p), ' '*4, end=' ')
        print('Recall: %5.1f%%' % (r), ' '*6, end=' ')
        print('F-Measure: %5.1f%%' % (f))
        print('Missing: {0} -> {1}'.format(len(m), m))
        print('Incorrect: {0} -> {1}'.format(len(i), i))

    return (a, p, r, f, i, m)

def detectotext(bigrams, phases, title, wordlist):
    if not bigrams:
       print('No ngrams available to detect query text')
    else:
       for bigram in bigrams:
           for word in wordlist:
               if str.lower(bigram[0][0]) == word or str.lower(bigram[0][1]) == word:
                   print(title.rstrip('.') + ', in {0}.'.format(word))
                   return True
               if args.verbose:
                   print("{0} not detected. This paper might be relevant to stuff other than {0}.".format(word))
    return False

def detectogenome(bigrams, phases, title):
    orig_title_tokens = nltk.word_tokenize(title)
    lc_title_tokens = nltk.word_tokenize(str.lower(title))

    #grammar = "NE: {<JJ>*<NN|JJ|NNP><NN|NNP><IN|TO>*}"
    #chunk_rule = ChunkRule("<.*>+", "Chunk all")
    another_chunk = ChunkRule("<JJ>*<NN|NNP>?<NN|NNP>+<IN|OF>*", "[Complete] genome sequence ish")
    chink_rule = ChinkRule("<CC|VB|IN|\.>", "Chink some")
    another_chink = ChinkRule("<NNP><NN>", "Chink some more")
    #split_rule = SplitRule("<JJ>*<NN|NNP><NN|NNP>", "Split some")

#    cp = nltk.RegexpParser(grammar)
    cp = nltk.RegexpChunkParser([another_chunk, chink_rule, another_chink], chunk_label='NP')
   
    o_chunks = chunk_entities(orig_title_tokens)
    o_tree = cp.parse(o_chunks)
    stats = score_chunks(cp, o_tree)
    
    maybe = False

    ## check abstract text
    for phase in phases:
        for sent in phase:
            if args.verbose:
                print('\nPhase: {0}'.format(sent))

            phase_tokens = nltk.word_tokenize(sent)
            phase_chunks = chunk_entities(phase_tokens)
            phase_tree = cp.parse(phase_chunks)
            stats = score_chunks(cp, phase_tree)

            if stats[0] > 50 and stats[1] > 50 and stats[2] > 25 and stats[3] > 0:
                print('Looks like we have ourselves a genome paper')
                maybe = True

    for bigram in bigrams:
        if str.lower(bigram[0][0]) == 'complete' and str.lower(bigram[0][1]).endswith('genome'):
            print('We DEFINITELY have a genome paper')
            return True
    return maybe
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query')
    parser.add_argument('-w', '--wordfile')
    parser.add_argument('-d', dest='detect', action='store_true')
    parser.add_argument('-g', dest='genome', action='store_true')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    if args.wordfile is None:
        if args.verbose:
            print('No wordfile. Defaulting to human stuff')
        wordlist = ['GCRh37', 'GCRh38', 'human']
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
    doi = None
    meshterms = set()
    papers = []

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
            if 'MeshHeadingList' in paper["MedlineCitation"]:
                for mesh in paper["MedlineCitation"]['MeshHeadingList']:
                    for term in nltk.word_tokenize(mesh["DescriptorName"]):
                        meshterms.add(wnl.lemmatize(str.lower(term), 'n'))
                if args.verbose:
                    print("Found MeSH terms:\n{0}".format(meshterms))
            abstract = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
            doi = (elem for elem in paper["PubmedData"]["ArticleIdList"] if elem.attributes["IdType"] == 'doi')
            break

    if not abstract:
        print('No abstracts available')
    else:
        ## tokenise into word context
        tokens = nltk.word_tokenize(abstract)

        ## identify named entities
        entities = chunk_entities(tokens)
        if args.verbose:
            print('Chunks:\n{0}'.format(entities))
        
        # extract base phases
        base_phases = extract_phases(tokens, wordlist)
        filtered_ngrams = find_and_filter_bigrams(tokens)
        if args.verbose:
            print('Located potential base phases:\n{0}'.format(base_phases))

        if args.detect:
            if detectotext(filtered_ngrams, base_phases, args.query, wordlist) == False:
                ## inject mesh terms and extract phases
                tokens.extend(meshterms)
                all_phases = extract_phases(tokens, wordlist)
                mesh_bigrams = find_and_filter_bigrams(tokens)
                if args.verbose:
                    print('Located potential MeSH phases:\n{0}'.format(all_phases))
                detectotext(mesh_bigrams, all_phases, args.query, wordlist)
        
        if args.genome:
            detected = detectogenome(filtered_ngrams, base_phases, args.query)
            if detected:
                content = grab_pmid(id_list[0])
                seq_list = ['GenBank', 'ID', 'No.', 'no.', 'ENA', 'SRA', 'NCBI', 'genome', 'accession', 'annotation', 'DDBJ', 'EMBL', 'deposited']
                full_tokens = nltk.word_tokenize(content)
                no_stop_tokens = [w.translate(translator) for w in full_tokens if not w in stopwords.words('english')]
                entities = chunk_entities(full_tokens)
                full_phases = extract_phases(full_tokens, seq_list)
                full_filtered_ngrams = find_and_filter_trigrams(no_stop_tokens)
                if args.verbose:
                    print('Located full-text phases:\n{0}'.format(full_phases))
                    print(full_filtered_ngrams)

                refseq_matches = set(refseq_pattern.findall(content))
                insdc_matches = set(insdc_pattern.findall(content))
                genbank_matches = set(genbank_pattern.findall(content))

                fallback_matches = fallback_json.findall(content)
                if not refseq_matches and not insdc_matches and not genbank_matches:
                    for fbm in fallback_matches:
                        refseq_matches = set(refseq_pattern.findall(fbm))
                        insdc_matches = set(insdc_pattern.findall(fbm))
                        genbank_matches = set(genbank_pattern.findall(fbm))

                print("Accessions:")
                print("Refseq: {0}".format(['https://www.ncbi.nlm.nih.gov/refseq/?term={}'.format(m) for m in refseq_matches]))
                print("INSDC: {0}".format(['https://www.ncbi.nlm.nih.gov/sra/?term={}'.format(m) for m in insdc_matches]))
                print("GenBank: {0}".format(['https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(m) for m in genbank_matches]))

