import re, argparse, sys, scholarly, urllib
from Bio import Entrez

email = 'your.human.email@you.com'
r = re.compile(r"GRCh3[0-9]|human|[Hh]omo sapiens")

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    results = search(args.query)
    id_list = results['IdList']
    match = None

    if not id_list:
        # try arxiv
        # url = 'http://export.arxiv.org/api/query?search_query=ti:%22{0}%22:&start=0&max_results=1'.format(args.query)
        # data = urllib.request.urlopen(url).read()
        # print(data)

        # try scholar
        search_query = scholarly.search_pubs_query(args.query)
        for scholar_result in search_query:
            match = r.search(scholar_result.bib["abstract"])
    else:
        papers = fetch_details(id_list)
        for i, paper in enumerate(papers["PubmedArticle"]):
            abstract = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
            match = r.search(abstract)
    if match:
        print(args.query.rstrip('.') + ", in humans.")
    else:
        print("No humans detected.")
        

