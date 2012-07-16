from lxml import etree
from urllib2 import urlopen, URLError, Request
from urllib import urlencode
import re

class FetchGene:
    "Fetch gene details from the NCBI Entrez database"
    
    def __ie(self, item):
        "Check if an item exists in list, if not, return none"
        if len(item) > 0:
            return item[0]
        else:
            return None
    
    def fetchEntrez(self, entrez_id):
        "Fetch a gene its entrez gene id"
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&tool=DigitalAgeingAtlas&retmode=xml&id="
        try:
            dom = etree.parse(urlopen(url+str(entrez_id))).getroot()
        except URLError:
            print "error in URL: "+URLError.reason
            
        gene = dom.getchildren()[0]
        
        extracted = {}
        
        extracted['entrez_id'] = self.__ie(gene.xpath('Entrezgene_track-info/Gene-track/Gene-track_geneid/text()'))
        extracted['symbol'] = self.__ie(gene.xpath('Entrezgene_gene/Gene-ref/Gene-ref_locus/text()'))
        extracted['name'] = self.__ie(gene.xpath('Entrezgene_gene/Gene-ref/Gene-ref_desc/text()'))
        extracted['description'] = self.__ie(gene.xpath('Entrezgene_summary/text()'))#gene.xpath('Entrezgene_prot/Prot-ref/Prot-ref_desc/text()'))
        
        extracted['unigene'] = self.__ie(gene.xpath('Entrezgene_comments/Gene-commentary/Gene-commentary_heading[.="Additional Links"]/../Gene-commentary_comment/Gene-commentary/Gene-commentary_text[.="UniGene"]/../Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_str/text()'))
        extracted['uniprot'] = self.__ie(gene.xpath('Entrezgene_comments/Gene-commentary/Gene-commentary_heading[. = "NCBI Reference Sequences (RefSeq)"]/../Gene-commentary_comment/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_comment/Gene-commentary/Gene-commentary_heading[.="UniProtKB"]/../Gene-commentary_comment/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_db[.="UniProtKB/Swiss-Prot"]/../Dbtag_tag/Object-id/Object-id_str/text()'))
        extracted['omim'] = self.__ie(gene.xpath('Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_db[.="MIM"]/../Dbtag_tag/Object-id/Object-id_id/text()'))
        extracted['ensembl'] = self.__ie(gene.xpath('Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_db[.="Ensembl"]/../Dbtag_tag/Object-id/Object-id_str/text()'))

        extracted['species'] = self.__ie(gene.xpath('Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_common/text()'))
        
        alias_list = gene.xpath('Entrezgene_gene/Gene-ref/Gene-ref_syn')
        if len(alias_list) > 0:
            for alias in alias_list:
                if alias.text.strip() != '':
                    extracted['alias'] += ' '+alias.text
                    
        '''
        go_root = gene.xpath('Entrezgene_properties/Gene-commentary/Gene-commentary_heading[.="GeneOntology"]/../Gene-commentary_comment/Gene-commentary')
        if len(go_root) > 0:
            go = {}
            extracted['go_terms'] = []
            for cat in go_root:
                go_type = cat.xpath('Gene-commentary_label/text()')[0]
                terms = cat.xpath('Gene-commentary_comment/Gene-commentary')
                prev_term_id = ''
                for t in terms:
                    go = {}
                    go['type'] = go_type.lower()
                    go['go_id'] = t.xpath('Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_id/text()')[0]
                    go['go_term'] = t.xpath('Gene-commentary_source/Other-source/Other-source_anchor/text()')[0]
                    if prev_term_id != go['go_id']:
                        extracted['go_terms'].append(go)
                    prev_term_id = go['go_id']
                
        homologene_id = gene.xpath('Entrezgene_homology/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_id/text()')
        
        if len(homologene_id) > 0:
            hurl = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=homologene&retmode=xml&id="
            try:
                hdom = etree.parse(urlopen(hurl+homologene_id[0])).getroot()
            except URLError:
                print "error in URL: "+URLError.reason
            
            genes = hdom.xpath('HG-Entry/HG-Entry_genes/HG-Gene')
            hg = {}
            extracted['homologs'] = []
            for g in genes:
                hg['organism'] = g.xpath('HG-Gene_taxid/text()')
                if len(hg['organism']) > 0:
                    hg['organism'] = hg['organism'][0]
                hg['homolog_entrez_id'] = g.xpath('HG-Gene_geneid/text()')
                if len(hg['homolog_entrez_id']) > 0:
                    hg['homolog_entrez_id'] = hg['homolog_entrez_id'][0]
                hg['symbol'] = g.xpath('HG-Gene_symbol/text()')
                if len(hg['symbol']) > 0:
                    hg['symbol'] = hg['symbol'][0]
                extracted['homologs'].append(hg)
        '''
        
        return extracted
    
    def fetchSymbol(self, symbol, species):
        "Fetch a gene its symbol and species"
        pass
        
class FetchReference:
    "Fetch a reference from the NCBI Pubmed database"
    
    def __ie(self, item, makeNum=False):
        "Check if an item exists in list, if not, return none"
        if len(item) > 0:
            if makeNum:
                return re.sub("[^0-9]", "", item[0])
            else:
                return item[0]
        else:
            return None
    
    def fetchPubmed(self, pubmed):
        "Use a pubmed ID to fetch a reference from the PubMed database"
        
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&tool=DigitalAgeingAtlas&retmode=xml&id="
        try:
            dom = etree.parse(urlopen(url+str(pubmed))).getroot()
        except URLError as e:
            print "error in URL: "+e.reason
            dom = None
        if dom is not None:
            article = dom.xpath('/PubmedArticleSet/PubmedArticle/MedlineCitation/Article')[0]
            
            extracted = {}
            extracted['pubmed'] = pubmed
            extracted['title'] = self.__ie(article.xpath('ArticleTitle/text()'))
            extracted['volume'] = self.__ie(article.xpath('Journal/JournalIssue/Volume/text()'))
            extracted['journal'] = self.__ie(article.xpath('Journal/ISOAbbreviation/text()'))
            
            if self.__ie(article.xpath('Journal/JournalIssue/PubDate/Year/text()')) is not None:
                extracted['year'] = article.xpath('Journal/JournalIssue/PubDate/Year/text()')[0]
            else:
                extracted['year'] = self.__ie(article.xpath('Journal/JournalIssue/PubDate/MedlineDate/text()'), makeNum=True)
            
            pages = article.xpath('string(Pagination/MedlinePgn/text())')
            if '-' in pages:
                pages = pages.split('-')
                extracted['pages'] = pages[0]+'-'+pages[0][0:(len(pages[0]) - len(pages[1]))]+pages[1]
            else:
                extracted['pages'] = pages
                
            authors = article.xpath('AuthorList')[0]
            if len(authors) > 2:
                extracted['author'] = authors[0].xpath('LastName/text()')[0]+' et al'
            elif len(authors) == 2:
                extracted['author'] = authors[0].xpath('LastName/text()')[0]+' and '+authors[1].xpath('LastName/text()')[0]
            else:
                if authors[0].xpath('ForeName/text()') > 0:
                    firstName = authors[0].xpath('ForeName/text()')[0]
                elif authors[0].xpath('FirstName/text()') > 0:
                    firstName = authors[0].xpath('FirstName/text()')[0]
                else:
                    firstName = "";
                extracted['author'] = authors[0].xpath('LastName/text()')[0]+' '+firstName
                
            return extracted
        return {}

class FetchDetails:
    
    def entrezToUniProtID(self, entrez_id):
        "Take an EntrezGene ID and convert to a corrosponding UniProt ID"
        url = 'http://www.uniprot.org/mapping/'
        params = {
            'from': 'P_ENTREZGENEID',
            'to': 'ID',
            'format': 'list',
            'query': str(entrez_id),
        }
        data = urlencode(params)
        request = Request(url, data)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None
        
        rsplit = response.read().split(None, 1)
        return rsplit[0]

    def fetchDetailsFromUniProt(self, uniprotid):
        "Fetch details from the UniProt database"
        details = {}
        
        url = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch'
        params = {
            'db': 'uniprotkb',
            'format': 'uniprotxml',
            'style': 'raw',
            'pageHTML': 'false',
            'id': uniprotid,
        }
        data = urlencode(params)
        request = Request(url, data)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None

        dom = etree.parse(response).getroot()
        details['function'] = dom.xpath('string(u:entry/u:comment[@type="function"]/u:text/text())', namespaces={'u': 'http://uniprot.org/uniprot'})

        return details

    def fetchDetailsFromEPD(self, acc_num):
        "Fetch details from EPD using an SRS query"
        details = {}

        return details

    def fetchDetailsFromRefSeq(self, entrez_id):
        "Fetch details from RefSeq"
        details = {}
        
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"#?db=refseq&tool=DigitalAgeingAtlas&retmode=xml&id="
        params = {
            'db': 'refseq',
            'tool': 'DigitalAgeingAtlas',
            'retmode': 'xml',
            'id': entrez_id,
        }
        data = urlencode(params)
        print data
        request = Request(url, data)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None
        print response.read()
