from lxml import etree
from urllib2 import urlopen, URLError, Request
from urllib import urlencode
import re
import ConfigParser
import json
import pprint

""" DEPRECIATED: PLEASE USE FetchDetails CLASS """
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
    
    def fetchPubmed(self, pubmed, withTerms=False):
        "Use a pubmed ID to fetch a reference from the PubMed database"
        
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&tool=FetchScript&retmode=xml&id="
        try:
            dom = etree.parse(urlopen(url+str(pubmed))).getroot()
        except URLError as e:
            print "error in URL: "+e.reason
            dom = None
        if dom is not None:
            article = dom.xpath('/PubmedArticleSet/PubmedArticle/MedlineCitation/Article')[0]
            
            extracted = {}
            extracted['pubmed'] = pubmed
            extracted['title'] = article.xpath('string(ArticleTitle/text())')
            extracted['volume'] = article.xpath('string(Journal/JournalIssue/Volume/text())')
            extracted['issue'] = article.xpath('string(Journal/JournalIssue/Issue/text())')
            if article.xpath('string(Journal/ISOAbbreviation/text())') == '':
                extracted['journal'] = article.xpath('string(Journal/Title/text())')
            else:
                extracted['journal'] = article.xpath('string(Journal/ISOAbbreviation/text())')

            if article.xpath('string(Journal/JournalIssue/PubDate/Year/text())') != '':
                extracted['year'] = article.xpath('string(Journal/JournalIssue/PubDate/Year/text())')
            else:
                raw_date = article.xpath('string(Journal/JournalIssue/PubDate/MedlineDate/text())')
                if raw_date != '':
                    find_date = re.search(r'([0-9]{4})', raw_date)
                    if find_date:
                        extracted['year'] = find_date.group(1)
                    else:
                        extracted['year'] = self.__ie(article.xpath('Journal/JournalIssue/PubDate/MedlineDate/text()'), makeNum=True)
                else:
                    extracted['year'] = self.__ie(article.xpath('Journal/JournalIssue/PubDate/MedlineDate/text()'), makeNum=True)
            
            pages = article.xpath('string(Pagination/MedlinePgn/text())')
            if '-' in pages:
                pages = pages.split('-')
                extracted['pages'] = pages[0]+'-'+pages[0][0:(len(pages[0]) - len(pages[1]))]+pages[1]
            else:
                extracted['pages'] = pages
                
            authors = article.xpath('AuthorList')[0]

            firstName1stAuthor = authors[0].xpath('string(ForeName/text())') if authors[0].xpath('string(ForeName/text())') != '' else authors[0].xpath('string(FirstName/text())')
            lastName1stAuthor = authors[0].xpath('string(LastName/text())') if authors[0].xpath('string(LastName/text())') != '' else authors[0].xpath('string(LastName/text())')
            initials1stAuthor = authors[0].xpath('string(Initials/text())')
            if len(authors) > 2:
                extracted['author'] = lastName1stAuthor+' et al.' #+', '+initials1stAuthor+' et al.'
                extracted['author_initials'] = lastName1stAuthor+', '+initials1stAuthor+' et al.'
            elif len(authors) == 2:
                firstName2ndAuthor = authors[1].xpath('string(ForeName/text())') if authors[1].xpath('string(ForeName/text())') != '' else authors[1].xpath('string(FirstName/text())')
                lastName2ndAuthor = authors[1].xpath('string(LastName/text())') if authors[1].xpath('string(LastName/text())') != '' else authors[1].xpath('string(LastName/text())')
                initials2ndAuthor = authors[1].xpath('string(Initials/text())')
                extracted['author'] = lastName1stAuthor+' and '+lastName2ndAuthor #+', '+initials1stAuthor+' and '+lastName2ndAuthor+', '+initials2ndAuthor
                extracted['author_initials'] = lastName1stAuthor+', '+'.'.join(list(initials1stAuthor))+'. and '+lastName2ndAuthor+', '+'.'.join(list(initials2ndAuthor))+'.'
            else:
                extracted['author'] = lastName1stAuthor #+', '+initials1stAuthor+'.'
                extracted['author_initials'] = lastName1stAuthor+', '+initials1stAuthor+'.'

            extracted['authors'] = []
            for a in authors:
                extracted['authors'].append({
                    'first_name': a.xpath('string(ForeName/text())'),
                    'last_name': a.xpath('string(LastName/text())'),
                    'initials': a.xpath('string(Initials/text())')
                })

            pubTypes = article.xpath('PublicationTypeList/PublicationType')
            for ptype in pubTypes:
                pt = ptype.xpath('string(text())')
                if pt == 'Review':
                    extracted['review'] = 1

            if withTerms:
                extracted['terms'] = []
                meshTerms = dom.xpath('/PubmedArticleSet/PubmedArticle/MedlineCitation/MeshHeadingList/MeshHeading')
                for t in meshTerms:
                    currentTerm = ''
                    for s in t:
                        currentTerm += s.text+'/'
                    extracted['terms'].append(currentTerm.rstrip('/'))

            affiliation = dom.xpath('string(/PubmedArticleSet/PubmedArticle/MedlineCitation/Article/Affiliation/text())')
            addresses = re.findall(r'\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}\b', affiliation, flags=re.IGNORECASE)
            if len(addresses) > 0:
                address_list = []
                for a in addresses:
                    address_list.append(a)
                extracted['contact_addresses'] = ",".join(address_list)
            else:
                extracted['contact_addresses'] = ''

            return extracted
        return {}
    
class FetchDetails:

    def __init__(self, config_location='fetch.cfg'):
        config = ConfigParser.ConfigParser()
        config.read(config_location)

        self.biogrid_key = None
        if config.has_option('access_keys', 'biogrid'):
            self.biogrid_key = config.get('access_keys', 'biogrid')

    def translateID(self, id, translate_from, translate_to):
        allowed_translations = (
            'P_ENTREZGENEID', # Entrez gene ID
            'ID', # UniProt ID - !! can only translate TO !!
            'ACC+ID', # UniProt acc or ID - !! can only translate FROM !!
            'P_REFSEQ_AC', # RefSeq protein acc
            'REFSEQ_NT_ID', # RefSeq nucleotide acc
        )
        url = 'http://www.uniprot.org/mapping/'
        if translate_from in allowed_translations and translate_to in allowed_translations:
            params = {
                'from': translate_from, #'P_ENTREZGENEID',
                'to': translate_to, #'ACC+ID',
                'format': 'list',
                'query': str(id),
            }
            data = urlencode(params)
            request = Request(url, data)
            try:
                response = urlopen(request)
            except URLError as e:
                print "URL Error: {0}".format(e.reason)
                return None
            
            try:
                rsplit = response.read().split(None, 1)
                return rsplit[0]
            except:
                return None
        else:
            return None

    def convertToEntrezGeneID(self, identifier, to='symbol', organism='human'):
        "Fetch a gene id based on an identifier, e.g. a symbol or accession number"
        if to == 'accession':
            term = '{0}[accession] AND {1}[organism]'.format(identifier, organism)
        else:
            term = '{0}[gene name] AND {1}[organism]'.format(identifier, organism)

        details = {}
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'gene',
            'tool': 'DigitalAgeingAtlas',
            'term': term, 
        }
        data = urlencode(params)
        request = Request(url, data)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None

        dom = etree.parse(response).getroot()
        tid = dom.xpath('string(IdList/Id/text())')
        return tid

    def symbolToEntrezGeneID(self, symbol, organism='human'):
        self.convertToEntrezGeneID(symbol, organism=organism)
        """
        "Fetch a gene id based on a gene symbol"
        details = {}
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'gene',
            'tool': 'DigitalAgeingAtlas',
            'term': '{0}[gene name] AND {1}[organism]'.format(symbol, organism)
        }
        data = urlencode(params)
        request = Request(url, data)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None

        dom = etree.parse(response).getroot()
        tid = dom.xpath('string(IdList/Id/text())')
        return tid
        """

    def fetchDetailsFromUniProt(self, uniprotid):
        "Fetch details from the UniProt database"
        details = {
            'fetch_db_name': 'uniprot',
        }
        
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
        details['unigene'] = dom.xpath('string(u:entry/u:dbReference[@type="UniGene"][1]/@id)', namespaces={'u': 'http://uniprot.org/uniprot'})
        details['name'] = dom.xpath('string(u:entry/u:protein/u:recommendedName/u:fullName/text())', namespaces={'u': 'http://uniprot.org/uniprot'})
        details['symbol'] = dom.xpath('string(u:entry/u:gene/u:name/text())', namespaces={'u': 'http://uniprot.org/uniprot'})
        details['uniprot_accession'] = dom.xpath('string(u:entry/u:accession[1]/text())', namespaces={'u': 'http://uniprot.org/uniprot'})
        details['uniprot'] = dom.xpath('string(u:entry/u:name/text())', namespaces={'u': 'http://uniprot.org/uniprot'})

        details['interactions'] = []
        for i in dom.xpath('u:entry/u:comment[@type="interaction"]/u:interactant/u:label/text()', namespaces={'u': 'http://uniprot.org/uniprot'}):
            details['interactions'].append("".join(i))

        details['GO'] = []#{'function': [], 'component': [], 'process': []}
        for g in dom.xpath('u:entry/u:dbReference[@type="GO"]', namespaces={'u': 'http://uniprot.org/uniprot'}):
            go_id = g.xpath('string(@id)')
            go_term = g.xpath('string(u:property[@type="term"]/@value)', namespaces={'u': 'http://uniprot.org/uniprot'}).split(':', 1)
            details['GO'].append({'id': go_id, 'term': go_term[1], 'type': go_term[0]}) 

        details['Pfam'] = []
        for g in dom.xpath('u:entry/u:dbReference[@type="Pfam"]', namespaces={'u': 'http://uniprot.org/uniprot'}):
            pfam_id = g.xpath('string(@id)')
            entry_name = g.xpath('string(u:property[@type="entry name"]/@value)', namespaces={'u': 'http://uniprot.org/uniprot'})
            details['Pfam'].append({'id': pfam_id, 'entry_name': entry_name}) 

        return details

    def fetchDetailsFromEntrez(self, entrez_id):
        "Fetch details from NCBI entrez gene (Braving the mess that is NCBI XML)"
        details = {
            'fetch_db_name': 'entrez',
        }
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            'db': 'gene',
            'tool': 'DigitalAgeingAtlas',
            'retmode': 'xml',
            'id': entrez_id,
        }
        data = urlencode(params)
        request = Request(url, data)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None

        dom = etree.parse(response).getroot().getchildren()[0]
        #gene = dom.getchildren()[0]

        details['symbol'] = dom.xpath('string(Entrezgene_gene/Gene-ref/Gene-ref_locus/text())')
        # Sometimes symbol isn't actually there, it's here: 
        if details['symbol'] == '':
            details['symbol'] = dom.xpath('string(Entrezgene_gene/Gene-ref/Gene-ref_locus-tag/text())')
        details['name'] = dom.xpath('string(Entrezgene_gene/Gene-ref/Gene-ref_desc/text())')
        if details['name'] == '':
            details['name'] = dom.xpath('string(Entrezgene_prot/Prot-ref/Prot-ref_name/Prot-ref_name_E/text())')
        details['description'] = dom.xpath('string(Entrezgene_summary/text())')

        details['chromosome_location'] = dom.xpath('string(Entrezgene_location/Maps/Maps_display-str/text())')
        details['location_start'] = dom.xpath('string(Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_from/text())')
        details['location_end'] = dom.xpath('string(Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_to/text())')
        orientation = dom.xpath('string(Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_strand/Na-strand/@value)')
        if orientation == 'plus':
            details['orientation'] = 1
        else:
            details['orientation'] = -1

        details['accession'] = dom.xpath('string(Entrezgene_locus/Gene-commentary/Gene-commentary_label[.="RefSeqGene"]/../Gene-commentary_accession/text())')

        details['unigene'] = dom.xpath('string(Entrezgene_comments/Gene-commentary/Gene-commentary_heading[.="Additional Links"]/../Gene-commentary_comment/Gene-commentary/Gene-commentary_text[.="UniGene"]/../Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_str/text())')[3:]
        details['omim'] = dom.xpath('string(Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_db[.="MIM"]/../Dbtag_tag/Object-id/Object-id_id/text())')
        details['ensembl'] = dom.xpath('string(Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_db[.="Ensembl"]/../Dbtag_tag/Object-id/Object-id_str/text())')
        details['hprd'] = dom.xpath('string(Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_db[.="HPRD"]/../Dbtag_tag/Object-id/Object-id_str/text())')

        details['species'] = dom.xpath('string(Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_common/text())')

        details['homologene'] = dom.xpath('string(Entrezgene_homology/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_id/text())')
        
        aliases = []
        alias_list = dom.xpath('Entrezgene_gene/Gene-ref/Gene-ref_syn/Gene-ref_syn_E')
        print alias_list
        if len(alias_list) > 0:
            for alias in alias_list:
                if alias.text.strip() != '':
                    aliases.append(alias.text.strip())
        details['alias'] = ' '.join(aliases)

        details['homologene'] = dom.xpath('string(Entrezgene_homology/Gene-commentary/Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_id/text())')

        details['go'] = []
        go_cat = dom.xpath('Entrezgene_properties/Gene-commentary/Gene-commentary_heading[.="GeneOntology"]/../Gene-commentary_comment/Gene-commentary')
        for cat in go_cat:
            category = cat.xpath('string(Gene-commentary_label/text())')
            terms = cat.xpath('Gene-commentary_comment/Gene-commentary')
            for t in terms: 
                gid = t.xpath('string(Gene-commentary_source/Other-source/Other-source_src/Dbtag/Dbtag_tag/Object-id/Object-id_id/text())')
                gname = t.xpath('string(Gene-commentary_source/Other-source/Other-source_anchor/text())')
                details['go'].append({'go': gid, 'name': gname, 'type': category})

        return details
        
    def fetchDetailsFromNucleotide(self, nuc_id):
        "Fetch sequence details from NCBI Nucleotide database"
        details = {
            'fetch_db_name': 'nucleotide',
        }
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        for get_type in ('fasta_cds_na', 'fasta_cds_aa'):
            params = {
                'db': 'nuccore',
                'tool': 'DigitalAgeingAtlas',
                'rettype': get_type,
                'id': nuc_id,
            }
            data = urlencode(params)
            request = Request(url, data)
            try:
                response = urlopen(request)
            except URLError as e:
                print "URL Error: {0}".format(e.reason)
                return None
            fasta_list = response.read().split('\n')
            identifiers = re.search('(NM_.*)\.[0-9]_cdsid_(NP_[0-9].*)\.[0-9]\s', fasta_list[0])
            details['acc_orf'] = identifiers.group(1)
            details['acc_cds'] = identifiers.group(2)
            if get_type == 'fasta_cds_na':
                details['seq_orf'] = "".join(fasta_list[1:-1])
            else:
                details['seq_cds'] = "".join(fasta_list[1:-1])

        return details

    def fetchDetailsFromHPRD(self, symbol):
        "Fetch details from HPRD using API on Alfred server"
        """
        THIS HAS BEEN DEPRECIATED AND WILL BE REMOVED. HPRD IS NO LONGER BEING UPDATED.
        """
        details = {
            'fetch_db_name': 'HPRD',
        }
        url = "http://alfred.liv.ac.uk/api/hprd/"
        params = {
            'symbol': symbol,
        }
        data = urlencode(params)
        request = Request(url, data)

        response = urlopen(request)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None
        return response.read()

    def fetchDetailsFromdbSNP(self, identifier):
        """Fetch details from NCBI dbSNP""" 
        details = {
            'fetch_db_name': 'dbSNP',
        }
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            'db': 'snp',
            'tool': 'DigitalAgeingAtlas',
            'retmode': 'xml',
            'id': identifier,
        }
        data = urlencode(params)
        request = Request(url, data)
        try:
            response = urlopen(request)
        except URLError as e:
            print "URL Error: {0}".format(e.reason)
            return None

        dom = etree.parse(response).getroot()#.getchildren()[0]

        details = {}

        fxn = dom.xpath('s:Rs/s:Assembly[@reference="true"]/s:Component/s:MapLoc/s:FxnSet', namespaces={'s': 'http://www.ncbi.nlm.nih.gov/SNP/docsum'})
        gene_id = None
        if len(fxn) > 0:
            gene_id = fxn[0].attrib['geneId']
        if gene_id is None:
            ps = dom.xpath('string(s:Rs/s:PrimarySequence/@accession)', namespaces={'s': 'http://www.ncbi.nlm.nih.gov/SNP/docsum'})
            if ps != '':
                gene_id = self.convertToEntrezGeneID(ps, to='accession')

        details['entrez_id'] = gene_id 

        details['from'] = dom.xpath('string(s:Rs/s:Assembly[@reference="true"]/s:Component/s:MapLoc/@asnFrom)', namespaces={'s': 'http://www.ncbi.nlm.nih.gov/SNP/docsum'})

        return details

    def fetchDetailsFromBioGrid(self, identifier, tax_id):
        details = {
            'fetch_db_name': 'BioGrid',
        }
        url = 'http://webservice.thebiogrid.org/interactions'
        if self.biogrid_key is not None:
            params = {
                'searchNames': True,
                'geneList': identifier,
                'includeInteractors': True,
                'includeInteractorInteractions': False,
                'additionalIdentifierTypes': 'ENTREZ_GENE',
                'taxId': tax_id,
                'accesskey': self.biogrid_key,
                'format': 'json',
            }
            data = urlencode(params)
            request = Request(url, data)
            try:
                response = urlopen(request)
            except URLError as e:
                print "URL Error: {0}".format(e.reason)
                return None
            details = json.loads(response.read())
        else:
            return None

        return details
