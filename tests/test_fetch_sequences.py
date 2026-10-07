#!/usr/bin/env python

# Test fetch_sequences module
from teemi.design.fetch_sequences import *
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

def test_retrieve_sequences_from_ncbi():
    # random acc_number
    acc_numbers = ['Q05001']
    # call the function: 
    assert retrieve_sequences_from_ncbi(acc_numbers, '../tests/files_for_testing/test_fetch.fasta') == None

    # check if it is the correct sequences
    Q05001_seq = SeqIO.read('../teemi/tests/files_for_testing/test_fetch.fasta', 'fasta')
    assert str(Q05001_seq.seq[:10]) == 'MDSSSEKLSP'
    assert Q05001_seq.id == 'sp|Q05001.1|NCPR_CATRO'
 
def test_read_fasta_files(): 
    sequences = read_fasta_files('../teemi/tests/files_for_testing/test_fetch.fasta')
    assert sequences[0].seq[:10] == 'MDSSSEKLSP'
    assert sequences[0].id == 'sp|Q05001.1|NCPR_CATRO'

def test_retrieve_sequences_from_PDB(): 
    acc_numbers = ['Q1PQK4']
    sequences = retrieve_sequences_from_PDB(acc_numbers)
    assert str(sequences[0][0].seq[:10]) == 'MQSTTSVKLS'
    assert str(sequences[0][0].id) == 'sp|A0A2U1LIM9|NCPR1_ARTAN'


def test_read_genbank_files():
    test_gb = read_genbank_files('../teemi/tests/files_for_testing/MIA-HA-1.gb')
    assert str(test_gb[0].seq[20:30]) == 'AGTTATATAG'
    assert test_gb[0].name == 'MIA-HA-1'


def test_concatenate_genbank_records():
    feature_a = SeqFeature(FeatureLocation(5, 10, strand=1), type='CDS')
    feature_a.qualifiers['locus_tag'] = ['GENE_A']
    feature_b = SeqFeature(FeatureLocation(2, 6, strand=-1), type='CDS')
    feature_b.qualifiers['locus_tag'] = ['GENE_B']

    record_a = SeqRecord(Seq('A' * 20), id='chrom1', name='chrom1', features=[feature_a])
    record_b = SeqRecord(Seq('T' * 15), id='chrom2', name='chrom2', features=[feature_b])

    combined_record = concatenate_genbank_records(
        [record_a, record_b],
        record_id='aspergillus_combined',
        record_name='aspergillus_combined',
    )

    assert len(combined_record.seq) == 35
    assert combined_record.id == 'aspergillus_combined'
    assert combined_record.name == 'aspergillus_combined'
    assert len(combined_record.features) == 2
    assert int(combined_record.features[0].location.start) == 5
    assert int(combined_record.features[0].location.end) == 10
    assert int(combined_record.features[1].location.start) == 22
    assert int(combined_record.features[1].location.end) == 26
    assert combined_record.features[1].qualifiers['locus_tag'] == ['GENE_B']


## TODO: I cannot make these work with github actions yet. It lacks some packages(intermine). Have to figure this out later. 
# def test_fetch_promoter(): 
#     cyc1 = fetch_promoter('CYC1')
#     assert cyc1[:20] == 'GAGGCACCAGCGTCAGCATT'

# def test_fetch_multiple_promoters(): 
#     list_of_promoters = ['YAR035C-A', 'YGR067C']
#     seqs = fetch_multiple_promoters(list_of_promoters)

#     assert str(seqs[0].seq[:10]) == 'CCCTGGTGGC'
#     assert str(seqs[1].seq[:10])== 'AGACAACCTA'
#     assert seqs[0].id == 'YAR035C-A'
#     assert seqs[1].id == 'YGR067C'


# ---------------------------------------------------------------------------
# Additional tests - all network access is mocked
# ---------------------------------------------------------------------------
import sys
import types

import teemi.design.fetch_sequences as fetch_sequences


class _FakeEntrezHandle:
    def __init__(self, text):
        self._text = text

    def read(self):
        return self._text


def test_retrieve_sequences_from_ncbi_writes_fasta(tmp_path, monkeypatch):
    fasta_records = {
        'ACC_1': '>ACC_1 first protein\nMDSSSEKLSP\n',
        'ACC_2': '>ACC_2 second protein\nMQSTTSVKLS\n',
    }
    calls = []

    def fake_efetch(db, id, rettype, retmode):
        calls.append({'db': db, 'id': id, 'rettype': rettype, 'retmode': retmode})
        return _FakeEntrezHandle(fasta_records[id])

    monkeypatch.setattr(fetch_sequences.Entrez, 'efetch', fake_efetch)
    out_file = tmp_path / 'ncbi_hits.fasta'

    assert retrieve_sequences_from_ncbi(['ACC_1', 'ACC_2'], str(out_file)) is None

    # one efetch call per accession number, all of them from the protein db
    assert [call['id'] for call in calls] == ['ACC_1', 'ACC_2']
    assert {call['db'] for call in calls} == {'protein'}
    assert {call['rettype'] for call in calls} == {'fasta'}
    assert {call['retmode'] for call in calls} == {'text'}
    assert fetch_sequences.Entrez.email == 'youremail@gmail.com'

    # the responses are concatenated into one fasta file
    assert out_file.read_text() == fasta_records['ACC_1'] + fasta_records['ACC_2']
    written = read_fasta_files(str(out_file))
    assert [record.id for record in written] == ['ACC_1', 'ACC_2']
    assert str(written[0].seq) == 'MDSSSEKLSP'
    assert str(written[1].seq) == 'MQSTTSVKLS'


def test_retrieve_sequences_from_ncbi_db_argument(tmp_path, monkeypatch):
    used_db = []

    def fake_efetch(db, id, rettype, retmode):
        used_db.append(db)
        return _FakeEntrezHandle('>%s\nACGT\n' % id)

    monkeypatch.setattr(fetch_sequences.Entrez, 'efetch', fake_efetch)
    out_file = tmp_path / 'nuc.fasta'

    retrieve_sequences_from_ncbi(['X1'], str(out_file), db='nucleotide')

    assert used_db == ['nucleotide']
    assert out_file.read_text() == '>X1\nACGT\n'


def test_retrieve_sequences_from_ncbi_handles_failure(tmp_path, monkeypatch, capsys):
    def fake_efetch(**kwargs):
        raise OSError('no connection')

    monkeypatch.setattr(fetch_sequences.Entrez, 'efetch', fake_efetch)
    out_file = tmp_path / 'broken.fasta'

    assert retrieve_sequences_from_ncbi(['BAD_ACC'], str(out_file)) is None
    assert 'An exception occurred' in capsys.readouterr().out
    assert out_file.read_text() == ''


def test_retrieve_sequences_from_PDB(monkeypatch):
    fasta = {
        'Q1PQK4': '>sp|Q1PQK4|TEST_PROT test protein\nMQSTTSVKLS\n',
        'Q05001': '>sp|Q05001|OTHER_PROT other protein\nMDSSSEKLSP\n',
    }
    requested_urls = []

    class _FakeResponse:
        def __init__(self, text):
            self.text = text

    def fake_post(url):
        requested_urls.append(url)
        accession = url.rsplit('/', 1)[-1].replace('.fasta', '')
        return _FakeResponse(fasta[accession])

    monkeypatch.setattr(fetch_sequences.r, 'post', fake_post)

    sequences = retrieve_sequences_from_PDB(['Q1PQK4', 'Q05001'])

    assert requested_urls == [
        'http://www.uniprot.org/uniprot/Q1PQK4.fasta',
        'http://www.uniprot.org/uniprot/Q05001.fasta',
    ]
    assert len(sequences) == 2
    assert str(sequences[0][0].seq) == 'MQSTTSVKLS'
    assert sequences[0][0].id == 'sp|Q1PQK4|TEST_PROT'
    assert str(sequences[1][0].seq) == 'MDSSSEKLSP'


class _FakeQuery:
    def __init__(self, rows):
        self._rows = rows
        self.views = []
        self.constraints = []

    def add_view(self, *views):
        self.views.extend(views)

    def add_constraint(self, *args, **kwargs):
        self.constraints.append((args, kwargs))

    def rows(self):
        return self._rows


class _FakeService:
    def __init__(self, url):
        self.url = url
        self.query = None

    def new_query(self, root):
        self.root = root
        self.query = _FakeQuery(self.rows)
        return self.query


def _install_fake_intermine(monkeypatch, rows):
    """Replace the (python2 only) intermine package with a stub."""
    created = {}

    class Service(_FakeService):
        rows = None

        def __init__(self, url):
            _FakeService.__init__(self, url)
            self.rows = rows
            created['service'] = self

    package = types.ModuleType('intermine')
    webservice = types.ModuleType('intermine.webservice')
    webservice.Service = Service
    package.webservice = webservice
    monkeypatch.setitem(sys.modules, 'intermine', package)
    monkeypatch.setitem(sys.modules, 'intermine.webservice', webservice)
    return created


def test_fetch_promoter(monkeypatch):
    rows = [
        {'flankingRegions.sequence.residues': 'AAAACCCCGGGGTTTT'},
        {'flankingRegions.sequence.residues': 'GAGGCACCAGCGTCAG'},
    ]
    created = _install_fake_intermine(monkeypatch, rows)

    promoter = fetch_promoter('CYC1')

    # the residues of the last returned row are used
    assert promoter == 'GAGGCACCAGCGTCAG'

    service = created['service']
    assert service.url == 'https://yeastmine.yeastgenome.org/yeastmine/service'
    assert service.root == 'Gene'
    assert service.query.views == [
        'secondaryIdentifier',
        'symbol',
        'length',
        'flankingRegions.direction',
        'flankingRegions.sequence.length',
        'flankingRegions.sequence.residues',
    ]
    assert service.query.constraints == [
        (('Gene', 'LOOKUP', 'CYC1', 'S. cerevisiae'), {'code': 'B'}),
        (('flankingRegions.direction', '=', 'upstream'), {'code': 'C'}),
        (('flankingRegions.distance', '=', '1.0kb'), {'code': 'A'}),
        (('flankingRegions.includeGene', '=', 'false'), {'code': 'D'}),
    ]


def test_fetch_promoter_without_hits(monkeypatch):
    _install_fake_intermine(monkeypatch, [])

    assert fetch_promoter('NOT_A_GENE') == ''


def test_fetch_multiple_promoters(monkeypatch):
    promoters = {'YAR035C-A': 'CCCTGGTGGC', 'YGR067C': 'AGACAACCTA'}
    monkeypatch.setattr(
        fetch_sequences, 'fetch_promoter', lambda name: promoters[name]
    )

    records = fetch_multiple_promoters(['YAR035C-A', 'YGR067C'])

    assert len(records) == 2
    assert str(records[0].seq) == 'CCCTGGTGGC'
    assert str(records[1].seq) == 'AGACAACCTA'
    assert records[0].id == 'YAR035C-A'
    assert records[1].id == 'YGR067C'
    assert records[0].name == 'YAR035C-A Promoter'
    assert records[1].name == 'YGR067C Promoter'
    assert records[0].description == (
        'Defined as being 1kb upstream of the TSS and fetched through Intermines API'
    )
