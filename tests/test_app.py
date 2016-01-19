import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from pybasespace.app import *


def test_json_depcopy():

    assert(json_deepcopy(metadatatemplate) == json_deepcopy(metadatatemplate))

    assert(json_deepcopy(metadatatemplate) is not json_deepcopy(metadatatemplate))


def test_parse_appsessionparams():

    testdir = os.path.dirname(os.path.abspath(__file__))

    appsession_location = testdir + '/data/input/AppSession.json'

    appsessionhref, appsessionparams = read_appsession(appsession_location)

    assert(appsessionhref == u'v1pre3/appsessions/31951397')

    assert(appsessionparams ==
 [{u'Content': u'False',
   u'Description': u'',
   u'Href': u'v1pre3/appsessions/31951397/properties/BaseSpace.Private.IsMultiNode',
   u'Name': u'BaseSpace.Private.IsMultiNode',
   u'Type': u'string'},
  {u'Content': u'Crispr-Eating 01/17/2016 10:21:02',
   u'Description': u'as analysis name',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.app-session-name',
   u'Name': u'Input.app-session-name',
   u'Type': u'string'},
  {u'Content': u'25',
   u'Description': u'Binding Interference Spacing',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.binding_interference_spacing',
   u'Name': u'Input.binding_interference_spacing',
   u'Type': u'string'},
  {u'Content': u'25',
   u'Description': u'Blast Chunk Size',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.blast_chunk_size',
   u'Name': u'Input.blast_chunk_size',
   u'Type': u'string'},
  {u'Content': u'1000',
   u'Description': u'Blast Max HSPs',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.blast_max_hsps',
   u'Name': u'Input.blast_max_hsps',
   u'Type': u'string'},
  {u'Content': u'Mouse8',
   u'Description': u'',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.genome_id',
   u'Name': u'Input.genome_id',
   u'Type': u'string'},
  {u'Content': u'chr6:136640001-136680000',
   u'Description': u'',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.genomic_coord',
   u'Name': u'Input.genomic_coord',
   u'Type': u'string'},
  {u'Description': u'',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.Projects',
   u'HrefItems': u'v1pre3/appsessions/31951397/properties/Input.Projects/items',
   u'Items': [{u'DateCreated': u'2016-01-14T22:02:43.0000000',
     u'DateModified': u'2016-01-14T22:02:43.0000000',
     u'Description': u'crispr-eating buffet project 1',
     u'HasCollaborators': False,
     u'Href': u'v1pre3/projects/27855842',
     u'Id': u'27855842',
     u'Name': u'crispr-eating buffet project 1',
     u'TotalSize': 0,
     u'UserOwnedBy': {u'GravatarUrl': u'https://secure.gravatar.com/avatar/c1de6e99641cad724d70888ea039b087.jpg?s=20&d=https%3a%2f%2fbasespace.illumina.com%2fpublic%2fimages%2fDefaultCustomerGravatar.png&r=PG',
      u'Href': u'v1pre3/users/7036038',
      u'Id': u'7036038',
      u'Name': u'Alain Domissy'}}],
   u'ItemsDisplayedCount': 1,
   u'ItemsTotalCount': 1,
   u'Name': u'Input.Projects',
   u'Type': u'project[]'},
  {u'Content': {u'DateCreated': u'2016-01-14T22:02:43.0000000',
    u'DateModified': u'2016-01-14T22:02:43.0000000',
    u'Description': u'crispr-eating buffet project 1',
    u'HasCollaborators': False,
    u'Href': u'v1pre3/projects/27855842',
    u'Id': u'27855842',
    u'Name': u'crispr-eating buffet project 1',
    u'TotalSize': 0,
    u'UserOwnedBy': {u'GravatarUrl': u'https://secure.gravatar.com/avatar/c1de6e99641cad724d70888ea039b087.jpg?s=20&d=https%3a%2f%2fbasespace.illumina.com%2fpublic%2fimages%2fDefaultCustomerGravatar.png&r=PG',
     u'Href': u'v1pre3/users/7036038',
     u'Id': u'7036038',
     u'Name': u'Alain Domissy'}},
   u'Description': u'under project',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.project_id',
   u'Name': u'Input.project_id',
   u'Type': u'project'},
  {u'Description': u'under project Attributes',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.project_id.attributes',
   u'HrefItems': u'v1pre3/appsessions/31951397/properties/Input.project_id.attributes/items',
   u'Items': [[{u'Key': u'FieldId', u'Values': [u'project_id']},
     {u'Key': u'ResourceType', u'Values': [u'project']},
     {u'Key': u'ResourceId', u'Values': [u'27855842']},
     {u'Key': u'ResourceHref', u'Values': [u'v1pre3/projects/27855842']}]],
   u'ItemsDisplayedCount': 1,
   u'ItemsTotalCount': 1,
   u'Name': u'Input.project_id.attributes',
   u'Type': u'map[]'},
  {u'Description': u'Restriction Enzymes',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.restriction_enzymes',
   u'HrefItems': u'v1pre3/appsessions/31951397/properties/Input.restriction_enzymes/items',
   u'Items': [u'BfaI', u'ScrFI', u'HpaII'],
   u'ItemsDisplayedCount': 3,
   u'ItemsTotalCount': 3,
   u'Name': u'Input.restriction_enzymes',
   u'Type': u'string[]'},
  {u'Description': u'',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.sample-id',
   u'HrefItems': u'v1pre3/appsessions/31951397/properties/Input.sample-id/items',
   u'Items': [{u'DateCreated': u'2015-11-08T06:11:08.0000000',
     u'ExperimentName': u'SRA Import',
     u'Genome': {u'Build': u'hg19',
      u'DisplayName': u'Homo Sapiens - UCSC (hg19)',
      u'Href': u'v1pre3/genomes/4',
      u'HrefAnnotations': u'v1pre3/genomes/4/annotations',
      u'Id': u'4',
      u'Source': u'UCSC',
      u'SpeciesName': u'Homo sapiens'},
     u'Href': u'v1pre3/samples/30427461',
     u'Id': u'30427461',
     u'IsMerged': False,
     u'IsPairedEnd': False,
     u'Name': u'SRR1015697',
     u'NumReadsPF': 81725433,
     u'NumReadsRaw': 81725433,
     u'Read1': 50,
     u'Read2': 0,
     u'SampleId': u'SRR1015697_GSM1249136-RTT5-iPS-15bi-Homo-sapiens-RNA-Seq',
     u'Status': u'Complete',
     u'StatusSummary': u'Application completed successfully',
     u'TotalClustersPF': 81725433,
     u'TotalClustersRaw': 81725433,
     u'TotalReadsPF': 81725433,
     u'TotalReadsRaw': 81725433,
     u'TotalSize': 3663526295,
     u'UserOwnedBy': {u'GravatarUrl': u'https://secure.gravatar.com/avatar/c1de6e99641cad724d70888ea039b087.jpg?s=20&d=https%3a%2f%2fbasespace.illumina.com%2fpublic%2fimages%2fDefaultCustomerGravatar.png&r=PG',
      u'Href': u'v1pre3/users/7036038',
      u'Id': u'7036038',
      u'Name': u'Alain Domissy'}},
    {u'DateCreated': u'2015-11-08T06:10:26.0000000',
     u'ExperimentName': u'SRA Import',
     u'Genome': {u'Build': u'hg19',
      u'DisplayName': u'Homo Sapiens - UCSC (hg19)',
      u'Href': u'v1pre3/genomes/4',
      u'HrefAnnotations': u'v1pre3/genomes/4/annotations',
      u'Id': u'4',
      u'Source': u'UCSC',
      u'SpeciesName': u'Homo sapiens'},
     u'Href': u'v1pre3/samples/30436467',
     u'Id': u'30436467',
     u'IsMerged': False,
     u'IsPairedEnd': False,
     u'Name': u'SRR1015696',
     u'NumReadsPF': 98343780,
     u'NumReadsRaw': 98343780,
     u'Read1': 50,
     u'Read2': 0,
     u'SampleId': u'SRR1015696_GSM1249135-RTT1-iPS-10bi-Homo-sapiens-RNA-Seq',
     u'Status': u'Complete',
     u'StatusSummary': u'Application completed successfully',
     u'TotalClustersPF': 98343780,
     u'TotalClustersRaw': 98343780,
     u'TotalReadsPF': 98343780,
     u'TotalReadsRaw': 98343780,
     u'TotalSize': 4606425168,
     u'UserOwnedBy': {u'GravatarUrl': u'https://secure.gravatar.com/avatar/c1de6e99641cad724d70888ea039b087.jpg?s=20&d=https%3a%2f%2fbasespace.illumina.com%2fpublic%2fimages%2fDefaultCustomerGravatar.png&r=PG',
      u'Href': u'v1pre3/users/7036038',
      u'Id': u'7036038',
      u'Name': u'Alain Domissy'}}],
   u'ItemsDisplayedCount': 2,
   u'ItemsTotalCount': 2,
   u'Name': u'Input.sample-id',
   u'Type': u'sample[]'},
  {u'Description': u' Attributes',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.sample-id.attributes',
   u'HrefItems': u'v1pre3/appsessions/31951397/properties/Input.sample-id.attributes/items',
   u'Items': [[{u'Key': u'FieldId', u'Values': [u'sample-id']},
     {u'Key': u'ResourceType', u'Values': [u'sample']},
     {u'Key': u'ResourceId', u'Values': [u'30427461']},
     {u'Key': u'ResourceHref', u'Values': [u'v1pre3/samples/30427461']}],
    [{u'Key': u'FieldId', u'Values': [u'sample-id']},
     {u'Key': u'ResourceType', u'Values': [u'sample']},
     {u'Key': u'ResourceId', u'Values': [u'30436467']},
     {u'Key': u'ResourceHref', u'Values': [u'v1pre3/samples/30436467']}]],
   u'ItemsDisplayedCount': 2,
   u'ItemsTotalCount': 2,
   u'Name': u'Input.sample-id.attributes',
   u'Type': u'map[]'},
  {u'Description': u'',
   u'Href': u'v1pre3/appsessions/31951397/properties/Input.Samples',
   u'HrefItems': u'v1pre3/appsessions/31951397/properties/Input.Samples/items',
   u'Items': [{u'DateCreated': u'2015-11-08T06:11:08.0000000',
     u'ExperimentName': u'SRA Import',
     u'Genome': {u'Build': u'hg19',
      u'DisplayName': u'Homo Sapiens - UCSC (hg19)',
      u'Href': u'v1pre3/genomes/4',
      u'HrefAnnotations': u'v1pre3/genomes/4/annotations',
      u'Id': u'4',
      u'Source': u'UCSC',
      u'SpeciesName': u'Homo sapiens'},
     u'Href': u'v1pre3/samples/30427461',
     u'Id': u'30427461',
     u'IsMerged': False,
     u'IsPairedEnd': False,
     u'Name': u'SRR1015697',
     u'NumReadsPF': 81725433,
     u'NumReadsRaw': 81725433,
     u'Read1': 50,
     u'Read2': 0,
     u'SampleId': u'SRR1015697_GSM1249136-RTT5-iPS-15bi-Homo-sapiens-RNA-Seq',
     u'Status': u'Complete',
     u'StatusSummary': u'Application completed successfully',
     u'TotalClustersPF': 81725433,
     u'TotalClustersRaw': 81725433,
     u'TotalReadsPF': 81725433,
     u'TotalReadsRaw': 81725433,
     u'TotalSize': 3663526295,
     u'UserOwnedBy': {u'GravatarUrl': u'https://secure.gravatar.com/avatar/c1de6e99641cad724d70888ea039b087.jpg?s=20&d=https%3a%2f%2fbasespace.illumina.com%2fpublic%2fimages%2fDefaultCustomerGravatar.png&r=PG',
      u'Href': u'v1pre3/users/7036038',
      u'Id': u'7036038',
      u'Name': u'Alain Domissy'}},
    {u'DateCreated': u'2015-11-08T06:10:26.0000000',
     u'ExperimentName': u'SRA Import',
     u'Genome': {u'Build': u'hg19',
      u'DisplayName': u'Homo Sapiens - UCSC (hg19)',
      u'Href': u'v1pre3/genomes/4',
      u'HrefAnnotations': u'v1pre3/genomes/4/annotations',
      u'Id': u'4',
      u'Source': u'UCSC',
      u'SpeciesName': u'Homo sapiens'},
     u'Href': u'v1pre3/samples/30436467',
     u'Id': u'30436467',
     u'IsMerged': False,
     u'IsPairedEnd': False,
     u'Name': u'SRR1015696',
     u'NumReadsPF': 98343780,
     u'NumReadsRaw': 98343780,
     u'Read1': 50,
     u'Read2': 0,
     u'SampleId': u'SRR1015696_GSM1249135-RTT1-iPS-10bi-Homo-sapiens-RNA-Seq',
     u'Status': u'Complete',
     u'StatusSummary': u'Application completed successfully',
     u'TotalClustersPF': 98343780,
     u'TotalClustersRaw': 98343780,
     u'TotalReadsPF': 98343780,
     u'TotalReadsRaw': 98343780,
     u'TotalSize': 4606425168,
     u'UserOwnedBy': {u'GravatarUrl': u'https://secure.gravatar.com/avatar/c1de6e99641cad724d70888ea039b087.jpg?s=20&d=https%3a%2f%2fbasespace.illumina.com%2fpublic%2fimages%2fDefaultCustomerGravatar.png&r=PG',
      u'Href': u'v1pre3/users/7036038',
      u'Id': u'7036038',
      u'Name': u'Alain Domissy'}}],
   u'ItemsDisplayedCount': 2,
   u'ItemsTotalCount': 2,
   u'Name': u'Input.Samples',
   u'Type': u'sample[]'},
  {u'Description': u'',
   u'Href': u'v1pre3/appsessions/31951397/properties/Output.Projects',
   u'HrefItems': u'v1pre3/appsessions/31951397/properties/Output.Projects/items',
   u'Items': [{u'DateCreated': u'2016-01-14T22:02:43.0000000',
     u'DateModified': u'2016-01-14T22:02:43.0000000',
     u'Description': u'crispr-eating buffet project 1',
     u'HasCollaborators': False,
     u'Href': u'v1pre3/projects/27855842',
     u'Id': u'27855842',
     u'Name': u'crispr-eating buffet project 1',
     u'TotalSize': 0,
     u'UserOwnedBy': {u'GravatarUrl': u'https://secure.gravatar.com/avatar/c1de6e99641cad724d70888ea039b087.jpg?s=20&d=https%3a%2f%2fbasespace.illumina.com%2fpublic%2fimages%2fDefaultCustomerGravatar.png&r=PG',
      u'Href': u'v1pre3/users/7036038',
      u'Id': u'7036038',
      u'Name': u'Alain Domissy'}}],
   u'ItemsDisplayedCount': 1,
   u'ItemsTotalCount': 1,
   u'Name': u'Output.Projects',
   u'Type': u'project[]'}]
           )

    assert(parse_appsessionparams(appsessionparams) == {
                u'input.binding_interference_spacing': u'25',
                u'input.blast_chunk_size': u'25',
                u'input.blast_max_hsps': u'1000',
                u'input.genome_id': u'Mouse8',
                u'input.genomic_coord': u'chr6:136640001-136680000',
                'input.projects_ids': u'27855842',
                u'input.restriction_enzymes': [u'BfaI', u'ScrFI', u'HpaII'],
                'input.samples': [
                    {'href': u'v1pre3/samples/30427461',
                        'id': u'30427461',
                        'name': u'SRR1015697'},
                    {'href': u'v1pre3/samples/30436467',
                        'id': u'30436467',
                        'name': u'SRR1015696'}],
                'output.projects_ids': u'27855842',
                'input.project_id': u'27855842'
                }
           )