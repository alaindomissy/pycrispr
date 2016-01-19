import os
import json

from pybasespace.payload import payload

appsession_location = os.environ.get('APPSESSION', '/data/input/AppSession.json')


def read_appsession(appsession_location):
    with open(appsession_location) as hdl:
        appsession = json.load(hdl)
    appsessionparams = appsession['Properties']['Items']
    appsessionhref = appsession['Href']
    return appsessionhref, appsessionparams


def parse_appsessionparams(appsessionparams):
    arguments_with_content = ['input.binding_interference_spacing',
                              'input.blast_chunk_size',
                              'input.blast_max_hsps',
                              'input.genome_id',
                              'input.genomic_coord']
    param_values = {}
    param_values.update(
        {param.get('Name').lower(): param.get('Content')
         for param in appsessionparams
         if param.get('Name').lower() in arguments_with_content
        }
    )
    param_values.update(
        {param.get('Name').lower(): param.get('Items')
         for param in appsessionparams
         if param.get('Name').lower()=='input.restriction_enzymes'
        }
    )
    param_values.update(
        {'input.project_id': param.get('Content').get('Id')
         for param in appsessionparams
         if param.get('Name') == 'Input.project_id'
        }
    )
    param_values.update(
        {'input.projects_ids': project.get('Id')
         for param in appsessionparams
         if param.get('Name')=='Input.Projects'
         for project in param['Items']
        }
    )
    param_values.update(
        {'output.projects_ids': project.get('Id')
         for param in appsessionparams
         if param.get('Name')=='Output.Projects'
         for project in param.get('Items')
        }
    )
    param_values.update(
        {'input.samples':
            [{'id': sample['Id'], 'href':sample['Href'], 'name': sample['Name']}
             for param in appsessionparams
             if param.get('Name')=='Input.Samples'
             for sample in param.get('Items')
            ]
        }
    )
    return param_values


metadatatemplate = {
     "Properties": [{"Items": [], "Type": "sample[]", "Name": "Input.Samples"}],
     "Name": "",
     "HrefAppSession": "",
     "Description": ""
}


def json_deepcopy(obj):
    return json.loads(json.dumps(obj))


def write_metadata(name, description, appsessionhref, sampleshrefs, output_dir):
    metadata = json_deepcopy(metadatatemplate)
    metadata['Name'] = name
    metadata['Description'] = description
    metadata['HrefAppSession'] = appsessionhref
    metadata['Properties'][0]['Items'].extend(sampleshrefs)
    with open(output_dir + '/_metadata.json', 'w') as out:
        json.dump(metadata, out)


# def process_sample(sample, sample_output_dir, param_values):
def process_sample(output_dir, param_values):
    ##############################
    result = payload(param_values)
    with open(output_dir + '/payload_result.txt','w') as out:
         out.write(str(result))
    ##############################
    # demonstration output of param_values table
    with open(output_dir + '/appsessionparams.csv','w') as out:
        for key, value in param_values.iteritems():
            if key != 'input.samples':
                out.write('%s,%s\n' % (key ,value))


def process_appsession(param_values):
    project_id = param_values.get('input.project_id')
    samples = param_values.get('input.samples')
    sampleshrefs = [sample['href'] for sample in samples]

    # for sample in samples:
    #     sample_output_dir = '/data/output/appresults/%s/%s' % (project_id, sample['name'])
    #     os.system('mkdir -p "%s"' % sample_output_dir)
    #     process_sample(sample, sample_output_dir, param_values)
    #     write_sample_metadata(sample['name'], 'Sample Description', appsessionhref, sampleshrefs, sample_output_dir)

    output_dir = '/data/output/appresults/%s/sessionsummary' % project_id
    os.system('mkdir -p "%s"' % output_dir)
    process_sample(output_dir, param_values)
    write_metadata('sessionsummary','Session Description', appsessionhref, sampleshrefs, output_dir)


# this file executed as script
##############################

if __name__ == '__main__':
    appsessionhref, appsessionparams = read_appsession(appsession_location)
    param_values = parse_appsessionparams(appsessionparams)
    process_appsession(param_values)
