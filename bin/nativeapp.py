import os
import json

appsession_location = '/media/ros/crispreating/crispreating/tests/data/input/AppSession.json'

with open(appsession_location) as hdl:
    appsession = json.load(hdl)

parameters = appsession['Properties']['Items']

arguments_with_content = ['Input.binding_interference_spacing',
                          'Input.blast_chunk_size',
                          'Input.blast_max_hsps',
                          'Input.genome-id',
                          'Input.genomic-coord']

param_values = {}
param_values.update(
    {param['Name'].lower(): param['Content']
     for param in parameters
     if param['Name']in arguments_with_content
    }
)
param_values.update(
    {param['Name'].lower(): param['Items']
     for param in parameters
     if param['Name']=='Input.restriction_enzymes'
    }
)
param_values.update(
    {param['Name'].lower(): param['Content']['Id']
     for param in parameters
     if param['Name']=='Input.project_id'
    }
)
param_values.update(
    {'input.projects_ids': project['Id']
     for param in parameters
     if param['Name']=='Input.Projects'
     for project in param['Items']
    }
)
param_values.update(
    {'output.projects_ids': project['Id']
     for param in parameters
     if param['Name']=='Output.Projects'
     for project in param['Items']
    }
)
param_values.update(
    {'input.samples':
        [{'id': sample['Id'], 'href':sample['Href'], 'name': sample['Name']}
         for param in parameters
         if param['Name']=='Input.Samples'
         for sample in param['Items']
        ]
    }
)


metadatatemplate = {"Properties": [{"Items": [], "Type": "sample[]", "Name": "Input.Samples"}],
                    "Name": "",
                    "HrefAppSession": "",
                    "Description": ""}


def json_deepcopy(obj):
    return json.loads(json.dumps(obj))



project_id = param_values['input.project_id']
samples = param_values['input.samples']
hrefs = [sample['href'] for sample in samples]


for sample in samples:

    sample_output_dir = '/data/output/appresults/%s/%s' % (project_id, sample['name'])

    os.system('mkdir -p "%s"' % sample_output_dir)
    #########################
    # RUN CRISPREATING HERE
    #########################
    with open(sample_output_dir + '/parameters.csv','w') as out:
        for key, value in param_values.iteritems():
            if key != 'input.samples':
                out.write('%s,%s\n' % (key ,value))

    metadata = json_deepcopy(metadatatemplate)
    metadata['Name'] = sample['id']
    metadata['Description'] = 'Sample Description'
    metadata['HrefAppSession'] = appsession['Href']
    metadata['Properties'][0]['Items'].extend(hrefs)

    with open(sample_output_dir + '/_metadata.json', 'w') as out:
        json.dump(metadata,out)
