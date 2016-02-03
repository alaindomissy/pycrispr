########################################################################################################################
#
# APP
#
# API functions:
#
# intended to be run as a scipt (if __name__ == '__main__': code)
#
########################################################################################################################

from __future__ import print_function
import os
import json
from datetime import datetime

from pybasespace.payload import payload
from crispr.config import APPSESS


def read_appsession(appsession_jsonfilename):
    with open(appsession_jsonfilename) as hdl:
        appsession = json.load(hdl)
    appsessionparams = appsession['Properties']['Items']
    appsessionhref = appsession['Href']
    return appsessionhref, appsessionparams


def parse_appsessionparams(appsessionparams):
    arguments_with_content = ['input.interference_spacing',
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
         if param.get('Name').lower() == 'input.restriction_enzymes'
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
         if param.get('Name') == 'Input.Projects'
         for project in param['Items']
        }
    )
    param_values.update(
        {'output.projects_ids': project.get('Id')
         for param in appsessionparams
         if param.get('Name') == 'Output.Projects'
         for project in param.get('Items')
        }
    )
    param_values.update(
        {'input.samples':
            [{'id': sample['Id'], 'href':sample['Href'], 'name': sample['Name']}
             for param in appsessionparams
             if param.get('Name') == 'Input.Samples'
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
        json.dump(metadata, out, indent=4, sort_keys=True)
    # TODO this is not pretty ptinting the object, fix
    # with open(output_dir + '/metadata.json', 'w') as out:
    #     json.dump(metadata, out, indent=4, sort_keys=True)


# def process_sample(sample, sample_output_dir, param_values):
def process_sample(output_dir, param_values):

    print("start process sample time : ",  datetime.now().isoformat('_'))

    # the application payload - CAREFUL this also creates the output_dir
    ####################################################################################################################
    result = payload(param_values, output_dir)
    ####################################################################################################################

    with open(output_dir + '/payload_result.txt','w') as out:
          out.write(str(result))
    # TODO get some results more interesting to save
    # print('\npayload_result.txt printed to: %s' % output_dir)


    # output of param_values table
    with open(output_dir + '/appsessionparams.csv','w') as out:
        for key, value in param_values.iteritems():
            if key != 'input.samples':
                out.write('%s\t\t%s\n' % (key,value))
    print('\nappsessionparams.csv printed to %s' % output_dir)


def process_appsession(param_values):
    project_id = param_values.get('input.project_id')
    samples = param_values.get('input.samples')
    sampleshrefs = [sample['href'] for sample in samples]

    # for samoutput_dirple in samples:
    #     sample_output_dir = '/data/output/appresults/%s/%s' % (project_id, sample['name'])
    #     os.system('mkdir -p "%s"' % sample_output_dir)
    #     process_sample(sample, sample_output_dir, param_values)
    #     write_sample_metadata(sample['name'], 'Sample Description', appsessionhref, sampleshrefs, sample_output_dir)

    output_dir = '/data/output/appresults/' + project_id + '/sessionsummary_' + datetime.now().isoformat('_') + '/'

    # see above output_dir get created by call to payload
    # os.system('mkdir -p "%s"' % output_dir)

    process_sample(output_dir, param_values)
    write_metadata('\nsessionsummary','Session Description', appsessionhref, sampleshrefs, output_dir)

    print("end process sample time : ",  datetime.now().isoformat('_'))


# this file executed as script
##############################

if __name__ == '__main__':
    appsessionhref, appsessionparams = read_appsession(APPSESS)
    param_values = parse_appsessionparams(appsessionparams)
    process_appsession(param_values)
