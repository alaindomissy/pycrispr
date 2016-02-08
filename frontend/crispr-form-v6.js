


{"$type": "Form",
 "fields": [
//////////
        {"$type": "SectionBreak","title": "Genome"},
        {"$type": "Select",
         "id": "genome_id",/*"label": "Reference Genome",*/
         "choices": [
          {"value":"hg38","text":"Human (UCSC hg38)","selected": "false"},
          {"value":"hg19","text":"Human (UCSC hg19)","selected": "false"},
          {"value":"hg18","text":"Human (UCSC hg18)","selected": "false"},
          {"value":"mm10","text":"Mouse (UCSC mm10)","selected": "false"},
          {"value":"mm9","text":"Mouse (UCSC mm9)","selected": "false"},
          {"value":"mm8","text":"Mouse (UCSC mm8)","selected": "true"},

          {"value":"phix","text": "Phi X (Illumina)","selected": "false"}

/*,       {"value":"xl","text": "Xenopus Laevis","selected": "false"},
          {"value":"dm6","text": "Drosophila Melanogaster (UCSC dm6)","selected": "false"},
          {"value":"ce10","text": "Caenorhabditis Elegans(UCSC ce10)","selected": "false"},
          {"value":"tair10","text": "Arabidopsis Thaliana (NCBI tair10)","selected": "false"},
          {"value":"yeast","text": "S. Cerevisia (UCSC sacCer2)","selected": "false"},
          {"value":"ecoli","text": "E. Coli DH10B (NCBI 2008-03-17)","selected": "false"},
          {"value":"ecoli2","text": "E. Coli MG1655 (NCBI 2001-10-15 )","selected": "false"},

          {"value":"B. Taurus","text": "B. Taurus  (Ensembl UMD3.1)","selected": "false"},
          {"value":"S. Aureus NCTC 8325","text": "S. Aureus NCTC 8325 (NCBI 2006-02-13)","selected": "false"},
          {"value":"Rat","text": "Rat (UCSC RN5)","selected": "false"},
          {"value":"Rhodobacter","text": "Rhodobacter  (NCBI 2005-10-07)","selected": "false"},

*/        ]},
//////////
    {"$type": "SectionBreak","title": "Target"},
    {"$type": "RadioButton",
        "id": "use-coordinates","label":"Target specification by",
        "value":1,
        "choices": [{"value":1,"label":"Coordinates"},
            {"value":4,"label":"File"},
            {"value":2,"label":"Sample"},
            {"value":3,"label":"App Result"}],
        "togglers": [{"toggleOnValue":1,"toggleFields":"genomic_coord"},
            {"toggleOnValue":4,"toggleFields":"file_id"},
            {"toggleOnValue":2,"toggleFields":"sample_id"},
            {"toggleOnValue":3,"toggleFields":"appresult_id"}]
    },


    {"$type": "TextBox",
          "id": "genomic_coord",/*"label":"Genomic Coordinates",*/
          "required":true,"requiredMessage": "Please enter genomic goordinates for your target.",
          "size": 400, "minLength": 0, "maxLength": 150,
          "value": "chr6:47000001-47000200"},
          //"value": "chr21:42905101-42965196"},
          //"value": "chr6:47599949-47640339"},


    {"$type":"FileChooser",
        "id":"file_id",/*"label":"File",*/
        "required":false,/*"requiredMessage": "Please choose a file",*/
	    "size":400,"valueType":"Input","multiselect":false,
	    "extensionFilters": ".fasta,.fastq,.fasta.gz,.fastq.gz",},
    {"$type":"SampleChooser",
        "id":"sample_id",/*"label":"Alternatively, select an existing Sample",*/
        "required":false,"size":400,"valueType": "Input","multiselect":true,
/*        "displayFields":[
            {"$type": "CheckBox","id": "is-sample-awesome","label": "Awesome?",
                "choices":[{"label":"Y","value":1},{"label":"Y+","value":2}]},
            {"$type": "RadioButton","id": "is-sample-wicked","label": "Wicked?",
                "choices":[{"label":"Y","value":1},{"label":"N","value":0},{"label":"M","value":2}]},
            {"$type": "Select","id": "is-sample-extraordinary","label": "Extraordinary?",
            "multiselect": "true","size": 3,
            "choices":[{"text":"Oh, Yes!","value":1},{"text":"Oh, No!","value":0},{"text":"Close Contender","value" 0.5}]}],
*/      "allowedPermissions":"read","rules":"sample-reader"},
    {"$type": "AppResultChooser",
        "id": "appresult_id",/*"label":"Alternatively, select an existing app Result",*/
        "required": false,"size":400,"valueType":"Input" ,"multiselect":"true" },

//////////
    {"$type":"SectionBreak","title":"Save results"},
    {"$type":"TextBox",
        "id":"app-session-name","label":"as analysis name",
        "required":true,"requiredMessage":"Please enter name for your app session.",
        "size":400,"minLength":0,"maxLength":150,
        "value": "CrisprEating Library [LocalDateTime]"},
    {"$type": "ProjectChooser",
        "id": "project_id","label": "under project",
        "required":true,"requiredMessage":"Please choose a project",
        "size":400, "valueType":"Output",
        "allowedPermissions": "owner","allowResourceCreation": true,"rules": "is-project-owner"},
/*  {"$type": "RadioButton",
        "id": "genome_id","label": "Reference Genome",
        "value": "hg38",
        "choices": [{"value": "hg38","label": "hg38"},
        	         {"value": "mm8","label": "mm8"},
                     {"value": "xl","label": "xl"}]},
    {"$type": "RadioButton",
        "id": "fileformat","label": "File Format",
        "value": "fastq.gz",
        "choices": [{"value": "fastq.gz","label": "fastq.gz"},
        	{"value": "fasta.gz","label": "fasta.gz"},
            {"value": "fastq","label": "fastq"},
        	{"value": "fasta","label": "fasta"}]},
*/
//////////
    {"$type": "SectionBreak","title": "Advanced"},
    {"$type": "FieldSet",
    "id": "digest-options","label":"Digest options",
    "isCollapsible": true,
    "isOpen": false,
    "helpText": "Digest options explanations.",
    "fields": [
        {"$type": "CheckBox",
        "id": "restriction_enzymes","label": "Restriction Enzymes",
        "choices": [{"value": "BfaI","label": "BfaI","checked": true},
        		    {"value": "ScrFI","label": "ScrFI","checked": true},
                    {"value": "HpaII","label": "HpaII","checked": true}],
        "helpText":"Restriction enzymes, also known as restriction endonucleases, are enzymes that cut a DNA molecule at a particular place. They are essential tools for recombinant DNA technology. The enzyme 'scans' a DNA molecule, looking for a particular sequence, usually of four to six nucleotides."},
    ]},
    {"$type": "FieldSet",
    "id": "blast-options","label":"Blast options",
    "isCollapsible": true,
    "isOpen": false,
    "helpText": "Blast options explanations.",
    "fields": [
        {"$type": "Numeric",
            "id": "blast_chunk_size","label": "Blast Chunk Size (10 to 100)",
            "size": 50, "required":false,
            "min": 10,"max": 1000,"value": 25,
            "numericType": "Integer",
            "rules": "is-it-a-pi"},
        {"$type": "Numeric",
            "id": "blast_max_hsps","label": "Blast Max HSPs (10 to 1000)",
            "size": 50,"required": false,
        	"min": 10,"max":1000,"value":10,
            "numericType":"Integer",
            "rules": "is-it-a-pi"}
    ]
    },
    {"$type": "FieldSet",
    "id": "score-options","label":"Score options",
    "isCollapsible": true,
    "isOpen": false,
    "helpText": "Score options explanations.",
    "fields": [
        {"$type": "Numeric",
            "id": "binding_interference_spacing","label": "Binding Interference Spacing (20 to 40)",
            "size": 50,"required":false,
            "min": 20,"max": 40,"value": 25,
            "numericType": "Integer",
             "rules": "is-it-a-pi"}
    ]
    },
    {"$type": "FieldSet",
    "id": "cluster-options","label":"Cluster options",
    "isCollapsible": true,
    "isOpen": false,
    "helpText": "Cluster options explanations.",
    "fields": [
    ]
    },
    {"$type": "FieldSet",
    "id": "primers-options","label":"Primers options",
    "isCollapsible": true,
    "isOpen": false,
    "helpText": "Primers options explanations.",
    "fields": [
    ]
    },
    {"$type": "FieldSet",
    "id": "logging-options","label":"Logging options",
    "isCollapsible": true,
    "isOpen": false,
    "helpText": "Logging options explanations.",
    "fields": [
    ]
    },
//////////
/*  {"$type": "SectionBreak","title": "Child Check"},
    {"$type": "CheckBox",
	    "id": "disable-child-check","label": "Disable Child Check",
	    "choices": [{"value":1,"label":"disable"}]},
    {"$type": "TextBox","id": "child-check","label": "Child Check",
        "required": false, "requiredMessage": "This is not reaaly required.",
	    "size": 400, "value": "it's a girl",
	    "helpText": "this help text will display next to the control in a pop over"}
*/
],
"rulesets":[
    {"$type": "PermissionValidationRule",
         "id": "sample-reader", "message": "You do not have read access to the selected sample",
         "permissions": "Read", "severity": "Error"},
    {"$type": "PermissionValidationRule",
         "id": "is-project-owner","message": "You aren't the owner of the selected project.",
         "permissions": "Own", "severity": "Error"}
    ]
}






// original form for list-genomes
/*
{
    "$type": "Form",
    "fields": [
        {
            "$type": "TextBox",
            "size": 400,
            "minLength": 0,
            "maxLength": 150,
            "value": "Example [LocalDateTime]",
            "label": "Analysis Name",
            "required": true,
            "requiredMessage": "Please enter name for your app session.",
            "id": "app-session-name"
        },
//
        {
            "$type": "SampleChooser",
            "size": 300,
            "valueType": "Input",
            "allowedPermissions": "read",
            "label": "Sample",
            "required": false,
            "id": "sample-id",
            "rules": "sample-reader"
        },
//
        {
            "$type": "ProjectChooser",
            "size": 300,
            "valueType": "Output",
            "allowedPermissions": "owner",
            "label": "Save Results To",
            "required": true,
            "requiredMessage": "Please choose a project",
            "id": "project-id",
            "allowResourceCreation": true,
            "rules": "is-project-owner"
        },
        {
            "$type": "SectionBreak"
        }
    ],
    "rulesets":[
//
        {
            "$type": "PermissionValidationRule",
            "permissions": "Read",
            "severity": "Error",
            "message": "You do not have read access to the selected sample",
            "id": "sample-reader"
        },
//
        {
            "$type": "PermissionValidationRule",
            "permissions": "Own",
            "severity": "Error",
            "message": "You aren't the owner of the selected project.",
            "id": "is-project-owner"
        }
    ]
}
*/