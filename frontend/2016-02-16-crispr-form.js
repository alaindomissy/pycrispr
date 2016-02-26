
{"$type": "Form",
 "fields": [
// GENOME
        {"$type": "SectionBreak","title": "Genome"},
        {"$type": "Select",
         "id": "refgenome_id",/*"label": "Reference Genome",*/
         "choices": [

          {"value":"hg18","text":"Human (UCSC hg18)","selected": "false"},

          {"value":"mm8","text":"Mouse (UCSC mm8)","selected": "false"},

          {"value":"---","text":" ","selected": "false", "disabled":"disabled"},

          {"value":"saccer3","text": "S. Cerevisia (UCSC saccer3)","selected": "false"},

          {"value":"---","text":" ","selected": "false", "disabled":"disabled"},

          {"value":"ecoli","text": "E. Coli K12 strain DH10B (NCBI 2008-03-17)","selected": "true"},
          //k12dh10b"

          {"value":"lambda","text": "Enterobacteriophage Lambda (NCBI 1993-04-28)","selected": "false"},
          {"value":"mycotube","text": "Mycobacterium Tuberculosis H37RV (NCBI 2001-09-07)","selected": "false"},
          {"value":"rhodobacter","text": "Rhodobacter  (NCBI 2005-10-07)","selected": "false"},
          {"value":"staphiaureus","text": "S. Aureus NCTC 8325 (NCBI 2006-02-13)","selected": "false"},
          {"value":"phix","text": "Phi X (Illumina)","selected": "false"},
/*,
          {"value":"hg38","text":"Human (UCSC hg38)","selected": "false"},
          {"value":"hg19","text":"Human (UCSC hg19)","selected": "false"},

          {"value":"mm10","text":"Mouse (UCSC mm10)","selected": "false"},
          {"value":"mm9","text":"Mouse (UCSC mm9)","selected": "false"},

          {"value":"B. Taurus","text": "B. Taurus  (Ensembl UMD3.1)","selected": "false"},
          {"value":"Rat","text": "Rat (UCSC RN5)","selected": "false"},

          {"value":"xl","text": "Xenopus Laevis","selected": "false"},
          {"value":"dm6","text": "Drosophila Melanogaster (UCSC dm6)","selected": "false"},
          {"value":"ce10","text": "Caenorhabditis Elegans(UCSC ce10)","selected": "false"},

          {"value":"tair10","text": "Arabidopsis Thaliana (NCBI tair10)","selected": "false"},

          {"value":"ecoli","text": "E. Coli K12 strain MG1655 (NCBI 2001-10-15 )","selected": "false"},
                   //k12mg1655

*/        ]},
// TARGET
    {"$type": "SectionBreak","title": "Target"},
    {"$type": "RadioButton",
        "id": "use-coordinates","label":"defined by:",
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

          // mm8
          //"value": "chr1:74080002-74084000",

          //"value": "chr3:92280001-92280100",
          //"value": "chr3:92280001-92320000",
          //"value": "chr1:74080002-74120000",



          // mm8 OLD
          //"value": "chr6:47000001-47000200"},
          //"value": "chr6:47599949-47640339"},

          // hg18

          //"value": "chr21:42711601-42715600",
          //"value": "chr21:42711601-4251600"

          //"value": "chr21:42711557-42751600",
          //"value": "chr21:42711558-42742369",
          //"value": "chr21:42462433-42497269",

          // hg18 OLD
          //"value": "chr21:42905101-42965196"}

          //ecoli
          "value":"chr:101001-105000"
          //"value":"chr:101001-121000"
          //"value":"chr:101001-141000"
    },
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

// SAVE AS ANALYSIS / PROJECT
    {"$type":"SectionBreak","title":"Results"},
    {"$type":"TextBox",
        "id":"app-session-name","label":"analysis name",
        "required":true,"requiredMessage":"Please enter name for your app session.",
        "size":400,"minLength":0,"maxLength":150,
        //"value": "CrisprEating Librarycd"},
        //"value": "CrisprEating Librarycd  [LocalDateTime]"
        },
    {"$type": "ProjectChooser",
        "id": "project_id","label": "under project",
        "required":true,"requiredMessage":"Please choose a project",
        "size":400, "valueType":"Output",
        "allowedPermissions": "owner","allowResourceCreation": true,"rules": "is-project-owner"},



    // options
    {"$type": "FieldSet",
    "id": "options","label":"Options",
    "isCollapsible": true,
    "isOpen": true,
    "helpText": "options explanations.",

    "fields": [
        //digest
        {"$type": "CheckBox",
        "id": "restr_enzymes","label": "Restriction Enzymes",
        "choices": [{"value": "BfaI","label": "BfaI","checked": true},
        		    {"value": "ScrFI","label": "ScrFI","checked": true},
                    {"value": "HpaII","label": "HpaII","checked": true}],
        "helpText":"Restriction enzymes, also known as restriction endonucleases,
        are enzymes that cut a DNA molecule at a particular place.
        They are essential tools for recombinant DNA technology.
        The enzyme 'scans' a DNA molecule,
        looking for a particular sequence, usually of four to six nucleotides."
        },
        //blast
        {"$type": "Numeric",
            "id": "blast_chunk_size","label": "Blast Chunk Size (10 to 100)",
            "size": 50, "required":false,
            "min": 1,"max": 100,"value":10,
            "numericType": "Integer",
            "rules": "is-it-a-pi"
        },
        {"$type": "Numeric",
            "id": "blast_max_hsps","label": "Blast Max HSPs (10 to 1000)",
            "size": 50,"required": false,
        	"min": 1,"max":100,"value":10,
            "numericType":"Integer",
            "rules": "is-it-a-pi"
        },
        //score
        {"$type": "Numeric",
            "id": "interference_gap","label": "Binding Interference Spacing (10 to 40)",
            "size": 50,"required":false,
            "min": 10,"max": 40,"value": 20,
            "numericType": "Integer",
             "rules": "is-it-a-pi"
        },
        //cluster

        //prime
        {"$type": "RadioButton",
        "id": "primer_screening","label":"primer screening method:",
        "value":"dumb",
        "choices": [{"value":"dumb","label":"Dumb"},
                    {"value":"epcr","label":"In Silico PCR"},
                   ]
        },

    ]},
    //logs
    {"$type": "FieldSet",
    "id": "logging-options","label":"Logging",
    "isCollapsible": true,
    "isOpen": true,
    "helpText": "",
    "fields": [
        {
            "$type": "CheckBox",
        	"id": "checkbox-logging",
        	"label": "Logging options",
        	"choices": [
                		{"value":"logdigeston","label":"digest","checked":true},
                		{"value":"logblaston","label":"blast","checked":true},
                        {"value":"logscoreon","label":"score","checked":true},
                		{"value":"logclusterton","label":"cluster","checked":true},
                        {"value":"logprimeton","label":"prime","checked":true},
	                   ]
        }
    ]
    },
    //caches
    {"$type": "FieldSet",
    "id": "reusing-options","label":"Reusing",
    "isCollapsible": true,
    "isOpen": true,
    "helpText": "",
    "fields": [
        {
            "$type": "CheckBox",
            "id": "checkbox-reusing",
        	"label": "Reusing options",
        	"choices": [
                		{"value":"logdigeston","label":"digest","checked":false},
                		{"value":"logblaston","label":"blast","checked":false},
                        {"value":"logscoreon","label":"score","checked":false},
                		{"value":"logclusterton","label":"cluster","checked":false},
                        {"value":"logprimeton","label":"prime","checked":false},
	                   ]
        }
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
