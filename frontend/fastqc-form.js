{
    "$type": "Form",
    "fields": [
        
        {
            "$type": "SectionBreak",
            "title": "Output"
        },
        {
            "$type": "TextBox",
            "size": 300,
            "value": "default value!",
            "label": "App Session Name",
        	"required": true,
        	"requiredMessage": "Please enter name for your app session.",
        	"id": "app-session-name"
        },
        {
            "$type": "ProjectChooser",
            "size": 250,
            "valueType": "Output",
            "allowedPermissions": "owner",
            "label": "Output Project",
            "required": true,
            "requiredMessage": "Please choose a project",
            "id": "project-id"
        },
        {
            "$type": "SectionBreak",
            "title": "Input"
        }, 
        {
            "$type": "SampleChooser",  
            "size": 300, 
            "allowedPermissions": "read", 
            "multiselect": false,
            "label": "Input Sample",
            "required": true,
            "requiredMessage": "Please choose a sample",
            "id": "sample-id"          
        },
        {
            "$type": "SectionBreak",
            "title": "Options"
        },
        {
            "$type": "Select",
            "id": "kmers",
            "label": "Kmer Size",
        	"multiselect": false,
            "helpText" : "Specifies the length of Kmer to look for in the 'Kmer content' module.  Note that large values may make the HTML reports unresponsive.",
        	"choices": [
        		{
        			"value": 2,
        			"text": "2",
        			"selected": false
        		},
                {
            		"value": 3,
        			"text": "3",
        			"selected": false
        		},
                {
            		"value": 4,
        			"text": "4",
        			"selected": false
        		},
                {
            		"value": 5,
        			"text": "5",
        			"selected": true
        		},
                {
            		"value": 6,
        			"text": "6",
        			"selected": false
        		},
                {
            		"value": 7,
        			"text": "7",
        			"selected": false
        		},
                {
            		"value": 8,
        			"text": "8",
        			"selected": false
        		},
                {
            		"value": 9,
        			"text": "9",
        			"selected": false
        		},
        		{
        			"value": 10,
        			"text": "10",
        			"selected": false
        		}
        	]
        },
        
        {
            "$type": "CheckBox",
            "id": "use-contaminants",
        	"label": "Use Contaminant Filter",
            "helpText" : "Screen overrepresented sequences against the default contaminant list",
        	"choices": [
        		{
        			"value": 1,
                    "checked": true
        		}
        	]
        },
                {
            "$type": "SectionBreak",
            "title": "BaseSpace Labs Disclaimer"
        },
        {
            "$type": "CheckBox",
        	"id": "basespace-labs-disclaimer",
            "required": true,
            "requiredMessage": "You must agree to the BaseSpace Labs Disclaimer before using this application.",
        	"label": "BaseSpace Labs",
        	"choices": [
        		{
        			"value": "Accepted",
        			"label": "I acknowledge and agree that (i) this is a BaseSpace Labs App, (ii) I am using it AS-IS without any warranty of any kind, (iii) Illumina has no obligation to provide any technical support for this App, and (iv) Illumina has no liability for my use of this App, including without limitation, any loss of data, incorrect results, or any costs, liabilities, or damages that result from use of this App."
        		}
        	]
        }
        
    ],
    "id": "form-container"
}