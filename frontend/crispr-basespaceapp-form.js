{   "$type": "Form",
    "fields": [
        {
                "$type": "TextBox", "id": "textbox-1", "label": "Textbox Text",
                "size": 100,
            	"value": "replace me",
            	"required": true,
            	"requiredMessage": "Please enter textbox text.",
            	"helpText": "this help text will display next to the control in a pop over"
        },
        {
                "$type": "Numeric", "id": "numeric-1",     "label": "Located in how many circles",
            	"size": 50,
            	"required": false,
            	"min": -6.28,
            	"max": 6.28,
            	"value": 2,
            	"numericType": "FloatingPoint",
            	"rules": "is-it-a-pi"
        },
        {
                "$type": "RadioButton",
                "id": "radio-1",
            	"label": "Stranded?",
            	"value": 1,
            	"choices": [
            		{ "value": 0,  "label": "No"    },
            		{ "value": 1,  "label": "Yes"   }
            	]
        },
        {
            "$type": "CheckBox", "id": "checkbox-1", "label": "checkbox control",

        	"choices": [
        		{ "value": 10, "label": "ten"	                  },
        		{ "value": 50, "label": "fifty", "checked": true  }
        	]
        },


        {
            "$type": "SectionBreak"
        },


        {   "$type": "TextBox", "label": "Analysis Name",
            "size": 400,
            "minLength": 0,
            "maxLength": 150,
            "value": "Example [LocalDateTime]",
            "required": true,
            "requiredMessage": "Please enter name for your app session.",
            "id": "app-session-name"
        },

        {   "$type": "SampleChooser", "label": "Sample", "id": "sample-id",

            "size": 300,
            "valueType": "Input",
            "allowedPermissions": "read",
            "required": false,
            "rules": "sample-reader"
        },

        {   "$type": "ProjectChooser", "id": "project-id", "label": "Save Results To",

            "size": 300,
            "valueType": "Output",
            "allowedPermissions": "owner",

            "required": true,
            "requiredMessage": "Please choose a project",
            "allowResourceCreation": true,
            "rules": "is-project-owner"
        },

        {
            "$type": "SectionBreak"
        }
    ],




    "rulesets":[
        {
            "$type": "PermissionValidationRule",
            "permissions": "Read",
            "severity": "Error",
            "message": "You do not have read access to the selected sample",
            "id": "sample-reader"
        },
        {
            "$type": "PermissionValidationRule",
            "permissions": "Own",
            "severity": "Error",
            "message": "You aren't the owner of the selected project.",
            "id": "is-project-owner"
        }
    ]
}
