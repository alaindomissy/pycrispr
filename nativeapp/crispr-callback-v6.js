
function launchSpec(dataProvider)
{
    var project = dataProvider.GetProperty("Input.project-id");
    var projectId = dataProvider.GetProperty("Input.project_id").Id;
    var appResultPathArg = "/data/output/appresults/" + project.Id + "/results";

    var retval = {
        //commandLine: ["/opt/illumina/list-genomes/app.sh", "/genomes", appResultPathArg],
        //commandLine: ["/opt/illumina/list-genomes/app.sh", "/data", appResultPathArg],

        //commandLine: ["echo", "Hello"],
        commandLine: ["/root/bin/appgithub" ],

        //containerImageId: "rwentzel/list-genomes"
        containerImageId: "alaindomissy/crispreating",

        Options: [ "bsfs.enabled=false" ]

    };

    return retval;
}




// original code for list-genomes
/*
function launchSpec(dataProvider)
{
    var project = dataProvider.GetProperty("Input.project-id");
    var appResultPathArg = "/data/output/appresults/" + project.Id + "/results";

    var retval = {
        //commandLine: ["/opt/illumina/list-genomes/app.sh", "/genomes", appResultPathArg],
        //commandLine: ["/opt/illumina/list-genomes/app.sh", "/data", appResultPathArg],

        commandLine: ["echo", "Hello"],

        //containerImageId: "rwentzel/list-genomes"
        containerImageId: "alaindomissy/crispreating"
    };

    return retval;
}

*/





// example multi-node launch spec
/*
function launchSpec(dataProvider)
{
    var ret = {
        nodes: []
    };

    ret.nodes.push({
        appSessionName: "Hello World 1",
        commandLine: [ "cat", "/illumina.txt" ],
        containerImageId: "basespace/demo",
        Options: [ "bsfs.enabled=true" ]
    });

    ret.nodes.push({
        appSessionName: "Hello World 2",
        commandLine: [ "cat", "/illumina.txt" ],
        containerImageId: "basespace/demo",
        Options: [ "bsfs.enabled=true" ]
    });

    return ret;
}
*/

/*
function billingSpec(dataProvider) {
    return [
    {
        "Id" : "insert product ID here",
        "Quantity": 1.0
    }];
}
*/