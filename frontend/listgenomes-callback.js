


function launchSpec(dataProvider)
{
    var project = dataProvider.GetProperty("Input.project-id");
    var appResultPathArg = "/data/output/appresults/" + project.Id + "/results";

    var retval = {
        //commandLine: ["/opt/illumina/list-genomes/app.sh", "/genomes", appResultPathArg],
        commandLine: ["/opt/illumina/list-genomes/app.sh", "/data", appResultPathArg],
        containerImageId: "rwentzel/list-genomes"
    };

    return retval;
}
