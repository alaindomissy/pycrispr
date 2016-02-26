
function launchSpec(dataProvider)
{

    var project = dataProvider.GetProperty("Input.project-id");

    var projectId = dataProvider.GetProperty("Input.project_id").Id;

    var appResultPathArg = "/data/output/appresults/" + project.Id + "/results";

    result = {

        commandLine: ["python", "-m", "basespaceapp.app", "/data/", "default_payload"],

        containerImageId: "alaindomissy/basespaceapp",

        Options: [ "bsfs.enabled=false" ]

    };

    return result

}
