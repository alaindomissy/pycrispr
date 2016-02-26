function launchSpec(dataProvider)
{
    var project = dataProvider.GetProperty("Input.project-id");
    var projectId = dataProvider.GetProperty("Input.project_id").Id;
    var appResultPathArg = "/data/output/appresults/" + project.Id + "/results";

    var retval = {


        //commandLine: ["/opt/illumina/list-genomes/app.sh", "/genomes", appResultPathArg],
        //commandLine: ["/opt/illumina/list-genomes/app.sh", "/data", appResultPathArg],
        //commandLine: ["/root/bin/listgenomes", "/genomes", appResultPathArg],
        //commandLine: ["/root/bin/listgenomes", "/data", appResultPathArg],



        commandLine: ["/BACKEND/appgithub.sh" ],

        //commandLine: ["python", "/BACKEND/demo.py],

        //commandLine: ["echo", "Hello"],


        containerImageId: "alaindomissy/docker-crispr",

        //containerImageId: "docker.illumina.com/alaindomissy/crispr-push-20160209-hg18-ecoli-phix",
        //containerImageId: "rwentzel/list-genomes",



        Options: [ "bsfs.enabled=false" ]

    };

    return retval;
}
