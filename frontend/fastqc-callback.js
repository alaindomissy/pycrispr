
function launchSpec(dataProvider)
{
    var retval = {
            commandLine: ["python", "home/apps/fastqc/fastqc.py"],
            containerImageId : "mtyagi/fastqc",
            Options: [ "bsfs.enabled=false" ]
            };

    return retval;
}