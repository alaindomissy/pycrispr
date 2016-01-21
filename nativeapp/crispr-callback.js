function launchSpec(dataProvider){
    var projectId = dataProvider.GetProperty("Input.project_id").Id;
    var ret = {
           //commandLine: ["/root/anaconda3/envs/buffet/bin/python","-m","pybasespace.app" ],
           commandLine: ["/root/app"],
           containerImageId: "alaindomissy/crispreating",
           Options: [ "bsfs.enabled=false" ]
          };
    return ret;
}



/*
function formUpdates(dataProvider)
{
    var checked = ('' + dataProvider.GetProperty("input.disable-child-check")) === "1";
    if (checked) {
        dataProvider.AttributeUpdates.Add({ ElementId: "child-check", AttributeName: "disabled", AttributeValue: "disabled" });
        dataProvider.AttributeUpdates.Add({ ElementId: "child-check", AttributeName: "title", AttributeValue: 'Uncheck "Disable Child" to enable' });
    } else {
        dataProvider.AttributeUpdates.Remove({ ElementId: "child-check", AttributeName: "disabled" });
        dataProvider.AttributeUpdates.Remove({ ElementId: "child-check", AttributeName: "title" });
    }
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