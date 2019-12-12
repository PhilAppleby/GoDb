db.runCommand( { dropIndexes: "samples", index: ["sample_id_1", "list_posn_1", "assaytype_1_list_posn_1"] })
//
db.samples.ensureIndex({"sample_id": 1})
//
db.samples.ensureIndex({"list_posn": 1})
//
db.samples.ensureIndex({"assaytype":1,"list_posn": 1})

