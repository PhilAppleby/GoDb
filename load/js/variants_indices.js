// Drop indexes - don't need to do this if the relevant collections are dropped
db.runCommand( { dropIndexes: "variants", index: ["rsid_1_assaytype_1", "chromosome_1_position_1_assaytype_1", "assaytype_1"] })

// (Re)build them
db.variants.ensureIndex({"rsid": 1, "assaytype":1})  // marker id, assaytype
db.variants.ensureIndex({"chromosome": 1, "position":1, "assaytype":1})  // chromosome and position 
db.variants.ensureIndex({"assaytype":1})  // assaytype
