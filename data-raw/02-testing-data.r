

# DRB1 --------------------------------------------------------------------


## DNA version ID 52  ##
HLA_DRB1_TEST_52 <- get_testdata(gene = "DRB1", nextype_basis_id = dna2basis("52"))
tdnew <- HLA_DRB1_TEST_52
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DRB1_TEST_52 <- tdnew2
devtools::use_data(HLA_DRB1_TEST_52, pkg = "~/Devel/HLAdata", overwrite = TRUE)

## DNA version ID 53  ##
HLA_DRB1_TEST_53 <- get_testdata(gene = "DRB1", nextype_basis_id = dna2basis("53"))
tdnew <- HLA_DRB1_TEST_53
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DRB1_TEST_53 <- tdnew2
devtools::use_data(HLA_DRB1_TEST_53, pkg = "~/Devel/HLAdata", overwrite = TRUE)

## DNA version ID 56  ##
HLA_DRB1_TEST_56 <- get_testdata(gene = "DRB1", nextype_basis_id = dna2basis("56"))
tdnew <- HLA_DRB1_TEST_56
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DRB1_TEST_56 <- tdnew2
devtools::use_data(HLA_DRB1_TEST_56, pkg = "~/Devel/HLAdata", overwrite = TRUE)

## DNA version ID 58  ##
HLA_DRB1_TEST_58 <- get_testdata(gene = "DRB1", nextype_basis_id = dna2basis("58"))
tdnew <- HLA_DRB1_TEST_58
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DRB1_TEST_58 <- tdnew2
devtools::use_data(HLA_DRB1_TEST_58, pkg = "~/Devel/HLAdata", overwrite = TRUE)


# DQB1 --------------------------------------------------------------------


## testdata DQB1; basis ID 1267  ##############################################
HLA_DQB1_TEST_1267 <- get_testdata(gene = "DQB1", nextype_basis_id = "1267")
tdnew <- HLA_DQB1_TEST_1267
tdnew <- update_testdata("ID11654815", "6432", 1:2, td = tdnew)
tdnew <- update_testdata("ID11824749", "6216", 2, td = tdnew)
tdnew <- update_testdata("ID11831442", "6366", 2, td = tdnew)
tdnew <- update_testdata("ID11834447", "6221", 2, td = tdnew)
tdnew <- update_testdata("ID11834500", "6226", 2, td = tdnew)
tdnew <- update_testdata("ID11837976", "6370", 1:3, td = tdnew)
tdnew <- update_testdata("ID11841576", "6240", 2, td = tdnew)
tdnew <- update_testdata("ID11842344", "6243", 2, td = tdnew)
tdnew <- update_testdata("ID11843480", "6264", 2, td = tdnew)
tdnew <- update_testdata("ID11845054", "6333", 2, td = tdnew)
tdnew <- update_testdata("ID11846181", "6267", 2, td = tdnew)
tdnew <- update_testdata("ID11851403", "6297", 2, td = tdnew)
tdnew <- update_testdata("ID11853936", "6311", 2, td = tdnew)
tdnew <- update_testdata("ID11855888", "6294", 2, td = tdnew)
tdnew <- update_testdata("ID11857038", "6312", 2, td = tdnew)
tdnew <- update_testdata("ID11858076", "6465", 1, td = tdnew)
tdnew <- update_testdata("ID11858126", "6465", 1, td = tdnew)
tdnew <- update_testdata("ID11859960", "6304", 2, td = tdnew)
tdnew <- update_testdata("ID11860544", "6417", 3, td = tdnew)
tdnew <- update_testdata("ID11861330", "6307", 2, td = tdnew)
tdnew <- update_testdata("ID11864273", "6424", 3, td = tdnew)
tdnew <- update_testdata("ID11864471", "6298", 3, td = tdnew)
tdnew <- update_testdata("ID11866234", "6328", 2, td = tdnew)
tdnew <- update_testdata("ID11870745", "6455", 1, td = tdnew)
tdnew <- update_testdata("ID11874953", "6455", 1, td = tdnew)
tdnew <- update_testdata("ID11875195", "6369", 3, td = tdnew)
tdnew <- update_testdata("ID11875467", "6369", 2, td = tdnew)
tdnew <- update_testdata("ID11875710", "6346", 2, td = tdnew)
tdnew <- update_testdata("ID11879254", "6438", 2, td = tdnew)
tdnew <- update_testdata("ID11880244", "6356", 3, td = tdnew)
tdnew <- update_testdata("ID11882611", "6358", 3, td = tdnew)
tdnew <- update_testdata("ID11883188", "6381", 3, td = tdnew)
tdnew <- update_testdata("ID11884730", "6395", 2, td = tdnew)
tdnew <- update_testdata("ID11886568", "6395", 2, td = tdnew)
tdnew <- update_testdata("ID11887582", "6439", 2, td = tdnew)
tdnew <- update_testdata("ID11888906", "6385", 2, td = tdnew)
tdnew <- update_testdata("ID11890489", "6384", 2, td = tdnew)
tdnew <- update_testdata("ID11897592", "6412", 2, td = tdnew)
tdnew <- update_testdata("ID11907450", "6440", 2, td = tdnew)
tdnew <- update_testdata("ID11907548", "6440", 2, td = tdnew)
tdnew2 <- clean_testing_data(td = tdnew)
testdata_dqb1_1267 <- tdnew2
devtools::use_data(testdata_dqb1_1267, overwrite = TRUE)

## testdata DQB1; basis ID 1412  ##############################################
HLA_DQB1_TEST_1412 <- get_testdata(gene = "DQB1", nextype_basis_id = "1412")
tdnew <- HLA_DQB1_TEST_1412
tdnew <- update_testdata("ID11639946", "7306", 3, td = tdnew)
tdnew <- update_testdata("ID12036005", "7046", 2, td = tdnew)
tdnew <- update_testdata("ID12036580", "7046", 1:4, td = tdnew)
tdnew <- update_testdata("ID12037063", "7047", 1:4, td = tdnew)
tdnew <- update_testdata("ID12037436", "7045", 2, td = tdnew)
tdnew <- update_testdata("ID12038406", "7048", 3, td = tdnew)
tdnew <- update_testdata("ID12043683", "7049", 3, td = tdnew)
tdnew <- update_testdata("ID12043854", "7050", 3, td = tdnew)
tdnew <- update_testdata("ID12050561", "7093", 2, td = tdnew)
tdnew <- update_testdata("ID12051117", "7094", 3, td = tdnew)
tdnew <- update_testdata("ID12051571", "7093", 2, td = tdnew)
tdnew <- update_testdata("ID12052307", "7052", 2, td = tdnew)
tdnew <- update_testdata("ID12052636", "7094", 3, td = tdnew)
tdnew <- update_testdata("ID12052643", "7094", 4, td = tdnew)
tdnew <- update_testdata("ID12052709", "7094", 2, td = tdnew)
tdnew <- update_testdata("ID12052721", "7094", 4, td = tdnew)
tdnew <- update_testdata("ID12053370", "7095", 2, td = tdnew)
tdnew <- update_testdata("ID12053519", "7095", 3, td = tdnew)
tdnew <- update_testdata("ID12053713", "7095", 2, td = tdnew)
tdnew <- update_testdata("ID12053749", "7095", 3, td = tdnew)
tdnew <- update_testdata("ID12054076", "7127", 2, td = tdnew)
tdnew <- update_testdata("ID12054794", "7052", 3, td = tdnew)
tdnew <- update_testdata("ID12070597", "7098", 3, td = tdnew)
tdnew <- update_testdata("ID12072709", "7099", 4, td = tdnew)
tdnew <- update_testdata("ID12075277", "7100", 2, td = tdnew)
tdnew <- update_testdata("ID12081562", "7129", 3, td = tdnew)
tdnew <- update_testdata("ID12083403", "7135", 2, td = tdnew)
tdnew <- update_testdata("ID12091015", "7138", 2, td = tdnew)
tdnew <- update_testdata("ID12094458", "7195", 3, td = tdnew)
tdnew <- update_testdata("ID12095966", "7139", 3, td = tdnew)
tdnew <- update_testdata("ID12097185", "7139", 4, td = tdnew)
tdnew <- update_testdata("ID12097712", "7160", 2, td = tdnew)
tdnew <- update_testdata("ID12100761", "7205", 2, td = tdnew)
tdnew <- update_testdata("ID12103452", "7197", 1:4, td = tdnew)
tdnew <- update_testdata("ID12105538", "7213", 4, td = tdnew)
tdnew <- update_testdata("ID12105994", "7261", 4, td = tdnew)
tdnew <- update_testdata("ID12106885", "7164", 4, td = tdnew)
tdnew <- update_testdata("ID12107223", "7202", 2, td = tdnew)
tdnew <- update_testdata("ID12123744", "7211", 1:4, td = tdnew)
tdnew <- update_testdata("ID12126165", "7211", 3, td = tdnew)
tdnew <- update_testdata("ID12127096", "7466", 2, td = tdnew)
tdnew <- update_testdata("ID12127591", "7200", 2, td = tdnew)
tdnew <- update_testdata("ID12129818", "7218", 3, td = tdnew)
tdnew <- update_testdata("ID12132085", "7270", 2, td = tdnew)
tdnew <- update_testdata("ID12132788", "7278", 2, td = tdnew)
tdnew <- update_testdata("ID12134142", "7219", 3, td = tdnew)
tdnew <- update_testdata("ID12134661", "7434", 2, td = tdnew)
tdnew <- update_testdata("ID12134774", "7220", 3, td = tdnew)
tdnew <- update_testdata("ID12135791", "7222", 2, td = tdnew)
tdnew <- update_testdata("ID12033536", "7044", 4, td = tdnew)
tdnew <- update_testdata("ID12095525", "7203", c(2,4), td = tdnew)
tdnew <- update_testdata("ID12140193", "7310", 1:4, td = tdnew)
tdnew <- update_testdata("ID12156880", "7341", 3, td = tdnew)
tdnew <- update_testdata("ID12162202", "7377", 3, td = tdnew)
tdnew <- update_testdata("ID12162437", "7579", 2, td = tdnew)
tdnew <- update_testdata("ID12163628", "7354", 3, td = tdnew)
tdnew <- update_testdata("ID12166975", "7380", 3, td = tdnew)
tdnew <- update_testdata("ID12167429", "7422", 2, td = tdnew)
tdnew <- update_testdata("ID12172300", "7423", 2, td = tdnew)
tdnew <- update_testdata("ID12179403", "7420", 1, td = tdnew)
tdnew <- update_testdata("ID12181348", "7411", 2, td = tdnew)
tdnew <- update_testdata("ID12183031", "7510", 2, td = tdnew)
tdnew <- update_testdata("ID12184440", "7415", 1:4, td = tdnew)
tdnew <- update_testdata("ID12186499", "7491", 4, td = tdnew)
tdnew <- update_testdata("ID12186603", "7452", 2, td = tdnew)
tdnew <- update_testdata("ID12186980", "7453", 3, td = tdnew)
tdnew <- update_testdata("ID12189756", "7509", 4, td = tdnew)
tdnew <- update_testdata("ID12190924", "7509", 2, td = tdnew)
tdnew <- update_testdata("ID12191441", "7492", 2, td = tdnew)
tdnew <- update_testdata("ID12192708", "7493", 2, td = tdnew)
tdnew <- update_testdata("ID12194767", "7500", 4, td = tdnew)
tdnew <- update_testdata("ID12195186", "7501", 3, td = tdnew)
tdnew <- update_testdata("ID12198838", "7503", 2, td = tdnew)
tdnew <- update_testdata("ID12207434", "7582", 2, td = tdnew)
tdnew <- update_testdata("ID12216504", "7595", 3, td = tdnew)
tdnew <- update_testdata("ID12216656", "7593", 2, td = tdnew)
tdnew <- update_testdata("ID12220103", "7599", 3, td = tdnew)
tdnew <- update_testdata("ID12220362", "7600", 2, td = tdnew)
tdnew <- update_testdata("ID12220794", "7600", 4, td = tdnew)
tdnew <- update_testdata("ID12223847", "7675", 2, td = tdnew)
tdnew <- update_testdata("ID12223879", "7675", 3, td = tdnew)
tdnew <- update_testdata("ID12229862", "7699", 3, td = tdnew)
tdnew <- update_testdata("ID12234406", "7687", 3, td = tdnew)
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DQB1_TEST_1412 <- tdnew2
devtools::use_data(HLA_DQB1_TEST_1412, pkg = "~/Devel/HLAdata", overwrite = TRUE)

## testdata DPB1; basis ID 1267  #############################################
HLA_DPB1_TEST_1267 <- get_testdata(gene = "DPB1", nextype_basis_id = "1267")
tdnew <- HLA_DPB1_TEST_1267
tdnew <- update_testdata("ID11836832", "6256", 2, td = tdnew)
tdnew <- update_testdata("ID11838031", "6235", 2, td = tdnew)
tdnew <- update_testdata("ID11840146", "6237", 2, td = tdnew)
tdnew <- update_testdata("ID11840166", "6237", 2, td = tdnew)
tdnew <- update_testdata("ID11840754", "6241", 2, td = tdnew)
tdnew <- update_testdata("ID11851483", "6297", 2, td = tdnew)
tdnew <- update_testdata("ID11851712", "6297", 2, td = tdnew)
tdnew <- update_testdata("ID11852451", "6297", 2, td = tdnew)
tdnew <- update_testdata("ID11852585", "6293", 2, td = tdnew)
tdnew <- update_testdata("ID11858326", "6295", 2, td = tdnew)
tdnew <- update_testdata("ID11859378", "6456", 3, td = tdnew)
tdnew <- update_testdata("ID11860387", "6299", 2, td = tdnew)
tdnew <- update_testdata("ID11861527", "6307", 2, td = tdnew)
tdnew <- update_testdata("ID11862782", "6377", 2, td = tdnew)
tdnew <- update_testdata("ID11863588", "6308", 2, td = tdnew)
tdnew <- update_testdata("ID11864270", "6424", 1:2, td = tdnew)
tdnew <- update_testdata("ID11867862", "6327", 2, td = tdnew)
tdnew <- update_testdata("ID11870649", "6326", 2, td = tdnew)
tdnew <- update_testdata("ID11870866", "6349", 2, td = tdnew)
tdnew <- update_testdata("ID11875145", "6369", c(2, 5:8), td = tdnew)
tdnew <- update_testdata("ID11876457", "6386", 2, td = tdnew)
tdnew <- update_testdata("ID11878388", "6392", 4, td = tdnew)
tdnew <- update_testdata("ID11879643", "6392", 3, td = tdnew)
tdnew <- update_testdata("ID11882665", "6358", 3, td = tdnew)
tdnew <- update_testdata("ID11882939", "6381", 2, td = tdnew)
tdnew <- update_testdata("ID11886416", "6344", 2, td = tdnew)
tdnew <- update_testdata("ID11886497", "6342", 2, td = tdnew)
tdnew <- update_testdata("ID11889373", "6367", 4, td = tdnew)
tdnew <- update_testdata("ID11889591", "6367", c(2,4), td = tdnew)
tdnew <- update_testdata("ID11889722", "6367", 4, td = tdnew)
tdnew <- update_testdata("ID11897300", "6447", 2, td = tdnew)
tdnew <- update_testdata("ID11899837", "6431", 2, td = tdnew)
tdnew <- update_testdata("ID11905280", "6452", 2, td = tdnew)
tdnew <- update_testdata("ID11840398", "6237", 2, td = tdnew)
tdnew <- update_testdata("ID11657200", "6330", 2, td = tdnew)
tdnew <- update_testdata("ID11812909", "6396", 2, td = tdnew)
tdnew <- update_testdata("ID11823252", "6228", 2, td = tdnew)
tdnew <- update_testdata("ID11845892", "6267", 2, td = tdnew)
tdnew <- update_testdata("ID11864818", "6299", 2, td = tdnew)
tdnew <- update_testdata("ID11871919", "6342", 2, td = tdnew)
tdnew <- update_testdata("ID11880686", "6357", 3, td = tdnew)
tdnew <- update_testdata("ID11880904", "6357", 3, td = tdnew)
tdnew <- update_testdata("ID11884391", "6393", 2, td = tdnew)
tdnew <- update_testdata("ID11884455", "6383", 2, td = tdnew)
tdnew <- update_testdata("ID11886862", "6438", 2, td = tdnew)
tdnew <- update_testdata("ID11890077", "6432", 1:3, td = tdnew)
tdnew <- update_testdata("ID11897192", "6412", 2, td = tdnew)
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DPB1_TEST_1267 <- tdnew2
devtools::use_data(HLA_DPB1_TEST_1267, pkg = "~/Devel/HLAdata", overwrite = TRUE)

## testdata DPB1; basis ID 1412  ##############################################
HLA_DPB1_TEST_1412 <- get_testdata(gene = "DPB1", nextype_basis_id = "1412")
tdnew <- HLA_DPB1_TEST_1412
tdnew <- update_testdata("ID12034763", "7045", 2, td = tdnew)
tdnew <- update_testdata("ID12036341", "7045", 4, td = tdnew)
tdnew <- update_testdata("ID12050950", "7094", 3, td = tdnew)
tdnew <- update_testdata("ID12052242", "7052", 2, td = tdnew)
tdnew <- update_testdata("ID12054548", "7052", 2, td = tdnew)
tdnew <- update_testdata("ID12055199", "7052", 4, td = tdnew)
tdnew <- update_testdata("ID12055367", "7053", 3, td = tdnew)
tdnew <- update_testdata("ID12087059", "7258", 3, td = tdnew)
tdnew <- update_testdata("ID12087274", "7158", 2, td = tdnew)
tdnew <- update_testdata("ID12096400", "7203", 2, td = tdnew)
tdnew <- update_testdata("ID12099145", "7161", 2, td = tdnew)
tdnew <- update_testdata("ID12099544", "7162", 2, td = tdnew)
tdnew <- update_testdata("ID12105078", "7197", 1, td = tdnew)
tdnew <- update_testdata("ID12113816", "7156", 4, td = tdnew)
tdnew <- update_testdata("ID12122625", "7209", 3, td = tdnew)
tdnew <- update_testdata("ID12038322", "7048", 4, td = tdnew)
tdnew <- update_testdata("ID12053576", "7095", 4, td = tdnew)
tdnew <- update_testdata("ID12053820", "7095", 4, td = tdnew)
tdnew <- update_testdata("ID12078848", "7128", 2, td = tdnew)
tdnew <- update_testdata("ID12105078", "7197", 1:2, td = tdnew)
tdnew <- update_testdata("ID12129708", "7217", 2, td = tdnew)
tdnew <- update_testdata("ID12133673", "7271", 2, td = tdnew)
tdnew <- update_testdata("ID12148331", "7343", 4, td = tdnew)
tdnew <- update_testdata("ID12148406", "7302", 2, td = tdnew)
tdnew <- update_testdata("ID12148591", "7343", 4, td = tdnew)
tdnew <- update_testdata("ID12148966", "7343", 2, td = tdnew)
tdnew <- update_testdata("ID12148982", "7343", 2, td = tdnew)
tdnew <- update_testdata("ID12150244", "7344", 3, td = tdnew)
tdnew <- update_testdata("ID12152167", "7351", 4, td = tdnew)
tdnew <- update_testdata("ID12153055", "7350", 3, td = tdnew)
tdnew <- update_testdata("ID12160886", "7377", 2, td = tdnew)
tdnew <- update_testdata("ID12162143", "7377", 4, td = tdnew)
tdnew <- update_testdata("ID12162478", "7579", 4, td = tdnew)
tdnew <- update_testdata("ID12176334", "7424", 3, td = tdnew)
tdnew <- update_testdata("ID12176343", "7424", 4, td = tdnew)
tdnew <- update_testdata("ID12177294", "7419", 4, td = tdnew)
tdnew <- update_testdata("ID12189046", "7454", 2, td = tdnew)
tdnew <- update_testdata("ID12189616", "7509", 4, td = tdnew)
tdnew <- update_testdata("ID12191201", "7455", 4, td = tdnew)
tdnew <- update_testdata("ID12194008", "7500", 4, td = tdnew)
tdnew <- update_testdata("ID12195655", "7500", 4, td = tdnew)
tdnew <- update_testdata("ID12195657", "7500", 4, td = tdnew)
tdnew <- update_testdata("ID12198996", "7580", 4, td = tdnew)
tdnew <- update_testdata("ID12201177", "7496", 2, td = tdnew)
tdnew <- update_testdata("ID12205580", "7505", 3, td = tdnew)
tdnew <- update_testdata("ID12210197", "7584", 4, td = tdnew)
tdnew <- update_testdata("ID12149140", "7343", 4, td = tdnew)
tdnew <- update_testdata("ID12223938", "7703", 2, td = tdnew)
tdnew <- update_testdata("ID12225510", "7697", 2, td = tdnew)
tdnew <- update_testdata("ID12231452", "7678", 2, td = tdnew)
tdnew <- update_testdata("ID12215122", "7597", 2, td = tdnew)
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DPB1_TEST_1412 <- tdnew2
devtools::use_data(HLA_DPB1_TEST_1412, pkg = "~/Devel/HLAdata", overwrite = TRUE)





## testdata A; basis ID 1267  ##############################################
testdata_a_1267 <- get_testdata(gene = "A", nextype_basis_id = "1267")
tdnew <- testdata_a_1267
tdnew2 <- clean_testing_data(td = tdnew)
testdata_a_1267 <- tdnew2
devtools::use_data(testdata_a_1267, overwrite = TRUE)

## testdata A; basis ID 1412  ##############################################
testdata_a_1412 <- get_testdata(gene = "A", nextype_basis_id = "1412")
tdnew <- testdata_a_1412
tdnew2 <- clean_testing_data(td = tdnew)
testdata_a_1412 <- tdnew2
devtools::use_data(testdata_a_1412, overwrite = TRUE)

## testdata B; basis ID 1267  ##############################################
testdata_b_1267 <- get_testdata(gene = "B", nextype_basis_id = "1267")
tdnew <- testdata_b_1267
tdnew2 <- clean_testing_data(td = tdnew)
testdata_b_1267 <- tdnew2
devtools::use_data(testdata_b_1267, overwrite = TRUE)

## testdata B; basis ID 1412  ##############################################
testdata_b_1412 <- get_testdata(gene = "B", nextype_basis_id = "1412")
tdnew <- testdata_b_1412
tdnew2 <- clean_testing_data(td = tdnew)
testdata_b_1412 <- tdnew2
devtools::use_data(testdata_b_1412, overwrite = TRUE)

## testdata C; basis ID 1267  ##############################################
testdata_c_1267 <- get_testdata(gene = "C", nextype_basis_id = "1267")
tdnew <- testdata_c_1267
tdnew2 <- clean_testing_data(td = tdnew)
testdata_c_1267 <- tdnew2
devtools::use_data(testdata_c_1267, overwrite = TRUE)

## testdata B; basis ID 1412  ##############################################
testdata_c_1412 <- get_testdata(gene = "C", nextype_basis_id = "1412")
tdnew <- testdata_c_1412
tdnew2 <- clean_testing_data(td = tdnew)
testdata_c_1412 <- tdnew2
devtools::use_data(testdata_c_1412, overwrite = TRUE)



